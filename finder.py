#!/usr/bin/env python

"""
motif_scraper: a tool for searching for nucleic acid motifs using degenerate nucleotide queries
"""

##########
# import #
##########
import multiprocessing as mp
import pyfaidx
import regex
import six
import sys
import tempfile
import time

####################
# Version and name #
####################
__script_path__ = sys.argv[0]
__script_name__ = __script_path__.split('/')[-1].split('\\')[-1]
__version__ = '1.0'


########
# fxns #
########	
def rev_comp(seq, molecule='dna'):
    """ DNA|RNA seq -> reverse complement
	"""

    if molecule == 'dna':
        nuc_dict = {"A": "T", "B": "V", "C": "G", "D": "H", "G": "C", "H": "D", "K": "M", "M": "K", "N": "N", "R": "Y",
                    "S": "S", "T": "A", "V": "B", "W": "W", "Y": "R"}
    elif molecule == 'rna':
        nuc_dict = {"A": "U", "B": "V", "C": "G", "D": "H", "G": "C", "H": "D", "K": "M", "M": "K", "N": "N", "R": "Y",
                    "S": "S", "U": "A", "V": "B", "W": "W", "Y": "R"}
    else:
        raise ValueError("rev_comp requires molecule to be dna or rna")

    if not isinstance(seq, six.string_types):
        raise TypeError("seq must be a string!")

    return ''.join([nuc_dict[c] for c in seq.upper()[::-1]])


def make_degenerate_regex(motif_seq, molecule='dna'):
    """ Degenerate sequence -> regex
	Example: NNYCGAARN -> [ACGT]{2}[CT]CGA{2}[AG][ACGT]
	"""

    if not isinstance(motif_seq, six.string_types):
        raise TypeError("motif_seq must be a string!")

    if molecule == 'dna':
        degenerate_code = {"A": "A", "B": "[CGT]", "C": "C", "D": "[AGT]", "G": "G", "H": "[ACT]", "K": "[GT]",
                           "M": "[AC]", "N": "[ACGT]", "R": "[AG]", "S": "[GC]", "T": "T", "V": "[ACG]", "W": "[AT]",
                           "Y": "[CT]"}
    elif molecule == 'rna':
        degenerate_code = {"A": "A", "B": "[CGU]", "C": "C", "D": "[AGU]", "G": "G", "H": "[ACU]", "K": "[GU]",
                           "M": "[AC]", "N": "[ACGU]", "R": "[AG]", "S": "[GC]", "U": "U", "V": "[ACG]", "W": "[AU]",
                           "Y": "[CU]"}
    else:
        raise ValueError("make_degenerate_regex requires molecule to be dna or rna")

    regex_string = ''

    idx = 0

    while idx < len(motif_seq):
        curr = motif_seq[idx]

        count = 1

        for next_idx in range(idx + 1, len(motif_seq)):
            next = motif_seq[next_idx]

            if next == curr:
                count += 1
            else:
                break

        regex_string += degenerate_code[curr]

        if count > 1:
            idx = idx + count - 1
            regex_string += "{%s}" % (count)

        idx += 1

    return regex_string

###########
# classes #
###########
class RepeatRegionStat:
    def __init__(self, motif, contig, positionStart, strand, molecule='dna'):
        self.contig = contig
        self.start = positionStart
        self.end = self.start
        self.motif = motif
        self.strand = strand
        self.motif_numb = 0
        self.insert_numb = 0
        self.insert_length = 0
        self.molecule = molecule

    def add(self, regexMatch, start):
        motif_length = len(regexMatch.group())
        motif_start = start + regexMatch.start()
        motif_end = start + regexMatch.end()
        insert_length = motif_start - self.end

        # if it is first motif
        if self.motif_numb == 0:
            self.start = motif_start
            self.end += insert_length + motif_length
            self.motif_numb += 1
        else:
            self.end += insert_length + motif_length
            self.motif_numb += 1
            if insert_length > 0:
                self.insert_numb += 1
                self.insert_length += insert_length

    def to_string(self):
        return "%s\t%s\t%s\t%s\t%s\t%s\t%s" % (
            self.contig, self.start, self.end, self.motif, abs(self.end - self.start), self.strand, self.motif_numb)

    def to_string_full(self):
        return "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (
            self.contig, self.start, self.end, self.motif, abs(self.end - self.start), self.strand, self.motif_numb,
            self.insert_numb, self.insert_length)


def get_seq(repeat_region, sequence, start):
    if repeat_region.strand == '+':
        return sequence[repeat_region.start - start:repeat_region.end - start]
    elif repeat_region.strand == '-':
        rev_seq = rev_comp(sequence[repeat_region.start - start:repeat_region.end - start], repeat_region.molecule)
        return rev_seq
    else:
        raise ValueError("SequenceMotif classes require strand to be either '+' or '-'")


######################
# fxns using classes #
######################
def fasta_motif_scan(fasta_fname, input_tuples, regex_ready=False, allow_overlaps=True,
                     molecule='dna'):
    """
	fasta_fname = string path to FASTA file
	input_tuples = tuple containing (1) motif sequence, (2) contig name, (3) start position*, (4) end position,
		(5) strand to search, (6) maximum distance between motif, (7) minimum number of motifs in repeat

	*start is expected to be 0-base, end is expected to be 1-base
	"""

    ###################
    # validity checks #
    ###################
    global motif_regex
    if not isinstance(fasta_fname, six.string_types):
        raise TypeError("In fasta_motif_scan, fasta_fname must be a string!")
    elif not isinstance(molecule, six.string_types):
        raise TypeError("In fasta_motif_scan, molecule must be a string!")
    elif isinstance(input_tuples, six.string_types) or not isinstance(input_tuples, tuple):
        raise (TypeError("In fasta_motif_scan, input_tuples should be a tuple!"))
    elif not type(regex_ready) is bool:
        raise (TypeError("In fasta_motif_scan, regex_ready should be a bool!"))

    #########################
    # setup some local vars #
    #########################
    motif_seq, contig, start, end, strand, distance, number = input_tuples

    if regex_ready == False:
        if strand == '+':
            motif_regex = make_degenerate_regex(motif_seq)
        elif strand == '-':
            motif_regex = make_degenerate_regex(rev_comp(motif_seq, molecule))
    else:
        if strand == '+':
            motif_regex = motif_seq
        elif strand == '-':
            motif_regex = rev_comp(motif_seq, molecule)

    regex_compiled = regex.compile(motif_regex)

    #######
    # run #
    #######
    site_count = 0

    site_list = []

    with pyfaidx.Fasta(fasta_fname, as_raw=True) as FAIDX:
        sequence = str(FAIDX[contig][start:end]).upper()

        # self, motif, contig, positionStart, strand, motifMatches molecule='dna'
        repeat_region = RepeatRegionStat(motif_seq, contig, start, strand, molecule)
        for m in regex_compiled.finditer(sequence, overlapped=allow_overlaps):
            if start + m.start() - repeat_region.end <= distance:
                repeat_region.add(m, start)
            else:
                # want only more then 1 motif in repeat
                if repeat_region.motif_numb >= number:
                    site_count += 1
                    site = get_seq(repeat_region, sequence, start)
                    site_list.append((repeat_region, site))

                repeat_region = RepeatRegionStat(motif_seq, contig, start, strand, molecule)
                repeat_region.add(m, start)

    if repeat_region.motif_numb >= number:
        site_count += 1
        site = get_seq(repeat_region, sequence, start)
        site_list.append((repeat_region, site))

    return input_tuples, site_list, site_count


########
# main #
########
if __name__ == "__main__":
    pass
