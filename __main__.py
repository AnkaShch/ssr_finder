import argparse
import logging
import multiprocessing as mp
import pyfaidx
import sys

from finder import __script_name__, __version__, fasta_motif_scan


def main():
    #############
    # arg parse #
    #############
    parser = argparse.ArgumentParser(prog=__script_name__, epilog="%s v%s" % (__script_name__, __version__))

    # FASTA index will be created if it does not exist when pyfaidx Fasta is initialized
    parser.add_argument('fastaFile', help="Path to FASTA file.")
    parser.add_argument('--motif', '-m', help="A degenerate sequence motif. Can be specified multiple times.",
                        action='append')
    parser.add_argument('--motif_file',
                        help="A file containing motifs to search file. Can be specified multiple times.",
                        action='append')
    parser.add_argument('--output_prefix', '-o', help="Prefix of generated bed files. Default is detected_ssrs",
                        default="detected_ssrs")
    parser.add_argument('--region', '-r',
                        help="Defines region to search for sites. Use 'contig:start-end' for regions, or 'contig' for whole contig. If no regions are specified, the entire FASTA file will be searched! 0-based format, start is included, end in not included",
                        action='append', required=False)
    parser.add_argument('--cores', '-p', help="Run search on multiple contigs / strands simultaneously", type=int,
                        default=1)
    parser.add_argument('--valid_regex',
                        help="Query is valid regex. *WILL NOT* reverse complement. Specify sequence and strand with careful consideration.",
                        default=False, action='store_true')
    parser.add_argument('--distance', '-d',
                        help="The maximum allowable distance between motifs in bp, within which it is assumed that it all falls into the repeat region. Default is 0.",
                        type=int, default=0)
    parser.add_argument('--number', '-n',
                        help="The minimum number of motives in one SSR. Default is 2.",
                        type=int, default=2)
    # molecule
    parser.add_argument('--strand', '-s', help="Default searches both strands, but can be set to only one.",
                        choices=['+', '-', 'both'], default='both')
    parser.add_argument('--motif_overlaps', help="If this option is set it turns on overlapping motif matches.",
                        action='store_true')
    parser.add_argument("--loglevel", choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'], default='INFO')

    args = parser.parse_args()

    #################
    # setup logging #
    #################
    logging.basicConfig(format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    logger = logging.getLogger(__script_name__)
    logger.setLevel(args.loglevel)

    ########################
    # fix up input parsing #
    ########################
    allow_overlapping_motifs = args.motif_overlaps

    if args.strand == "both":
        run_strands = ('+', '-')
    elif args.strand == '-':
        run_strands = ('-',)
    else:
        run_strands = ('+',)

    if args.motif is None and args.motif_file is None:
        logger.critical("Must specify at least one --motif or path to a --motif_file")
        sys.exit(1)

    if args.motif is None:
        args.motif = []

    if args.motif_file is not None:
        for path in args.motif_file:
            with open(path, 'r') as MOTIF_IN:
                for line in MOTIF_IN:
                    line = line.strip()

                    if len(line) == 0:
                        continue
                    elif line[0] == '#':
                        continue
                    args.motif.append(line.upper())

    args.motif = list(set([seq.upper() for seq in args.motif]))  # uniquify motifs in case they show up multiple times

    if len(args.motif) == 0:
        logger.critical("No motifs specified!")
        sys.exit(1)

    ######################
    # set options string #
    ######################
    options = "\n%s v%s\n\nOptions\n=======\n" % (__script_name__, __version__)
    options += "FASTA: %s\n" % (args.fastaFile)
    options += "Motifs to search: %s\n" % (str(args.motif))
    options += "Strands of FASTA to search: %s\n" % (str(run_strands))
    options += "Max distance between motifs: %s\n" % (str(args.distance))
    options += "Min number of motifs in SSR: %s\n" % (str(args.number))
    options += "Output files: %s.bed %s_full.bed\n" % (args.output_prefix, args.output_prefix)
    options += "Regions: %s\n" % ("All contigs" if args.region is None else str(args.region))
    options += "Motifs are valid regex (no processing): %s\n" % (str(args.valid_regex))
    options += "Allow overlapping motifs?: %s\n" % (str(allow_overlapping_motifs))
    options += "Processes: %s\n" % (args.cores)  # cores
    options += "Log level: %s\n" % (str(args.loglevel))

    logger.info(options)

    ###################
    # Open FASTA file #
    ###################
    region_list = []

    with pyfaidx.Fasta(args.fastaFile, as_raw=True) as FAIDX:
        if args.region is None:
            args.region = FAIDX.keys()

        for reg in args.region:
            contig, start, end = pyfaidx.ucsc_split(reg)  # returns [0,1) coordinates with NoneType if not start or end

            # no need for checking for contig in FASTA file, as this is done by pyfaidx
            if start is None:
                start = 0  # pyfaidx using 0-base start. Remember for custom regions!!!
            if end is None:
                end = len(FAIDX[contig])

            for curr_motif in args.motif:
                for curr_strand in run_strands:
                    region_list.append((curr_motif, contig, start, end, curr_strand, args.distance, args.number))

                    logger.debug("%s\n" % (str(region_list[-1])))

    ################
    # process data #
    #############################
    # async process the contigs #
    #############################
    work_pool = mp.Pool(processes=args.cores)
    working_data = [work_pool.apply_async(fasta_motif_scan, args=(
        args.fastaFile, x, args.valid_regex, allow_overlapping_motifs)) for x in region_list]
    work_output = [x.get() for x in working_data]

    ##############################
    # figure order out           #
    # b/c executed asyncronously #
    ##############################
    site_count = 0

    result_order = []

    for reg in region_list:
        index = 0

        for out, dict, current_site_count in work_output:
            if out == reg:
                logger.debug("Returned result %s matches %s" % (index, ' '.join([str(x) for x in reg])))

                result_order.append(index)
                site_count += current_site_count

                break
            else:
                index += 1

    if len(result_order) != len(region_list):
        logger.critical("Did not find results for all requested regions")
        sys.exit(1)

    #################
    # write outputs #
    #################
    output_file = args.output_prefix + ".bed"
    output_file_full = args.output_prefix + "_full.bed"

    out = open(output_file, 'w')
    out_full = open(output_file_full, 'w')

    out.write(
        "#Contig\tStart\tEnd\tMotif\tLength\tStrand\tNumber of motifs\tSequence\n")
    out_full.write(
        "#Contig\tStart\tEnd\tMotif\tLength\tStrand\tNumber of motifs\tNumber of inserts\tLength of inserts\tSequence\n")

    for idx in result_order:
        for (repeat_region, seq) in work_output[idx][1]:
            out.write("%s\t%s\n" % (repeat_region.to_string(), seq))
            out_full.write("%s\t%s\n" % (repeat_region.to_string_full(), seq))

    out.close()
    out_full.close()

    logger.info("%s total sites found" % site_count)


if __name__ == "__main__":
    main()
