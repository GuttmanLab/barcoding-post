import argparse
from contact import Contacts

def main():
    args = parse_arguments()

    contacts = Contacts(build = args.assembly,
                        chromosome = args.chromosome,
                        resolution = args.resolution,
                        downweighting = args.downweighting)

    contacts.get_raw_contacts(clusters_file = args.clusters,
                              min_cluster_size = args.min_cluster_size,
                              max_cluster_size = args.max_cluster_size)

    contacts.write_contacts_to_file(outfile = args.raw_contacts, fmt = "%1f")

    contacts.ice_raw_contacts(raw_contacts_file = args.raw_contacts,
                              bias_file = args.biases,
                              iterations = args.iterations,
                              hicorrector_path = args.hicorrector)

    contacts.write_contacts_to_file(outfile = args.iced, fmt = "%1f")

    contacts.truncate_to_median_diagonal_value()

    contacts.write_contacts_to_file(outfile = args.output, fmt = "%1f")


def parse_arguments():
    parser = argparse.ArgumentParser(
        description = 'Generates an ICEd matrix of intra- or ' +
                'interchromosomal contacts from a cluster file.')
    parser.add_argument('--clusters',
                        metavar = "FILE",
                        action = "store",
                        help = "The input clusters file.")
    parser.add_argument('--raw_contacts',
                        metavar = "FILE",
                        action = "store",
                        help = "An intermediate file of raw (un-ICEd) " + 
                               "contacts.")
    parser.add_argument('--biases',
                        metavar = "FILE",
                        action = "store",
                        help = "An intermediate file of biases generated by " +
                               "hicorrector.")
    parser.add_argument('--iced',
                        metavar = "FILE",
                        action = "store",
                        help = "An intermediate file of ICEd contacts.")
    parser.add_argument('-o', '--output',
                        metavar = "FILE",
                        action = "store",
                        help = "The output ICEd and scaled contacts file.")
    parser.add_argument('--assembly',
                        metavar = "ASSEMBLY",
                        action = "store",
                        choices = ["mm9", "hg19"],
                        default = "mm9",
                        help = "The genome assembly. (default mm9)")
    parser.add_argument('--chromosome',
                        metavar = "CHROM",
                        action = "store",
                        help = "The chromosome of interest. Input 'genome' " +
                               "for an interchromosomal matrix.")
    parser.add_argument('--max_cluster_size',
                        metavar = 'MAX',
                        type = int,
                        action = 'store',
                        default = 1000,
                        help = "Skip read-clusters with more reads than MAX. (default 1000)")
    parser.add_argument('--min_cluster_size',
                        metavar = 'MIN',
                        type = int,
                        action = 'store',
                        default = 2,
                        help = "Skip read-clusters with fewer reads than MIN. (default 2)")
    parser.add_argument('--resolution',
                        metavar = 'INT',
                        type = int,
                        action = 'store',
                        default = 1000000,
                        help = "Contact matrix resolution in bp. (default 1000000 = 1 Mb)")
    parser.add_argument('--downweighting',
                        metavar = 'DW',
                        action = 'store',
                        default = "none",
                        choices = ["none", "n_minus_one"],
                        help = "Downweighting strategy")
    parser.add_argument('--hicorrector',
                        metavar = "FILE",
                        action = 'store',
                        default = "/storage/Software/hicorrector/1.2/bin/ic",
                        help = "Path to hicorrector ic. (default /storage/Software/hicorrector/1.2/bin/ic)")
    parser.add_argument('--iterations',
                        metavar = 'INT',
                        type = int,
                        action = 'store',
                        default = 100,
                        help = "Number of ICE iterations (default 100)")
    return parser.parse_args()


if __name__ == "__main__":
    main()
