import argparse
import cluster as c

def main():
    args = parse_arguments()
    clusters = c.get_clusters(args.input, args.num_barcodes)
    c.write_clusters_to_file(clusters, args.output)

def parse_arguments():
    parser = argparse.ArgumentParser(
        description = 'Generates a clusters file from a BAM file.')
    parser.add_argument('-i', '--input',
                        metavar = "FILE",
                        action = "store",
                        help = "The input BAM file.")
    parser.add_argument('-o', '--output',
                        metavar = "FILE",
                        action = "store",
                        help = "The output clusters file.")
    parser.add_argument('-n', '--num_barcodes',
                        metavar = 'INT',
                        type = int,
                        action = 'store',
                        help = "The number of barcodes contained in the name " +
                               "of each BAM record.")
    return parser.parse_args()

if __name__ == "__main__":
    main()
