from chromsizes import ChromosomeSizes
from itertools import combinations
import numpy as np
import subprocess

class Contacts:
    def __init__(self, chromosome, assembly = "mm9", resolution = 1000000):

        self._chromosome = chromosome
        self._resolution = resolution
        self._assembly = assembly


    def get_raw_contacts(self, clusters_file, min_cluster_size = 2,
                max_cluster_size = 1000):

        self.init_chromosome_matrix()

        with open(clusters_file, 'r') as f:

            for line in f:
                reads = line.split()[1:]
                if not min_cluster_size <= len(reads) <= max_cluster_size:
                    continue
                bins = set()

                for read in reads:
                    chrom, position = read.split(':')
                    if chrom == self._chromosome:
                        read_bin = int(position) // self._resolution
                        bins.add(read_bin)

                for bin1, bin2 in combinations(bins, 2):
                    self._contacts[bin1][bin2] += 1
                    self._contacts[bin2][bin1] += 1


    def init_chromosome_matrix(self):
        chromosome_size = ChromosomeSizes.assembly[self._assembly][self._chromosome]
        num_bins = -(-chromosome_size // self._resolution)
        self._contacts = np.zeros((num_bins, num_bins))


    def write_contacts_to_file(self, outfile, fmt):
        np.savetxt(outfile, self._contacts, delimiter = "\t", fmt = fmt)


    def ice_raw_contacts(self, raw_contacts_file, bias_file, iterations,
                hicorrector_path):

        biases = self.calculate_bias_factors(raw_contacts_file = raw_contacts_file,
                bias_file = bias_file, hicorrector = hicorrector_path,
                iterations = iterations)

        median_diagonal_value = self.get_median_diagonal_value()

        for row in xrange(self._contacts.shape[0]):
            for col in xrange(self._contacts.shape[1]):
                val = float(self._contacts[row][col])
                if val > 0:
                    val /= (biases[row] * biases[col])
                    val = (1 if val >= median_diagonal_value
                           else val / median_diagonal_value)
                    self._contacts[row][col] = val
                
    
    def calculate_bias_factors(self, raw_contacts_file, bias_file, hicorrector,
                iterations):
        skip_first_row = "0"    # 0 == don't skip
        skip_first_column = "0"
        num_lines = self._contacts.shape[0]
        subprocess.check_call([hicorrector, raw_contacts_file, str(num_lines),
                str(iterations), skip_first_row, skip_first_column, bias_file])
        return self.parse_bias_file(bias_file)


    def parse_bias_file(self, bias_file):
        biases = []
        with open(bias_file) as f:
            for line in f:
                biases.append(float(line.strip()))
        return biases


    def get_median_diagonal_value(self):
        diagonal_values = []
        self._contacts

        for i in xrange(self._contacts.shape[0] - 1):
            diagonal_values.append(self._contacts[i + 1][i])
            diagonal_values.append(self._contacts[i][i + 1])

        return np.median(diagonal_values)
