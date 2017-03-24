from enum import Enum
from itertools import combinations
from random import random
import assembly
import numpy as np
import subprocess

class Downweighting(Enum):
    NONE = 1
    N_MINUS_ONE = 2
    N_OVER_TWO = 3
    UNKNOWN = 4


class Contacts:
    def __init__(self, chromosome, build = "mm9", resolution = 1000000,
                 downweighting = "none"):

        self._chromosome = chromosome
        self._resolution = resolution
        self._assembly = assembly.build(build)
        self._assembly.init_offsets(self._resolution)

        if downweighting == "none":
            self._downweighting = Downweighting.NONE
        elif downweighting == "n_minus_one":
            self._downweighting = Downweighting.N_MINUS_ONE
        elif downweighting == "n_over_two":
            self._downweighting = Downweighting.N_OVER_TWO
        else:
            self._downweighting = Downweighting.UNKNOWN

        if self._chromosome == "genome":
            self.init_genome_matrix()
        else:
            self.init_chromosome_matrix()


    def get_raw_contacts_from_clusters_file(self, clusters_file,
                min_cluster_size = 2, max_cluster_size = 1000):

        if self._chromosome.startswith("chr"):
            self.get_raw_contacts_over_chromosome(clusters_file,
                    min_cluster_size, max_cluster_size)
        elif self._chromosome == "genome":
            self.get_raw_contacts_over_genome(clusters_file,
                    min_cluster_size, max_cluster_size)
        else:
            raise Exception, ("Chromosome ID must start with 'chr' or " + 
                              "equal 'genome'")


    def get_raw_contacts_over_genome(self, clusters_file,
            min_cluster_size, max_cluster_size):

        with open(clusters_file, 'r') as f:

            for line in f:
                reads = line.split()[1:]
                if not min_cluster_size <= len(reads) <= max_cluster_size:
                    continue
                bins = set()

                for read in reads:
                    chrom, position = read.split(':')
                    read_bin = int(position) // self._resolution
                    offset = self._assembly.get_offset(chrom)
                    if offset is not None:  # offset == None if chrom not in dict
                        bins.add(read_bin + offset)

                self.add_bins_to_contacts(bins)


    def get_raw_contacts_over_chromosome(self, clusters_file,
                min_cluster_size, max_cluster_size):

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

                self.add_bins_to_contacts(bins)                


    def get_raw_contacts_from_hic_file(self, hic_file):

        with open(hic_file, 'r') as f:
            for line in f:
                line = line.rstrip()
                pos1, pos2, count = line.split()
                pos1 = int(pos1) // self._resolution
                pos2 = int(pos2) // self._resolution
                count = int(float(count))
                self._contacts[pos1][pos2] = count
                self._contacts[pos2][pos1] = count


    def add_bins_to_contacts(self, bins):

        if len(bins) > 1:
            if self._downweighting == Downweighting.N_OVER_TWO:
                inc = 2.0 / len(bins)
            elif self._downweighting == Downweighting.N_MINUS_ONE:
                inc = 1.0 / (len(bins) - 1)
            else:
                assert self._downweighting == Downweighting.NONE
                inc = 1.0

            for bin1, bin2 in combinations(bins, 2):
                self._contacts[bin1][bin2] += inc
                self._contacts[bin2][bin1] += inc


    def zero_diagonal_entries(self):
        for i in xrange(len(self._contacts)):
            self._contacts[i][i] = 0


    def init_chromosome_matrix(self):
        chromosome_size = self._assembly.get_size(self._chromosome)
        num_bins = -(-chromosome_size // self._resolution)
        self._contacts = np.zeros((num_bins, num_bins))


    def init_genome_matrix(self):
        num_bins = 0
        for chromosome_size in self._assembly._chromsizes.itervalues():
            num_bins += -(-chromosome_size // self._resolution)
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
                val = self._contacts[row][col]
                if val > 0:
                    val /= (biases[row] * biases[col])
                    self._contacts[row][col] = val


    def truncate_to_median_diagonal_value(self):
        median_diagonal_value = self.get_median_diagonal_value()

        for row in xrange(self._contacts.shape[0]):
            for col in xrange(self._contacts.shape[1]):
                val = self._contacts[row][col]
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

        for i in xrange(self._contacts.shape[0] - 1):
            diagonal_values.append(self._contacts[i + 1][i])
            diagonal_values.append(self._contacts[i][i + 1])

        return np.median(diagonal_values)


    def downsample(self, target):
        dim = len(self._contacts)

        total_contacts = 0

        # Only sum contacts from diagonal and upper-triangle
        # Otherwise, will double-count.
        for i in xrange(dim):
            for j in xrange(i + 1):
                total_contacts += self._contacts[i][j]

        downsample_ratio = float(target) / total_contacts

        for i in xrange(dim):
            for j in xrange(i + 1):
                num_contacts = self._contacts[i][j]
                for contact in xrange(num_contacts):
                    if random() < downsample_ratio:
                        num_contacts -= 1
                self._contacts[i][j] = num_contacts
                self._contacts[j][i] = num_contacts


