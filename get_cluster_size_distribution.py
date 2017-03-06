from collections import Counter
import sys

def main():
    counts = get_cluster_sizes(sys.argv[1])
    print_cluster_sizes(counts, sys.argv[2])

def get_cluster_sizes(in_file):
    counts = Counter()

    with open(in_file, 'r') as f:
        for line in f:
            fields = line.split()
            counts[len(fields) - 1] += (len(fields) - 1)

    return counts


def print_cluster_sizes(counts, out_file):
    with open(out_file, 'w') as f:
        for size, count in counts.iteritems():
            if size == 1:
                category = "SINGLE"
            elif 2 <= size <= 10:
                category = "TWO_TO_TEN"
            elif 11 <= size <= 100:
                category = "ELEVEN_TO_ONE_HUNDRED"
            elif 101 <= size <= 1000:
                category = "ONE_HUNDRED_ONE_TO_ONE_THOUSAND"
            else:
                category = "OVER_ONE_THOUSAND"
            f.write(str(size) + "\t" + str(count) + "\t" + category + "\n")

if __name__ == "__main__":
    main()
