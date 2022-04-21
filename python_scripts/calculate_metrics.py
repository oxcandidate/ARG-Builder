
from audioop import avg
from statistics import mean, median
import sys
import getopt
from typing import Sequence
import numpy as np

from clean_binary_sequences import read_sequences

verbose = False

def pairwise_hamming_old(sequences):
    num_sequences = len(sequences)
    num_cols = len(sequences[0])

    distances = np.zeros((num_sequences, num_sequences), dtype=int)
    diffs_data = []

    for i in range(num_sequences):
        if i % 30 == 0 and verbose:
            print(i)
        for j in range(i+1, num_sequences):
            diffs = sum([1 for c in range(num_cols) if sequences[i][c] != sequences[j][c]])
            # for c in range(num_cols):
            #     if sequences[i][c] != sequences[j][c]:
            #         diffs += 1
            distances[i][j] = diffs
            distances[j][i] = diffs
            diffs_data.append(diffs)
    if verbose:
        print()
    print("hamming distance")
    print("avg", mean(diffs_data), " max", max(diffs_data), " median", median(diffs_data))
    return (mean(diffs_data), max(diffs_data), median(diffs_data))


def pairwise_hamming(sequences):
    num_sequences = len(sequences)
    num_cols = len(sequences[0])

    seq_descriptions = [[site for site in range(num_cols) if seq[site] == 1] for seq in sequences]

    distances = np.zeros((num_sequences, num_sequences), dtype=int)
    diffs_data = []

    for i in range(num_sequences):
        if i % 30 == 0 and verbose:
            print(i)
        for j in range(i+1, num_sequences):
            seqA = seq_descriptions[i]
            seqB = seq_descriptions[j]
            sames = 0
            a = 0
            b = 0
            while(a < len(seqA) and b < len(seqB)):
                if seqA[a] == seqB[b]:
                    a += 1
                    b += 1
                    sames += 1
                elif seqA[a] < seqB[b]:
                    a += 1
                else:
                    b += 1
            
            diffs = len(seqA) + len(seqB) - 2*sames

            distances[i][j] = diffs
            distances[j][i] = diffs
            diffs_data.append(diffs)
    if verbose:
        print()
    print("hamming distance")
    print("avg", mean(diffs_data), " max", max(diffs_data), " median", median(diffs_data))
    return (mean(diffs_data), max(diffs_data), median(diffs_data))

def mutation_sums(sequences):
    print("mutations per sequence")
    muts = [sum(seq) for seq in sequences]
    print("avg", mean(muts), " max", max(muts), " median", median(muts))

    print("mutations per site")
    site_muts = [ sum([seq[site] for seq in sequences]) for site in range(len(sequences[0]))]
    print("avg", mean(site_muts), " max", max(site_muts), " median", median(site_muts))

    return (mean(muts), max(muts), median(muts), mean(site_muts), max(site_muts), median(site_muts))

def output(outfile, rowname, values):
    file = open(outfile, 'a')

    file.write(rowname + "," + ",".join([str(v) for v in values]) + "\n")

    file.close()

def main(argv):
    inputfile = ''
    outputfile = ''

    try:
        opts, args = getopt.getopt(argv, "hi:o:", ["ifile=", "ofile="])
    except getopt.GetoptError:
        print('clean_binary_sequences.py -i <input file> -o <output file>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('clean_binary_sequences.py -i <input file> -o <output file>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg

    print('Input file is: ', inputfile)

    (_, sequences, _) = read_sequences(inputfile)

    if len(sequences) == 0 or len(sequences[0]) == 0:
        print("no data to give metrics for")
        return

    (h_mean, h_max, h_median) = pairwise_hamming(sequences)
    (seq_mean, seq_max, seq_median, site_mean, site_max, site_median) = mutation_sums(sequences)

    output(outputfile, inputfile.split(".")[0], [h_mean, h_max, h_median, seq_mean, seq_max, seq_median, site_mean, site_max, site_median])

    



if __name__ == "__main__":
    main(sys.argv[1:])