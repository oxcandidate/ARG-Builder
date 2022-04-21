import numpy as np
import sys
import getopt


remove_multiallelic = True
keep_seq_labels = True
# Reference should be first sequence

verbose = False

bases = ['A', 'G', 'C', 'T']


def read_problematic_sites(maskfile):
    sites_to_mask = []

    file = open(maskfile, "r")
    for line in file.readlines():
        if line[0] != '#':
            sites_to_mask.append(int(line.split('\t')[1]))
    
    return sites_to_mask

def read_sequences(inputfile):
    sequences = []
    sequence_labels = []

    file = open(inputfile, "r")
    seq = ""
    for line in file.readlines():
        if line[0] == '>':
            sequence_labels.append(line[1:].replace('\n', ''))
            if(seq != ""):
                sequences.append(list(seq.upper()))
                seq = ""
        else:
            seq += line.replace('\n', '')
    sequences.append(list(seq.upper()))
    file.close()

    num_sequences = len(sequences)
    num_cols = len(sequences[0])

    print("seqs: ", num_sequences)
    print("cols: ", num_cols)

    return (sequence_labels, sequences)

def mask_sites(sequences, sites_to_mask):
    ref_seq = sequences[0]
    seq_length = len(ref_seq)

    multiallelic_sites = []
    
    for site in range(seq_length):
        bases_seen = []

        for seq in sequences:
            if seq[site] not in bases_seen:
                bases_seen.append(seq[site])
        
        if 'N' in bases_seen:
            bases_seen.remove('N')

        mask_site = False
        if len(bases_seen) > 2:
            if verbose:
                print(site, " is multiallelic")
            mask_site = True
            multiallelic_sites.append(site)
        elif site in sites_to_mask:
            mask_site = True
        
        if mask_site:
            for seq in sequences:
                seq[site] = 'N'
    
    return (sequences, multiallelic_sites)

def make_binary(sequences):
    ref_seq = sequences[0]
    seq_length = len(ref_seq)
    num_sequences = len(sequences)
    binary_sequences = [[0] * seq_length for seq in sequences]

    for site in range(seq_length):
        ref_base = ref_seq[site]
        
        if ref_base != 'N':
            for i in range(num_sequences):
                if sequences[i][site] in ['N', ref_base]:
                    binary_sequences[i][site] = 0
                else:
                    binary_sequences[i][site] = 1
    
    return binary_sequences


def remove_zero_columns(binary_sequences):
    seq_length = len(binary_sequences[0])
    num_sequences = len(binary_sequences)

    keep_cols = []

    for site in range(seq_length):
        for seq in binary_sequences:
            if seq[site] == 1:
                keep_cols.append(site)
                break
    
    new_sequences = [[seq[site] for site in keep_cols] for seq in binary_sequences]
    return new_sequences

def output(binary_sequences, sequence_labels, outfile):
    file = open(outfile, 'w')
    for j in range(len(binary_sequences)):
        if keep_seq_labels:
            file.write(">" + sequence_labels[j] + "\n")
        line = ""
        for v in binary_sequences[j]:
            line += str(v)
        file.write(line + "\n")

    file.close()


def main(argv):
    inputfile = ''
    outputfile = ''
    problematicsitesfile = ''

    try:
        opts, args = getopt.getopt(argv, "hm:i:o:", ["mfile=", "ifile=", "ofile="])
    except getopt.GetoptError:
        print('process_fasta_to_binary.py -m <problematic sites file (a vcf)> -i <input fasta file> -o <outputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('process_fasta_to_binary.py -m <problematic sites file (a vcf)> -i <input fasta file> -o <outputfile>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
        elif opt in ("-m", "--mfile"):
            problematicsitesfile = arg

    print('Input file is: ', inputfile)
    print('Output file is: ', outputfile)
    print('mask file is: ', problematicsitesfile)

    sites_to_mask = read_problematic_sites(problematicsitesfile)
    (sequence_labels, sequences) = read_sequences(inputfile)
    mask_sites(sequences, sites_to_mask)
    bin_seqs = make_binary(sequences)
    # bin_seqs = remove_zero_columns(bin_seqs)
    output(bin_seqs, sequence_labels, outputfile)



if __name__ == "__main__":
    main(sys.argv[1:])