import sys
import getopt
from statistics import median


# Reads from csv containing many runs and create a grid for optimal args

def split(inputfile, num_sets, is_growing):
    file_root = inputfile.split(".")[0]

    files = [open(file_root + "_" + str(k+1) + ".csv", "w") for k in range(num_sets)]

    file = open(inputfile, "r")
    i = 0
    for line in file.readlines():
        if len(line) == 0:
            continue
        elif line[0] == 'r':
            # this is the header which starts recombination
            for f in files:
                f.write(line)
            continue
        
        if is_growing:
            for j in range(i, num_sets):
                files[j].write(line)
        else:
            files[i].write(line)

        i = (i + 1) % num_sets


    file.close()

    for f in files:
        f.close()


def main(argv):
    inputfile = ''
    num_sets = 0
    is_growing = False

    try:
        opts, args = getopt.getopt(argv, "hi:k:g", ["ifile="])
    except getopt.GetoptError:
        print('split_records.py -i <input file> -k <how many subsets to split into> -g [make subsets growing]')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('split_records.py -i <input file> -k <how many subsets to split into> -g [make subsets growing]')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-k"):
            num_sets = int(arg)
        elif opt in ("-g"):
            is_growing = True

    print('Input file is: ', inputfile, " num_sets: ", num_sets, " is growing: ", is_growing)

    split(inputfile, num_sets, is_growing)


if __name__ == "__main__":
    main(sys.argv[1:])
