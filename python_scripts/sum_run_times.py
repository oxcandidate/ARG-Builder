import sys
import getopt
from statistics import median


# Reads from csv containing many runs and create a grid for optimal args

def sum_runtimes(inputfile):
    total_time = 0.0
    times = []

    file = open(inputfile, "r")
    for line in file.readlines():
        if len(line) == 0:
            continue
        elif line[0] == 'r':
            continue  # this is the header which starts recombination

        values = line.split(',')
        run_time = float(values[4])
        times.append(run_time)
        total_time += run_time

    file.close()

    print("total: ", total_time, " max: ", max(times), " min: ", min(times), " median: ", median(times), " mean: ", sum(times) / len(times))
    return times


def main(argv):
    inputfile = ''

    try:
        opts, args = getopt.getopt(argv, "hi:", ["ifile="])
    except getopt.GetoptError:
        print('sum_run_times.py -i <input file>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('sum_run_times.py -i <input file>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg

    print('Input file is: ', inputfile)

    sum_runtimes(inputfile)


if __name__ == "__main__":
    main(sys.argv[1:])
