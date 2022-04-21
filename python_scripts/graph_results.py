from cmath import inf
import re
import numpy as np
import sys
import getopt
import matplotlib as mpl
from matplotlib import pyplot as plt
from typing import List, Dict, Tuple


np.set_printoptions(edgeitems=30, linewidth=100000,
                    formatter=dict(float=lambda x: "%.3g" % x))

plt.rc('font', size=12) 

# Reads from csv containing many runs and create a grid for optimal args


def read_args(inputfile):
    arg_runs = []

    file = open(inputfile, "r")
    for line in file.readlines():
        if len(line) == 0:
            continue
        elif line[0] == 'r':
            continue  # this is the header which starts recombination

        values = line.split(',')
        recombinations = int(values[0])
        back_mutations = int(values[1]) # On kwarg this is sequential errors
        recurrent_mutations = int(values[2])
        if (recombinations >= 0 and back_mutations >= 0 and recurrent_mutations >= 0):
            # Negative values are used for runs which finished in error
            arg_runs.append((recombinations, back_mutations, recurrent_mutations))

    file.close()

    return arg_runs

def read_lower_bounds(inputfile):
    lower_bounds = {}

    file = open(inputfile, "r")
    for line in file.readlines():
        if len(line) == 0:
            continue
        elif line[0] == 'B':
            continue  # this is the header ("Bound Type..")

        bound_type, bound_value, bound_parameter = line.strip('\n').split(',')

        if bound_type in ["hudson_overlap", "hudson_nonoverlap", "RecMin", "RecMinRobust"]:
            lower_bounds[bound_type] = int(bound_value)
        else:
            if bound_type in lower_bounds:
                lower_bounds[bound_type].append((int(bound_parameter), int(bound_value)))
            else:
                lower_bounds[bound_type] = [(int(bound_parameter), int(bound_value))]

    file.close()

    return lower_bounds


def runs_to_dict(arg_runs):
    fewest_rms = {}

    max_recombs = 0
    max_bms = 0

    for (recombs, bms, rms) in arg_runs:
        max_recombs = max(max_recombs, recombs)
        max_bms = max(max_bms, bms)

        if (recombs, bms) in fewest_rms:
            current_rms = fewest_rms[(recombs, bms)]
            if (rms < current_rms):
                fewest_rms[(recombs, bms)] = rms
        else:
            fewest_rms[(recombs, bms)] = rms

    return (fewest_rms, max_recombs, max_bms)


def calculate_optimal_grid(fewest_rms, max_recombs, max_bms):
    grid = np.full((max_recombs+1, max_bms+1), sys.maxsize, dtype=int)

    # Fill out grid with known args first
    for ((recombs, bms), rms) in fewest_rms.items():
        grid[recombs, bms] = rms

    # Now fill in grid
    for i in range(max_recombs+1):
        for j in range(max_bms+1):
            if i > 0 and j > 0:
                grid[i, j] = min(grid[i, j], grid[i-1, j], grid[i, j-1])
            elif i > 0:
                grid[i, j] = min(grid[i, j], grid[i-1, j])
            elif j > 0:
                grid[i, j] = min(grid[i, j], grid[i, j-1])
    
    for i in range(max_recombs+1):
        for j in range(max_bms+1):
            if grid[i, j] == sys.maxsize:
                grid[i, j] = -1

    return grid


def graph_grid(grid, outputfile):
    fig = plt.figure(1)

    ax1 = fig.add_axes([0.15, 0.12, 0.65, 0.75])

    n = np.max(grid)
    colours = ['green', 'blue', 'red']
    nodes = [0.0, 0.5, 1.0]
    bounds = range(0,n+2)

    if (np.min(grid) == -1):
        colours = ['black', 'green', 'blue', 'red']
        nodes = [0.0, 1.0/float(n), 0.5, 1.0]
        bounds = range(-1,n+2)

    cmap = mpl.colors.LinearSegmentedColormap.from_list("mycmap", list(zip(nodes, colours)), N = len(bounds) - 1)

    img = plt.imshow(grid.transpose(), interpolation='nearest',
                     cmap=cmap, aspect='auto', origin='lower')

    plt.title(
        'Fewest recurrent-mutations required \n given recombination and back-mutations')
    plt.xlabel('Recombinations')
    plt.ylabel('Back-mutations')

    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    ax2 = fig.add_axes([0.85, 0.12, 0.05, 0.75])
    fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
                cax=ax2, orientation='vertical',
                label="recurrent mutations")

    fig.savefig(outputfile)


trade_rms_for_recombs = False
def calculate_optimal_curve(arg_runs):
    recombinations = []
    rare_muts = []
    max_recombs = 0
    min_recombs = sys.maxsize
    for (recombs, bms, rms) in arg_runs:
        recombinations.append(recombs)
        rare_muts.append(bms + rms)
        max_recombs = max(max_recombs, recombs)
        min_recombs = min(min_recombs, recombs)

    fewest_rare_muts = [sys.maxsize for i in range(max_recombs+1)]
    for (recombs, bms, rms) in arg_runs:
        fewest_rare_muts[recombs] = min(fewest_rare_muts[recombs], bms+rms)

    for i in range(1, max_recombs+1):
        fewest_rare_muts[i] = min(fewest_rare_muts[i-1], fewest_rare_muts[i])
    
    if trade_rms_for_recombs:
        # Can remove a recurrent mutation by using two recombinations
        for i in range(2, max_recombs+1):
            fewest_rare_muts[i] = max(0, min(fewest_rare_muts[i-2] - 1, fewest_rare_muts[i]))
        
        if (fewest_rare_muts[max_recombs] > 0):
            # Can extend so as to reach all recombinations
            new_max_recombs = max_recombs + 2*fewest_rare_muts[max_recombs]
            for i in range(max_recombs+1, new_max_recombs+1):
                fewest_rare_muts.append(max(0, fewest_rare_muts[i-2] - 1))
            
            max_recombs = new_max_recombs
    
    # Now remove tail only if there is something other than tail
    if fewest_rare_muts[min_recombs] == 0:
        max_recombs = min(max_recombs, min_recombs+10)
    else:
        for i in range(min_recombs, max_recombs):
            if fewest_rare_muts[i] == 0:
                max_recombs = i
                break
    
    fewest_rare_muts = fewest_rare_muts[:max_recombs+1]

    
    return (fewest_rare_muts, min_recombs, max_recombs, recombinations, rare_muts)

def graph_lower_bounds(lower_bound, max_recombs):
    for (bound_type, value) in lower_bound.items():
        if bound_type == "RecMin":
            plt.scatter(x=int(value), y=0, marker='*', s=300, alpha=0.7, label='RecMin')
        elif bound_type == "hudson_overlap":
            plt.scatter(x=int(value), y=0, marker='8', s=300, alpha=0.7, label='Hudson-Kaplan')
        elif bound_type == "RecMinRobust":
            plt.plot([0, int(value)+0.1], [int(value)+0.1, 0], label='RecMin Robust')
        elif bound_type == "hudson_nonoverlap":
            plt.plot([0, int(value)-0.1], [int(value)-0.1, 0], label='1-robust HK')

        elif bound_type == "strong_hudson":
            line = [0 for _ in range(max_recombs)]
            for (parameter, bound) in value:
                for i in range(bound+1):
                    line[i] = max(line[i], (bound - i) * parameter)
            
            end = 0
            while (end < max_recombs and line[end] != 0):
                end += 1

            plt.plot(range(end+1), line[:end+1], '--', alpha=0.7, label='k-strong HK')
        
        elif bound_type == "robust_hudson":
            line = [0 for _ in range(max_recombs)]
            for (parameter, bound) in value:
                for i in range(bound+1):
                    line[i] = max(line[i], (bound - i) * parameter)
            
            end = 0
            while (end < max_recombs and line[end] != 0):
                end += 1

            plt.plot(range(end+1), line[:end+1], ':', label='k-robust HK')
        
        elif bound_type == "random_adversarial_hudson":
            line_dict = {}
            for (p,b) in value:
                if b in line_dict:
                    line_dict[b] = min(p, line_dict[b])
                else:
                    line_dict[b] = p

            plt.plot(list(line_dict.keys()),list(line_dict.values()), '-.', alpha=0.7, label="RNG-adversarial")
        
        elif bound_type == "adversarial_hudson":
            line_dict = {}
            for (p,b) in value:
                if b in line_dict:
                    line_dict[b] = min(p, line_dict[b])
                else:
                    line_dict[b] = p

            plt.plot(list(line_dict.keys()),list(line_dict.values()), '-.', alpha=0.7, label="Adversarial HK")

def graph_all_recombs_vs_rare_mutations(arg_runs, lower_bound, outputfile):
    fig = plt.figure(2, figsize=[8, 6])

    (fewest_rare_muts, min_recombs, max_recombs, recombinations, rare_muts) = calculate_optimal_curve(arg_runs)

    print(fewest_rare_muts)


    plt.plot(range(min_recombs, max_recombs+1), fewest_rare_muts[min_recombs:],
             color='red', alpha=0.7, label="optimal")
    plt.scatter(recombinations, rare_muts, color='blue', alpha=0.7,
                label='networks created')

    
    graph_lower_bounds(lower_bound, max_recombs)


    plt.title('Recombinations vs Recurrent mutations \n for all args created')
    plt.xlabel('Recombinations')
    plt.ylabel('recurrent mutations (all kinds)')
    plt.legend(loc="upper right")

    # plt.show()

    fig.savefig(outputfile, dpi=200, bbox_inches='tight')

def graph_comparison(arg_runs_list, labels, lower_bound, outputfile):
    fig = plt.figure(3, figsize=[4*1.5, 3*1.5])

    lowest_max_recombs = 100000

    for arg_runs, label in zip(arg_runs_list, labels):
        (fewest_rare_muts, min_recombs, max_recombs, recombinations, rare_muts) = calculate_optimal_curve(arg_runs)
        lowest_max_recombs = min(lowest_max_recombs, max_recombs)
        
        plt.plot(range(min_recombs, max_recombs+1), fewest_rare_muts[min_recombs:], label=label, alpha=0.7)
    
    if len(lower_bound) > 0:
        graph_lower_bounds(lower_bound, lowest_max_recombs)


    # plt.title('Recombinations vs Recurrent Mutations')
    plt.xlabel('Recombinations')
    plt.ylabel('recurrent mutations (all kinds)')
    plt.legend(loc="upper right")

    fig.savefig(outputfile, dpi=200, bbox_inches='tight')


def main(argv):
    inputfile = ''
    lowerboundfile = ''
    outputfile = ''
    labels_text = ''
    grid_outputfile = ''
    global trade_rms_for_recombs

    try:
        opts, args = getopt.getopt(argv, "hi:b:o:l:g:t", ["ifile=", "ofile=", "labels=", "gridfile="])
    except getopt.GetoptError:
        print('graph_results.py -i <input file(s)> -b <lowerbound file(s)> -o <output file> -l <labels> -g <output file for grid> -t (smooth curve by trading rms for recombs)')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('graph_results.py -i <input file(s)> -b <lowerbound file(s)> -o <output file> -l <labels> -g <output file for grid> -t (smooth curve by trading rms for recombs)')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-b"):
            lowerboundfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
        elif opt in ("-l", "--labels"):
            labels_text = arg
        elif opt in ("-g", "--gridfile"):
            grid_outputfile = arg
        elif opt in ("-t"):
            trade_rms_for_recombs = True

    print('Input file is: ', inputfile)
    print('lowerbound file is: ', lowerboundfile)
    print('Output file is: ', outputfile)
    print('Labels are: ', labels_text)
    print('Output grid file: ', grid_outputfile)

    if ',' in inputfile:
        # Perform a comparison
        inputfiles = inputfile.split(',')
        arg_runs_list = [read_args(file) for file in inputfiles]

        labels = labels_text.split(',')

        if (len(labels) != len(arg_runs_list)):
            print("number of labels different to number of inputs")
            return

        lower_bound = {}
        if lowerboundfile != '':
            lower_bound = read_lower_bounds(lowerboundfile)
        
        graph_comparison(arg_runs_list, labels_text.split(','), lower_bound, outputfile)
    else:
        arg_runs = read_args(inputfile)
        (fewest_rms, max_recombs, max_bms) = runs_to_dict(arg_runs)

        print(max_recombs)
        print(max_bms)

        lower_bound = {}
        if lowerboundfile != '':
            lower_bound = read_lower_bounds(lowerboundfile)
        
        graph_all_recombs_vs_rare_mutations(arg_runs, lower_bound, outputfile)


        if (grid_outputfile != ''):
            grid = calculate_optimal_grid(fewest_rms, max_recombs, max_bms)
            graph_grid(grid, grid_outputfile)
            print(grid)


        


if __name__ == "__main__":
    main(sys.argv[1:])
