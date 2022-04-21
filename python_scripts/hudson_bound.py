from re import sub
from xmlrpc.client import Unmarshaller
import numpy as np
import sys
from os.path import exists
import getopt
import random

from haplotype_bound import RecMin, RecMinWithRobust

from clean_binary_sequences import read_sequences

genes = np.zeros(0)
incomp_grid = np.zeros(0)
num_sequences = 0
num_cols = 0


def read_genes(inputfile):
    global genes
    global num_sequences
    global num_cols

    (_, sequences, num_cols) = read_sequences(inputfile)
    num_sequences = len(sequences)

    print("seqs: ", num_sequences, " cols: ", num_cols)

    genes = np.zeros((num_sequences, num_cols), dtype=int)
    for i in range(num_sequences):
        for j in range(num_cols):
            if(sequences[i][j] == '-'):
                genes[i][j] = -1
            else:
                genes[i][j] = int(sequences[i][j])


def incomp(i, j):
    global genes
    global num_sequences
    # want to check is col i and col j are incompatible
    b00 = False
    b01 = False
    b10 = False
    b11 = False

    for k in range(num_sequences):
        a = genes[k][i]
        b = genes[k][j]
        if(not b00 and a == 0 and b == 0):
            b00 = True
            if b00 and b01 and b10 and b11:
                break
        if(not b10 and a == 1 and b == 0):
            b10 = True
            if b00 and b01 and b10 and b11:
                break
        if(not b01 and a == 0 and b == 1):
            b01 = True
            if b00 and b01 and b10 and b11:
                break
        if(not b11 and a == 1 and b == 1):
            b11 = True
            if b00 and b01 and b10 and b11:
                break

    return b00 and b01 and b10 and b11


def make_incomp_grid():
    global num_cols
    global incomp_grid

    incomp_grid = np.zeros((num_cols, num_cols))
    for i in range(num_cols):
        for j in range(num_cols):
            incomp_grid[i][j] = incomp(i, j)


def hudson(overlap=True):
    regions = []
    i = 0
    j = 1
    while j < num_cols:
        for a in range(i, j):
            if incomp(a, j):
                regions.append((a, j))
                i = j
                if not overlap:
                    i += 1
                    j += 1
                break
        j += 1
    return (len(regions), regions)


def two_hudson():
    # will have a and b incomp with i and j.
    # Can have some overlap: a1-b1=i1 a2 j1 b2=i2-j2, but can't allow overlap of numbers
    # Idea is that to remove any of these recombinations needs two recurrent mutations
    regions = []
    left = 0
    i = 2
    prev_j = -1
    while i < num_cols - 1:
        ab_incomps = []
        for a in range(left, i):
            if a == prev_j:
                continue  # Can't have a or b same as previous j
            if(incomp(a, i)):
                ab_incomps.append(a)

        if(len(ab_incomps) >= 2):
            for j in range(i+1, num_cols):
                x = 0
                a = -1
                b = -1

                # First try an find a, this is allowed to be before previous j
                while x < len(ab_incomps):
                    if(incomp(ab_incomps[x], j)):
                        a = ab_incomps[x]
                        x += 1
                        break
                    x += 1

                # Now look for b which must be after previous j
                while (x < len(ab_incomps) and ab_incomps[x] < prev_j):
                    x += 1

                while x < len(ab_incomps):
                    if(incomp(ab_incomps[x], j)):
                        b = ab_incomps[x]
                        break
                    x += 1

                if(a != -1 and b != -1):
                    regions.append((a, b, i, j))
                    left = i + 1
                    i = j + 1  # i will increase once more at end of iteration
                    prev_j = j
                    break

        i += 1

    return (len(regions), regions)


def two_hudson2():
    # This version aims to make j minimal rather than i
    regions = []
    left = 0
    j = 3
    prev_j = -1
    while j < num_cols:
        ab_incomps = [a for a in range(
            left, j-1) if a != prev_j and incomp(a, j)]

        if(len(ab_incomps) >= 2):
            for i in range(j-1, max(prev_j + 2, ab_incomps[1]+1), -1):
                x = 0
                a = -1
                b = -1

                # First try an find a, this is allowed to be before previous j
                while x < len(ab_incomps):
                    if(incomp(ab_incomps[x], i)):
                        a = ab_incomps[x]
                        x += 1
                        break
                    x += 1

                # Now look for b which must be after previous j
                while (x < len(ab_incomps) and ab_incomps[x] < prev_j):
                    x += 1

                while x < len(ab_incomps):
                    if(incomp(ab_incomps[x], i)):
                        b = ab_incomps[x]
                        break
                    x += 1

                if(a != -1 and b != -1):
                    regions.append((a, b, i, j))
                    left = i + 1
                    prev_j = j
                    j += 2  # j will increase once more at end of iteration
                    break

        j += 1

    return (len(regions), regions)


def find_br_down(k, r, start_search, current_incomps):
    # have found bk..b_(r+1) (b_r = b_index) and know current_incomps are incompatible with all of them

    for br in range(start_search, 0, -1):
        subset = [i for i in current_incomps if incomp(i, br)]
        if(len(subset) >= k and subset[k-1] <= br - r):
            if(r == 1):
                # done
                return ([br], subset[0:k])
            else:
                (bs, ays) = find_br_down(k, r-1, br - 1, subset)
                if(len(bs) > 0):
                    bs.append(br)
                    return(bs, ays)

    return ([], [])


def k_strong_hudson(k):
    # will have a1 to ak incompatible with b1 to bk.
    # Each region must be completely separate
    # To make them smallest we actually base it on bk
    regions = []
    left = 0
    bk = 2*k - 1
    while bk < num_cols:
        a_incomps = [a for a in range(left, bk) if incomp(a, bk)]
        if(len(a_incomps) >= k and a_incomps[k-1] <= bk - k):
            # Checks if there is enough of a gap to fit b1,...,b(k-1)

            (b_list, a_list) = find_br_down(k, k-1, bk-1, a_incomps)
            if(len(b_list) > 0):
                b_list.append(bk)
                regions.append((a_list, b_list))
                left = bk+1
                bk = left + 2*k - 2

        bk += 1

    return (len(regions), regions)


def find_disjoint_pairs(all_pairs, current_pairs, number_needed):
    if number_needed == 0:
        return current_pairs

    if len(all_pairs) == 0:
        return []

    current_points = {i for pair in current_pairs for i in pair}
    
    (i, j) = all_pairs[0]
    if i in current_points or j in current_points:
        return find_disjoint_pairs(all_pairs[1:], current_pairs[:], number_needed)
    else:
        result = find_disjoint_pairs(all_pairs[1:], current_pairs + [(i,j)], number_needed - 1)
        if result == []:
            result = find_disjoint_pairs(all_pairs[1:], current_pairs[:], number_needed)
        
        return result

def find_closest_disjoint_pairs(all_pairs, number_needed):
    result = find_disjoint_pairs(all_pairs, [], number_needed)
    if result == []:
        return []

    best_result = result
    best_right = max([j for (i,j) in result])

    while (True):
        all_pairs = [(i,j) for (i,j) in all_pairs if j < best_right]
        result = find_disjoint_pairs(all_pairs, [], number_needed)
        
        if result != []:
            best_result = result
            best_right = max([j for (i,j) in result])
        else:
            break
    
    return best_result


def k_robust_hudson(k: int):
    global num_cols

    # Searches for pairs 
    regions = [] # region formed of (start, end, [pairs])

    left = 0
    middle = 0 # pairs will all go across middle

    while middle < num_cols - k:
        best_right : int = int(num_cols)
        best_pairs = []
        active_incomps = []

        while middle < best_right - k:
            # remove incomps that are no longer active
            active_incomps = [(i, j) for (i, j) in active_incomps if (j > middle and j < best_right)]

            # add new incomps
            new_incomps = [(middle, j) for j in range(middle+1, best_right) if incomp(middle, j)]
            active_incomps.extend(new_incomps)

            if (len({i for (i,j) in active_incomps}) >= k):
                # Now check if k disjoint pairs can be found 
                pairs = find_closest_disjoint_pairs(active_incomps, k)

                if pairs != []:
                    best_right = max([j for (i,j) in pairs])
                    best_pairs = pairs

            middle += 1
        
        if best_pairs != []:
            regions.append((left, best_right, best_pairs))

        left = best_right + 1
        middle = left
    
    return (len(regions), regions)



# Tail overlap section is problematic as it allows too much


def masked_hudson(mask):
    global num_cols
    global incomp_grid

    unmasked = [a for a in range(num_cols) if not (a in mask)]

    regions = []
    left = 0

    for j in unmasked:
        for i in range(left, j):
            if i in mask:
                continue
            if incomp_grid[i][j]:
                regions.append((i, j))
                left = j
                break

    return (len(regions), regions)


def split_out_dups(sites_used):
    sites_used_no_dups = []
    sites_used_dups = []
    # how many site masking events could remove two recombinations
    potential_dups_removable = 0

    prev_was_removable_dup = False
    i = 1 
    while i < len(sites_used):
        if sites_used[i-1] == sites_used[i]:
            sites_used_dups.append(sites_used[i])
            i += 2
            if not prev_was_removable_dup:
                potential_dups_removable += 1
                prev_was_removable_dup = True
            else:
                prev_was_removable_dup = False
        else:
            sites_used_no_dups.append(sites_used[i-1])
            prev_was_removable_dup = False
            i += 1
    sites_used_no_dups.append(sites_used[-1])

    print("no_dups: ", sites_used_no_dups)
    print("dups: ", sites_used_dups)
    print("removable: ", potential_dups_removable)


_target_ratio = 2.0

def crude_adversarial_hudson_iter(current_mask, m, ignored_sites):
    global _target_ratio
    # can place m more recurrent mutations
    (num_regions, regions) = masked_hudson(current_mask)
    if (m == 0 or num_regions == 0):
        return (num_regions, regions, current_mask)

    sites_used = [i for ij in regions for i in ij]
    # may have duplicates but not a problem

    best_num_regions = sys.maxsize
    best_regions = []
    best_mask = []

    target = max(0, float(num_regions) - (2 * m)) # Best that can be done is to 
    for site in sites_used:
        if site in ignored_sites:
            continue

        (new_num_regions, new_regions, new_mask) = crude_adversarial_hudson_iter(current_mask + [site], m-1, ignored_sites[:])
        if (new_num_regions < best_num_regions):
            # Adverserial so want to minimize
            best_num_regions = new_num_regions
            best_regions = new_regions
            best_mask = new_mask

            if best_num_regions <= target:
                break

        ignored_sites.append(site)

    return [best_num_regions, best_regions, best_mask]


def crude_adversarial_hudson(m):
    # adversarial is allowed to select m sites to place recurrent mutations on
    return crude_adversarial_hudson_iter([], m, [])


def target_adversarial_hudson_iter(current_mask, m, ignored_sites):
    global _target_ratio
    # can place m more recurrent mutations
    (num_regions, regions) = masked_hudson(current_mask)
    if (m == 0 or num_regions == 0):
        return (num_regions, regions, current_mask)

    sites_used = [i for ij in regions for i in ij]
    # may have duplicates but not a problem

    sites_used_no_dups = []
    sites_used_dups = []
    # how many site masking events could remove two recombinations
    potential_dups_removable = 0

    prev_was_removable_dup = False
    i = 1
    while i < len(sites_used):
        if sites_used[i-1] == sites_used[i]:
            sites_used_dups.append(sites_used[i])
            i += 2
            if not prev_was_removable_dup:
                potential_dups_removable += 1
                prev_was_removable_dup = True
            else:
                prev_was_removable_dup = False
        else:
            sites_used_no_dups.append(sites_used[i-1])
            prev_was_removable_dup = False
            i += 1
    sites_used_no_dups.append(sites_used[-1])

    best_num_regions = sys.maxsize
    best_regions = []
    best_mask = []

    target = num_regions - (2 * min(m, potential_dups_removable)) - max(0, m - potential_dups_removable)
    sites_considered = 0
    for site in sites_used_dups + sites_used_no_dups:
        if site in ignored_sites:
            continue

        (new_num_regions, new_regions, new_mask) = target_adversarial_hudson_iter(current_mask + [site], m-1, ignored_sites[:])
        if (new_num_regions < best_num_regions):
            # Adverserial so want to minimize
            best_num_regions = new_num_regions
            best_regions = new_regions
            best_mask = new_mask

            if best_num_regions <= target:
                break

        ignored_sites.append(site)

    return [best_num_regions, best_regions, best_mask]


def target_adversarial_hudson(m):
    # adversarial is allowed to select m sites to place recurrent mutations on
    return target_adversarial_hudson_iter([], m, [])


def random_adversarial_hudson_iter(current_mask, m):
    # can place m more recurrent mutations
    (num_regions, regions) = masked_hudson(current_mask)
    if (m == 0 or num_regions == 0):
        return (num_regions, regions, current_mask)

    sites_used = [i for ij in regions for i in ij]

    site = random.choice(sites_used)

    return random_adversarial_hudson_iter(current_mask + [site], m-1)


def random_adversarial_hudson(m):
    # adversarial is allowed to select m sites to place recurrent mutations on

    best_num_regions = sys.maxsize
    best_regions = []
    bset_mask = []

    for i in range(300):
        (num_regions, regions, mask) = random_adversarial_hudson_iter([], m)
        if (num_regions < best_num_regions):
            best_num_regions = num_regions
            best_regions = regions
            best_mask = mask

    return [best_num_regions, best_regions, best_mask]


def output_standard(outfile, calculate_hudson, robust_k_max, strong_k_max, adversarial_max_m, random_adversarial_max_m, calculate_recmin):
    global genes

    if not exists(outfile):
        file = open(outfile, 'w')
        file.write("Bound type, bound value, parameter (k/m)\n")
        file.close()
    
    file = open(outfile, 'a')

    if calculate_hudson:
        hudson_bound, _ = hudson(True)
        file.write("hudson_overlap" + "," + str(hudson_bound) + "," + "\n")
        print("hudson found")
    
    if calculate_recmin:
        haplotype_bound, robust_haplotype_bound = RecMinWithRobust(genes, 12, 10)
        file.write("RecMin" + "," + str(haplotype_bound) + "," + "\n")
        file.write("RecMinRobust" + "," + str(robust_haplotype_bound) + "," + "\n")
        print("RecMin found")

    if calculate_hudson:
        hudson_bound, _ = hudson(False)
        file.write("hudson_nonoverlap" + "," + str(hudson_bound) + "," + "\n")
        print("hudson non overlap found")
    
    for k in range(2, robust_k_max+1):
        robust_bound, _ = k_robust_hudson(k)
        file.write("robust_hudson" + "," + str(robust_bound) + "," + str(k) + "\n")
        if robust_bound == 0:
            break
    print("robust hudsons found")
    
    for k in range(2, strong_k_max+1):
        strong_bound, _ = k_strong_hudson(k)
        file.write("strong_hudson" + "," + str(strong_bound) + "," + str(k) + "\n")
        if strong_bound == 0:
            break
    print("strong hudsons found")
    
    if adversarial_max_m >= 0 or random_adversarial_max_m >= 0:
        make_incomp_grid()

    for m in range(0, adversarial_max_m+1):
        adversarial_bound, _, _ = crude_adversarial_hudson(m)
        file.write("adversarial_hudson" + "," + str(adversarial_bound) + "," + str(m) + "\n")
        if adversarial_bound == 0:
            break
    print("adversarial hudsons found")

    for m in range(0, random_adversarial_max_m+1):
        adversarial_bound, _, _ = random_adversarial_hudson(m)
        file.write("random_adversarial_hudson" + "," + str(adversarial_bound) + "," + str(m) + "\n")
        if adversarial_bound == 0:
            break
    print("random adversarial hudsons found")

    file.close()

def main(argv):
    inputfile = ''
    outputfile = ''
    calculate_hudson = False
    robust_max = -1
    strong_max = -1
    adversarial_max = -1
    random_adversarial_max = -1
    calculate_recmin = False

    try:
        opts, args = getopt.getopt(argv, "i:o:hr:s:a:b:m", ["ifile=","ofile="])
    except getopt.GetoptError:
        print('hudson_bound.py -i <input file> -o <output file> -h (calculate hudson) -r <robust max k> -s <strong max k> -a <adversarial max m> -b <random adversarial max m> -m (calculate recmin)')
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
        elif opt in ("-h"):
            calculate_hudson = True
        elif opt in ("-r"):
            robust_max = int(arg)
        elif opt in ("-s"):
            strong_max = int(arg)
        elif opt in ("-a"):
            adversarial_max = int(arg)
        elif opt in ("-b"):
            random_adversarial_max = int(arg)
        elif opt in ("-m"):
            calculate_recmin = True

    print('Input file is: ', inputfile)
    print('Output file is: ', outputfile)

    read_genes(inputfile)
    # Input gets read to global variables

    if (outputfile != ''):
        output_standard(outputfile, calculate_hudson, robust_max, strong_max, adversarial_max, random_adversarial_max, calculate_recmin)
        return


if __name__ == "__main__":
    main(sys.argv[1:])
