/*******************************************************************
 *
 *    kwarg_logic.cpp
 *
 *    Implementation of functions to compute a heuristic method to
 *    obtain a (near-)minimal ARG in the presence of both recombination
 *    and recurrent mutation (KwARG).
 *
 ********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include <memory>
#include <iostream>

#include "gene.h"
#include "common.h"
#include "mergesort.h"
#include "bounds.h"
#include "llist.h"
#include "kwarg_logic.h"
#include "beagle_logic.h"
#include "hashtable.h"
#include "bitfunctions.h"
#include "backtrack.h"
#include "common.h"

/* Update global quantities with contribution from g
 */
static double _maxam;
static void update_maxam(Genes *g)
{
    int am = ancestral_material(g);

    if (am > _maxam)
        _maxam = am;
}

/* Score computation for each state in the neighbourhood */
static double sc_min = DBL_MAX, sc_max = 0;
double scoring_function(Genes *g, double step_cost, bool choice_fixed, const RunSettings &run_settings, RunData &run_data, int &lower_bound)
{
    double score;

    // If we have already reached the end, we still cycle through all the possible
    // choices of last step and select the cheapest.
    // We set the score to -(cost of step) if it resolves the last incompatibility,
    // otherwise set the score to -(very big number). The random_select function will
    // pick the move with the least negative score in this case, as needed.
    if (choice_fixed)
    {
        int sign = (run_settings.temp < 0) - (run_settings.temp > 0) - (run_settings.temp == 0);
        if (no_recombinations_required(g, run_data))
            score = sign * step_cost;
        else
            score = sign * DBL_MAX;
    }
    // If we have not reached the end, score the move as usual.
    else
    {
        if (_maxam < 75)
            lower_bound = noexp_rmin(run_data);
        //         else if(_am < 150) {
        //             lower_bound = _eagl(g);
        //         }
        else if (_maxam < 200)
            lower_bound = hb(g, run_data);
        else
            lower_bound = hudson_kaplan_genes(g);

        score = (step_cost + lower_bound) * _maxam + get_am();
        // low score better here

        if (score < sc_min)
            sc_min = score;
        if (score > sc_max)
            sc_max = score;
    }
    return score;
}

double scoring_function(Genes *g, double step_cost, bool choice_fixed, const RunSettings &run_settings, RunData &run_data)
{
    int lb = 0;
    return scoring_function(g, step_cost, choice_fixed, run_settings, run_data, lb);
}

/* Once scores have been computed, renormalise and apply annealing */
double score_renormalise(Genes *g, double score, double step_cost, bool choice_fixed, const RunSettings &run_settings, RunData &run_data)
{
    int sign;

    if (choice_fixed)
    {
        sign = (run_settings.temp < 0) - (run_settings.temp > 0) - (run_settings.temp == 0);
        if (no_recombinations_required(g, run_data))
        {
            score = sign * step_cost;
        }
        else
        {
            score = sign * DBL_MAX;
        }
    }
    else
    {
        if (sc_max != sc_min)
        {
            if (run_settings.temp != -1)
            {
                score = exp(run_settings.temp * (1 - (score - sc_min) / (sc_max - sc_min)));
                // high score is better (flipped from before normalising)
            }
        }
        else
        {
            score = 1;
        }
    }

    return score;
}

/* Update the lookup list of SE/RM and recombination numbers
 * This is of length rec_max, and keeps track of the maximum number of RM events seen
 * for each given number of recombinations already proposed. For example, if we have
 * seen a solution with 5 recombinations and 10 RMs, and the current solution reaches
 * 5 recombinations and 10 RMs but has not yet resolved all incompatibilities, then
 * this solution will be sub-optimal and can be abandoned.
 */
void update_lookup(std::vector<int> &lku, int recombs, int rms)
{
    if (lku.size() <= recombs)
        return;
    if (lku[recombs] <= rms)
        return;
    // Let S = number of SE + RM in the solution
    // Let R = number of recombinations in the solution
    // We use the fact that a rm can be replaced by two recombs
    // Then lookup[R] = S, lookup[R + 1 : R + 2*S] <= S, lookup[R + 2*S : end] = 0
    int k = (lku.size() - 1 > recombs + 2 * rms ? recombs + 2 * rms : lku.size() - 1);
    lku[recombs] = rms;
    for (int i = recombs + 1; i <= k; i++)
    {
        if (lku[i] > rms)
        {
            lku[i] = rms;
        }
    }
    for (int i = k + 1; i < lku.size(); i++)
    {
        lku[i] = 0;
    }
}

const char *names[5] = {
    "Coalescence",
    "Sequencing error",
    "Recurrent mutation",
    "Single recombination",
    "Double recombination"};

static bool _generate_predecessors(std::vector<std::unique_ptr<HistoryFragment>> &predecessors, Genes *g, FILE *print_progress, RunData &path_run_data, const RunSettings &run_settings)
{
    bool choice_fixed = false;
    Action ac;
    _maxam = 0;

    int preds = 0;
    int nbdsize = 0;
    /* Determine interesting recombination ranges */
    Index *start = maximumsubsumedprefixs(g);
    Index *end = maximumsubsumedpostfixs(g);

    auto store_to_fragment = [&](Genes *g, RunData run_data)
    {
        auto f = std::make_unique<HistoryFragment>();

        /* Wrap configuration and events leading to it in a HistoryFragment */
        f->events = run_data.eventlist;
        f->g = g;
        f->step_cost = run_data.current_step_cost;
        f->elements = std::move(run_data.sequence_labels);
        f->sites = std::move(run_data.site_labels);
        f->action = ac;
        if (!f->elements.empty() && g->n != 0 && g->n != f->elements.size())
        {
            fprintf(stderr, "Error: number of sequence labels in sequence_labels [%d] not equal to current size of dataset [%d]. Event type: %.1f", f->elements.size(), g->n, run_data.current_step_cost);
            exit(0);
        }
        if (!f->elements.empty() && g->length > 0 && g->length != f->sites.size())
        {
            fprintf(stderr, "Error: number of site labels in sites not equal to current size of dataset.");
            exit(0);
        }
        if (!choice_fixed && no_recombinations_required(g, run_data))
        {
            /* Found a path to the MRCA - choose it */
            choice_fixed = true;
            //         greedy_choice = f;
        }

        predecessors.push_back(std::move(f));
        update_maxam(g);
    };

    if (g_howverbose > 0)
    {
        fprintf(print_progress, "-------------------------------------------------------------------------------------\n");
        fprintf(print_progress, "Searching possible predecessors:\n");
    }

    path_run_data.current_step_cost = 0;
    ac = COAL;
    coalesce_compatibleandentangled_map(g, path_run_data, store_to_fragment);
    preds = predecessors.size() - nbdsize;
    nbdsize = predecessors.size();
    if (g_howverbose > 0)
    {
        fprintf(print_progress, "%-40s %3d\n", "Coalescing entangled: ", preds);
    }

    if (run_settings.se_cost != -1)
    {
        path_run_data.current_step_cost = run_settings.se_cost;
        ac = SE;

        seqerror_flips(g, path_run_data, store_to_fragment, run_settings);
        preds = predecessors.size() - nbdsize;
        nbdsize = predecessors.size();
        if (g_howverbose > 0)
        {
            fprintf(print_progress, "%-40s %3d\n", "Sequencing errors: ", preds);
        }
    }

    if (run_settings.rm_cost != -1)
    {
        // g_step_cost = run_settings.rm_cost;
        path_run_data.current_step_cost = run_settings.rm_cost;
        ac = RM;

        recmut_flips(g, path_run_data, store_to_fragment, run_settings);

        preds = predecessors.size() - nbdsize;
        nbdsize = predecessors.size();
        if (g_howverbose > 0)
        {
            fprintf(print_progress, "%-40s %3d\n", "Recurrent mutations: ", preds);
        }
    }

    /* Try all sensible events with one split */
    if (run_settings.r_cost != -1)
    {
        // g_step_cost = run_settings.r_cost;
        path_run_data.current_step_cost = run_settings.r_cost;
        ac = RECOMB1;

        maximal_prefix_coalesces_map(g, start, end, path_run_data, store_to_fragment);
        preds = predecessors.size() - nbdsize;
        nbdsize = predecessors.size();
        if (g_howverbose > 0)
        {
            fprintf(print_progress, "%-40s %3d\n", "Prefix recombinations: ", preds);
        }

        maximal_postfix_coalesces_map(g, start, end, path_run_data, store_to_fragment);
        preds = predecessors.size() - nbdsize;
        nbdsize = predecessors.size();
        if (g_howverbose > 0)
        {
            fprintf(print_progress, "%-40s %3d\n", "Postfix recombinations: ", preds);
        }
    }

    /* Try all sensible events with two splits */
    if (run_settings.rr_cost != -1)
    {
        // g_step_cost = run_settings.rr_cost;
        path_run_data.current_step_cost = run_settings.rr_cost;
        ac = RECOMB2;

        maximal_infix_coalesces_map(g, start, end, path_run_data, store_to_fragment);
        preds = predecessors.size() - nbdsize;
        nbdsize = predecessors.size();
        if (g_howverbose > 0)
        {
            fprintf(print_progress, "%-40s %3d\n", "Two recombinations (infix): ", preds);
        }

        maximal_overlap_coalesces_map(g, start, end, path_run_data, store_to_fragment);
        preds = predecessors.size() - nbdsize;
        nbdsize = predecessors.size();
        if (g_howverbose > 0)
        {
            fprintf(print_progress, "%-40s %3d\n", "Two recombinations (overlap): ", preds);
        }
    }

    if (g_howverbose > 0)
    {
        fprintf(print_progress, "%-40s %3d\n", "Finished constructing predecessors.", predecessors.size());
        fprintf(print_progress, "-------------------------------------------------------------------------------------\n");
    }

    free(start);
    free(end);

    return choice_fixed;
}

/* Main function of kwarg implementing neighbourhood search.
 */
Result run_kwarg(Genes *g, FILE *print_progress, int (*select)(double), void (*reset)(void),
                 RunSettings run_settings, RunData &main_path_run_data, double max_run_time)
{
    int nbdsize = 0, total_nbdsize = 0, seflips = 0, rmflips = 0, recombs = 0, preds;
    bool bad_soln = false;
    double r = 0;
    Action ac;
    bool choice_fixed = false;

    /* Store HistoryFragments of possible predecessors in predecessors */
    auto predecessors = std::vector<std::unique_ptr<HistoryFragment>>{};

#ifdef ENABLE_VERBOSE
    int v = verbose();
    set_verbose(0);
#endif

    if (run_settings.rm_max < INT_MAX)
    {
        update_lookup(g_lookup, 0, run_settings.rm_max);
    }

    /* Create working copy of g */
    g = copy_genes(g);

    if (g_howverbose > 0)
    {
        fprintf(print_progress, "Input data:\n");
        if (g_howverbose == 2)
        {
            output_genes(g, print_progress, NULL);
        }
        fprintf(print_progress, "%d sequences with %d sites\n", g->n, g->length);
    }

    // Reduce the dataset
    implode_genes(g, main_path_run_data);
    if (g_howverbose > 0)
    {
        if (g_howverbose == 2)
        {
            output_genes(g, print_progress, NULL);
        }
        printf("%d sequences with %d sites after reducing\n", g->n, g->length);
    }
    if (!g_lookup.empty())
    {
        int site_count = g->n * g->length;
        // Can make history purely with recurrent mutations by making every history direct to the root
        update_lookup(g_lookup, 0, site_count);

        if (site_count < run_settings.rec_max)
        {
            // Can also do history where we have zero root and all-1 child. Everything else comes from recombinations of these two.
            run_settings.rec_max = site_count;
            update_lookup(g_lookup, site_count, 0);
        }
    }

    /* Repeatedly choose an event back in time, until data set has been
     * explained.
     */
    predecessors.clear();

    choice_fixed = no_recombinations_required(g, main_path_run_data);
    if (choice_fixed)
        /* Data set can be explained without recombinations */
        free_genes(g);

    // greedy_choice will point to History fragment to be selected
    
    clock_t tic = clock();

    while (!choice_fixed)
    {
        double timer = (double)(clock() - tic) / CLOCKS_PER_SEC;
        std::cout << "time elapsed on search: " << timer << "\r";
        std::cout.flush();

        if ((max_run_time > 0) && (timer > max_run_time))
        {
            std::cout << "Timed out.                           \n";
            bad_soln = true;
            break;
        }

        auto greedy_choice = std::make_unique<HistoryFragment>();

        choice_fixed = _generate_predecessors(predecessors, g, print_progress, main_path_run_data, run_settings);

        /* Finalise choice and prepare for next iteration */
        free_genes(g);

        /* Still looking for path to MRCA */

        /* So far we have only enumerated putative predecessors -
         * score these and choose one.
         */

        // Set the tracking lists to NULL for the score computation, and destroy the old g_sequence_labels/sites
        reset();

        nbdsize = predecessors.size(); // number of predecessors we score
        if (nbdsize == 0)
        {
            fprintf(stderr, "No neighbours left to search but MRCA not reached.");
            bad_soln = true;
            break;
        }
        total_nbdsize = total_nbdsize + nbdsize;

        // Calculate all the scores and store in an array
        // Update sc_min and sc_max for renormalising the score later
        double score_array[nbdsize];
        if (!choice_fixed)
        {
            sc_min = DBL_MAX, sc_max = 0;
            int i = 0;
            for (auto &f : predecessors)
            {
                RunData scoring_data(true); // Empty RunData that won't track anything
                // output_genes(f->g, stderr, "\nGenes of f:\n");
                reset_beagle_builtins(f->g, scoring_data); // set f to be _greedy_currentstate

                // Calculate all the scores and update the min and max

                score_array[i] = scoring_function(f->g, f->step_cost, choice_fixed, run_settings, scoring_data);

                i++;
            }
        }

        // Now consider each predecessor one by one, score, and set as the new choice if the score is lower
        int i = 0;
        for (auto &f : predecessors)
        {
            RunData scoring_data(true);
            reset_beagle_builtins(f->g, scoring_data); // set _greedy_currentstate to be f->g

            // g_step_cost = f->recombinations;
            double printscore = score_renormalise(f->g, score_array[i], f->step_cost, choice_fixed, run_settings, scoring_data);
            if (print_progress != NULL && g_howverbose == 2)
            {
                fprintf(print_progress, "Predecessor %d obtained with event cost %.1f:\n", i + 1, f->step_cost);
                output_genes(f->g, print_progress, NULL);
                print_int_vector(f->elements, "Sequences: ");
                print_int_vector(f->sites, "Sites: ");
                fprintf(print_progress, "Predecessor score: %.0f \n\n",
                        (printscore == -DBL_MAX ? -INFINITY : (printscore == DBL_MAX ? INFINITY : printscore)));
                fflush(print_progress);
            }
            if (select(printscore))
            {
                // compute score and check if better than that of greedy_choice
                /* Set f to be new choice */
                greedy_choice = std::move(f);
                // output_genes(greedy_choice->g, stderr, "greedy_choice update:\n");
            }

            i++;
        }

        predecessors.clear();

        // greedy_choice is a unique pointer and everything in it will be reset so need to copy out elements.
        if (!choice_fixed)
            g = copy_genes(greedy_choice->g);

        main_path_run_data.sequence_labels = greedy_choice->elements;
        main_path_run_data.site_labels = greedy_choice->sites;

        main_path_run_data.current_step_cost = greedy_choice->step_cost;

        switch (greedy_choice->action)
        {
        case COAL:
            break;
        case SE:
            seflips += greedy_choice->step_cost / run_settings.se_cost;
            break;
        case RM:
            rmflips += greedy_choice->step_cost / run_settings.rm_cost;
            break;
        case RECOMB1:
            recombs++;
            break;
        case RECOMB2:
            recombs += 2;
            break;
        }

        if (print_progress != NULL && g_howverbose == 2)
        {
            fprintf(print_progress, "%s completed at cost of %.3f.\n", names[greedy_choice->action], greedy_choice->step_cost);
            fprintf(print_progress, "-------------------------------------------------------------------------------------\n");
            fprintf(print_progress, "Current data:\n");
            output_genes(greedy_choice->g, print_progress, NULL);
            fflush(print_progress);
        }
        if (print_progress != NULL && g_howverbose == 1)
        {
            fprintf(print_progress, "%s at cost %.3f \n", names[greedy_choice->action], greedy_choice->step_cost);
            fflush(print_progress);
        }

        if (g_use_eventlist && main_path_run_data.eventlist.in_use)
        {
            main_path_run_data.eventlist.append(greedy_choice->events);
        }

        r += greedy_choice->step_cost;

        /* Clean up */

        // Can abandon the run if the number of recombinations already exceeds rec_max
        if (recombs > run_settings.rec_max)
        {
            bad_soln = true;
            std::cout << "rec_max exceeded                                          \n";
            break;
        }

        // Can also abandon the run if the number of SE+RM when we have r recombinations is greater than what we've
        // seen in earlier solutions.
        if (run_settings.rec_max != INT_MAX && !g_lookup.empty())
        {
            if (seflips + rmflips > g_lookup[recombs])
            {
                bad_soln = true;
                std::cout << "  g_lookup exceeded exceeded. rms: " << seflips + rmflips << ", recombs: " << recombs << ", lookup: " << g_lookup[recombs] << "\n";
                break;
            }
        }
    }

    // If we exited the loop because of a sub-optimal solution, record this
    if (bad_soln)
    {
        if (run_settings.run_reference > 0)
        {
            fprintf(print_progress, "%10d %13.0f %6.1f %8.2f %8.2f %8.2f %8.2f %3d %3d %3d %10d ", run_settings.run_reference, run_settings.run_seed, run_settings.temp, run_settings.se_cost, run_settings.rm_cost, run_settings.r_cost, run_settings.rr_cost, -seflips, -rmflips, -recombs, total_nbdsize);
        }
        else
        {
            fprintf(print_progress, "%13.0f %6.1f %8.2f %8.2f %8.2f %8.2f %3d %3d %3d %10d ", run_settings.run_seed, run_settings.temp, run_settings.se_cost, run_settings.rm_cost, run_settings.r_cost, run_settings.rr_cost, -seflips, -rmflips, -recombs, total_nbdsize);
        }

        Result result{
        .seflips = -seflips,
        .rmflips = -rmflips,
        .recombs = -recombs,
        .depth = -1};
        result.run_settings = run_settings;

        return result;
    }
    else
    {
        // Otherwise, record the result
        if (print_progress != NULL && g_howverbose > 0)
        {
            fprintf(print_progress, "\nTotal number of states considered: %d\n", total_nbdsize);
            fprintf(print_progress, "Total event cost: %.1f\n", r);
            if (run_settings.run_reference > 0)
            {
                fprintf(print_progress, "%10s %13s %6s %8s %8s %8s %8s %3s %3s %3s %10s %15s\n", "Ref", "Seed", "Temp", "SE_cost", "RM_cost", "R_cost", "RR_cost",
                        "SE", "RM", "R", "N_states", "Time");
            }
            else
            {
                fprintf(print_progress, "%13s %6s %8s %8s %8s %8s %3s %3s %3s %10s %15s\n", "Seed", "Temp", "SE_cost", "RM_cost", "R_cost", "RR_cost", "SE", "RM", "R", "N_states", "Time");
            }
        }
        if (run_settings.run_reference > 0)
        {
            fprintf(print_progress, "%10d %13.0f %6.1f %8.2f %8.2f %8.2f %8.2f %3d %3d %3d %10d ", run_settings.run_reference, run_settings.run_seed, run_settings.temp, run_settings.se_cost, run_settings.rm_cost, run_settings.r_cost, run_settings.rr_cost, seflips, rmflips, recombs, total_nbdsize);
        }
        else
        {
            fprintf(print_progress, "%13.0f %6.1f %8.2f %8.2f %8.2f %8.2f %3d %3d %3d %10d ", run_settings.run_seed, run_settings.temp, run_settings.se_cost, run_settings.rm_cost, run_settings.r_cost, run_settings.rr_cost, seflips, rmflips, recombs, total_nbdsize);
        }
        if (!g_lookup.empty())
        {
            if (seflips + rmflips < g_lookup[recombs])
            {
                // If found a better bound r < rec_max for Rmin, update.
                if (seflips + rmflips == 0)
                {
                    run_settings.rec_max = recombs;
                }
                update_lookup(g_lookup, recombs, seflips + rmflips);
            }
        }

        Result result{
        .seflips = seflips,
        .rmflips = rmflips,
        .recombs = recombs,
        .depth = 0};
        result.run_settings = run_settings;

        return result;
    }
}

static void _print_depth_indent(FILE *print_progress, int depth)
{
    fprintf(print_progress, (std::string(3 * depth, ' ')).c_str());
}

std::vector<Result> results{};

static void add_result(RunSettings &run_settings, int seflips, int rmflips, int recombs, int depth)
{
    Result result{
        .seflips = seflips,
        .rmflips = rmflips,
        .recombs = recombs,
        .depth = depth};
    result.run_settings = run_settings;
    results.push_back(result);
    

    if (!g_lookup.empty())
    {
        if (seflips + rmflips < g_lookup[recombs])
        {
            // If found a better bound r < rec_max for Rmin, update.
            if (seflips + rmflips == 0)
            {
                run_settings.rec_max = recombs;
            }
            update_lookup(g_lookup, recombs, seflips + rmflips);
        }
    }
}

double _mass_run_tic;
double _max_run_time;
static void _mass_run_kwarg_recursion(Genes *g, FILE *print_progress, std::vector<int> (*take_sample)(std::vector<double>, int),
                                      RunSettings &run_settings, RunData &path_run_data, int num_samples, int seflips, int rmflips, int recombs, int depth)
{
    if (g_howverbose >= 1)
    {
        _print_depth_indent(print_progress, depth);
        fprintf(print_progress, "Recurse with (se, rm, r): (%2d, %2d, %2d) at depth: %2d \n", seflips, rmflips, recombs, depth);
    }

    double timer = (double)(clock() - _mass_run_tic) / CLOCKS_PER_SEC;
    std::cout << "time elapsed on search: " << timer << "\r";
    std::cout.flush();

    if ((_max_run_time > 0) && (timer > _max_run_time))
    {
        std::cout << "Timed out.                           \n";
        return;
    }
    

    // Note we assume that the g passed in belongs to a unique pointer to a HistoryFragment. So we do not free it.

    int nbdsize = 0, total_nbdsize = 0, preds;
    bool bad_soln = false;
    double r = 0;
    Action ac;
    bool choice_fixed = false;

    // Check if we have got to the end (no more recombinations required)
    choice_fixed = no_recombinations_required(g, path_run_data);
    if (choice_fixed)
    {
        if (g_howverbose >= 1)
        {
            _print_depth_indent(print_progress, depth);
            fprintf(print_progress, "Found result(1) with (se, rm, r): (%2d, %2d, %2d) at depth: %2d \n", seflips, rmflips, recombs, depth);
        }
        add_result(run_settings, seflips, rmflips, recombs, depth);
        return;
    }

    // 3. Create nbhd for current state
    /* Store HistoryFragments of possible predecessors in predecessors */
    auto predecessors = std::vector<std::unique_ptr<HistoryFragment>>{};

    choice_fixed = _generate_predecessors(predecessors, g, print_progress, path_run_data, run_settings);

    // 4. Score all predecessors
    nbdsize = predecessors.size(); // number of predecessors we score
    if (nbdsize == 0)
    {
        fprintf(stderr, "No neighbours left to search but MRCA not reached.");
        return;
    }
    total_nbdsize = total_nbdsize + nbdsize;

    // Calculate all the scores and store in an array
    // Update sc_min and sc_max for renormalising the score later
    std::vector<double> score_array(nbdsize);
    std::vector<double> lb_array(nbdsize);
    if (!choice_fixed)
    {
        sc_min = DBL_MAX, sc_max = 0;
        int i = 0;
        for (auto &f : predecessors)
        {
            RunData scoring_data(true); // Empty RunData that won't track anything
            // output_genes(f->g, stderr, "\nGenes of f:\n");
            reset_beagle_builtins(f->g, scoring_data); // set f to be _greedy_currentstate

            // Calculate all the scores and update the min and max

            int lower_bound; // In worse case this is the hudson-kaplan bound
            score_array[i] = scoring_function(f->g, f->step_cost, choice_fixed, run_settings, scoring_data, lower_bound);
            lb_array[i] = lower_bound;

            i++;
        }
    }

    // Now consider each predecessor one by one, score, and set as the new choice if the score is lower
    int i = 0;
    for (auto &f : predecessors)
    {
        RunData scoring_data(true);
        reset_beagle_builtins(f->g, scoring_data); // set _greedy_currentstate to be f->g

        // g_step_cost = f->recombinations;
        double printscore = score_renormalise(f->g, score_array[i], f->step_cost, choice_fixed, run_settings, scoring_data);
        if (print_progress != NULL && g_howverbose == 2)
        {
            fprintf(print_progress, "Predecessor %d obtained with event cost %.1f:\n", i + 1, f->step_cost);
            output_genes(f->g, print_progress, NULL);
            print_int_vector(f->elements, "Sequences: ");
            print_int_vector(f->sites, "Sites: ");
            fprintf(print_progress, "Predecessor score: %.0f \n\n",
                    (printscore == -DBL_MAX ? -INFINITY : (printscore == DBL_MAX ? INFINITY : printscore)));
            fflush(print_progress);
        }
        score_array[i] = printscore;

        i++;
    }

    if (!choice_fixed)
    {
        // Now take samples
        auto samples = take_sample(score_array, num_samples);

        for (int j = 0; j < nbdsize; j++)
        {
            if (samples[j] > 0)
            {
                RunData new_path_data = path_run_data;
                new_path_data.sequence_labels = predecessors[j]->elements;
                new_path_data.site_labels = predecessors[j]->sites;
                new_path_data.current_step_cost = predecessors[j]->step_cost;

                int new_seflips = seflips, new_rmflips = rmflips, new_recombs = recombs;
                switch (predecessors[j]->action)
                {
                case COAL:
                    break;
                case SE:
                    new_seflips += predecessors[j]->step_cost / run_settings.se_cost;
                    break;
                case RM:
                    new_rmflips += predecessors[j]->step_cost / run_settings.rm_cost;
                    break;
                case RECOMB1:
                    new_recombs++;
                    break;
                case RECOMB2:
                    new_recombs += 2;
                    break;
                }

                // Can abandon the run if the number of recombinations already exceeds rec_max
                if (new_recombs > run_settings.rec_max)
                {
                    // bad_soln = true;
                    if (g_howverbose >= 1)
                    {
                        _print_depth_indent(print_progress, depth);
                        fprintf(print_progress, "Abandon run as new_recombs: %3d > rec_max: %3d\n", new_recombs, run_settings.rec_max);
                    }
                    continue;
                }

                // Can also abandon the run if the number of SE+RM when we have r recombinations is greater than what we've
                // seen in earlier solutions.
                if (run_settings.rec_max != INT_MAX && !g_lookup.empty())
                {
                    // if (new_seflips + new_rmflips >= g_lookup[new_recombs])
                    // {
                    //     // bad_soln = true;
                    //     if (g_howverbose >= 1)
                    //     {
                    //         _print_depth_indent(print_progress, depth);
                    //         fprintf(print_progress, "Abandon run as new_seflips + new_rmflips: %3d + %3d >= g_lookup[new_recombs]: %3d\n", new_seflips, new_rmflips, g_lookup[new_recombs]);
                    //     }
                    //     continue;
                    // }

                    if (new_seflips + new_rmflips + (lb_array[j] / 2) > g_lookup[new_recombs])
                    {
                        // bad_soln = true;
                        if (g_howverbose >= 1)
                        {
                            _print_depth_indent(print_progress, depth);
                            fprintf(print_progress, "Abandon run as new_seflips + new_rmflips + (lower_bound/2): %3d + %3d + (%3d/2) >= g_lookup[new_recombs]: (%2d, %3d)\n", new_seflips, new_rmflips, lb_array[j], new_recombs, g_lookup[new_recombs]);
                            _print_depth_indent(print_progress, depth);
                            fprintf(print_progress, "g_lookup: (0, %2d), (1, %2d), (2, %2d), (3, %2d), (4, %2d), (5, %2d), (6, %2d), (7, %2d)\n",
                                g_lookup[0], g_lookup[1], g_lookup[2], g_lookup[3], g_lookup[4], g_lookup[5], g_lookup[6], g_lookup[7]);
                        }
                        continue;
                    }
                }

                _mass_run_kwarg_recursion(predecessors[j]->g, print_progress, take_sample, run_settings, new_path_data, samples[j], new_seflips, new_rmflips, new_recombs, depth + 1);
            }
        }
    }
    else
    {
        // lowest magnitude score is the best option
        double lowest_score = -DBL_MAX;
        int lowest_index = -1;
        for (int i = 0; i < nbdsize; i++)
        {
            if (score_array[i] < 0 && score_array[i] > lowest_score)
            {
                lowest_score = score_array[i];
                lowest_index = i;
            }
        }

        switch (predecessors[lowest_index]->action)
        {
        case COAL:
            break;
        case SE:
            seflips += predecessors[lowest_index]->step_cost / run_settings.se_cost;
            break;
        case RM:
            rmflips += predecessors[lowest_index]->step_cost / run_settings.rm_cost;
            break;
        case RECOMB1:
            recombs++;
            break;
        case RECOMB2:
            recombs += 2;
            break;
        }

        if (recombs > run_settings.rec_max)
        {
            if (g_howverbose >= 1)
            {
                _print_depth_indent(print_progress, depth);
                fprintf(print_progress, "Abandon run as recombs: %3d > max: %3d\n", recombs, run_settings.rec_max);
            }
            return;
        }

        // Can also abandon the run if the number of SE+RM when we have r recombinations is greater than what we've
        // seen in earlier solutions.
        if (run_settings.rec_max != INT_MAX && !g_lookup.empty())
        {
            if (seflips + rmflips > g_lookup[recombs])
            {
                // bad_soln = true;
                if (g_howverbose >= 1)
                {
                    _print_depth_indent(print_progress, depth);
                    fprintf(print_progress, "Abandon run as seflips + rmflips: %3d + %3d >= g_lookup[recombs]: %3d\n", seflips, rmflips, g_lookup[recombs]);
                }
                return;
            }
        }

        if (g_howverbose >= 1)
        {
            _print_depth_indent(print_progress, depth);
            fprintf(print_progress, "Found result(2) with (se, rm, r): (%2d, %2d, %2d) at depth: %2d \n", seflips, rmflips, recombs, depth + 1);
        }
        add_result(run_settings, seflips, rmflips, recombs, depth + 1);
    }
}

std::vector<Result> mass_run_kwarg(Genes *g, FILE *print_progress, std::vector<int> (*take_sample)(std::vector<double>, int),
                                   RunSettings run_settings, RunData &main_path_run_data, int num_samples, double max_run_time)
{
    // Roughly same result as running kwarg num_samples times. Each time a nbhd is created we will take
    // many samples from it rather than just 1. We then take a dfs approach
    results.clear();

    _max_run_time = max_run_time;

#ifdef ENABLE_VERBOSE
    int v = verbose();
    set_verbose(0);
#endif

    // 1. Clean input g
    if (run_settings.rm_max < INT_MAX)
    {
        update_lookup(g_lookup, 0, run_settings.rm_max);
    }

    /* Create working copy of g */
    g = copy_genes(g);

    if (g_howverbose >= 1)
    {
        fprintf(print_progress, "Input data:\n");
        if (g_howverbose == 2)
        {
            output_genes(g, print_progress, NULL);
        }
        fprintf(print_progress, "%d sequences with %d sites\n", g->n, g->length);
    }

    // Reduce the dataset
    implode_genes(g, main_path_run_data);
    if (g_howverbose >= 1)
    {
        if (g_howverbose == 2)
        {
            output_genes(g, print_progress, NULL);
        }
        printf("%d sequences with %d sites after reducing\n", g->n, g->length);
    }
    if (!g_lookup.empty())
    {
        int site_count = g->n * g->length;
        // Can make history purely with recurrent mutations by making every history direct to the root
        update_lookup(g_lookup, 0, site_count);

        if (site_count < run_settings.rec_max)
        {
            // Can also do history where we have zero root and all-1 child. Everything else comes from recombinations of these two.
            run_settings.rec_max = site_count;
            update_lookup(g_lookup, site_count, 0);
        }
    }

    // 2. Start recursion
    _mass_run_tic = clock();
    _mass_run_kwarg_recursion(g, print_progress, take_sample, run_settings, main_path_run_data, num_samples, 0, 0, 0, 0);

    free_genes(g);

    // Now can read from results

    if (print_progress != NULL && g_howverbose >= 1)
    {
        fprintf(print_progress, "Printing %3d Results from DFS search:\n", results.size());
        fprintf(print_progress, "%13s %6s %8s %8s %8s %8s %3s %3s %3s %10s %15s\n", "Seed", "Temp", "SE_cost", "RM_cost", "R_cost", "RR_cost", "SE", "RM", "R", "N_states", "Depth");
    }

    for (Result result : results)
    {
        fprintf(print_progress, "%13.0f %6.1f %8.2f %8.2f %8.2f %8.2f %3d %3d %3d %10d\n",
                run_settings.run_seed, run_settings.temp, run_settings.se_cost, run_settings.rm_cost, run_settings.r_cost, run_settings.rr_cost, result.seflips, result.rmflips, result.recombs, -1, result.depth);
    }

    return results;
}
