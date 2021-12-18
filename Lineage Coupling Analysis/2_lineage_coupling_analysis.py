import codecs
from contextlib import closing
import requests
import csv
import itertools
import statistics
import random
import collections
import copy
import math
import pprint as pp
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
import argparse


def parse_args():
    """Parse the arguments specified by the user."""
    min_num_shared_clones = 1
    min_clonesize = 1
    vmin = -5.0
    vmax = 5.0
    num_shufflings = 100
    csv_data_url = ('https://...')
    csv_data_file = 'Datasets/TIJ_HARMONY_Clusters.csv'
    csv_metric_values_matrix_output_file = (
                    'Results/metric_values_per_state_pair_matrix.csv')
    csv_lineage_coupling_scores_matrix_output_file = (
                                'Results/lineage_coupling_scores_matrix.csv')
    csv_lineage_coupling_scores_correlation_matrix_output_file = (
                    'Results/lineage_coupling_scores_correlation_matrix.csv')
    csv_num_shared_clones_matrix_output_file = (
                                    'Results/num_shared_clones_matrix.csv')
    csv_num_cells_of_shared_clones_matrix_output_file = (
                            'Results/num_cells_of_shared_clones_matrix.csv')
    lineage_coupling_scores_clustermap_output_file = (
                            'Results/lineage_coupling_scores_clustermap.pdf')
    lineage_coupling_scores_correlation_clustermap_output_file = (
                'Results/lineage_coupling_scores_correlation_clustermap.pdf')
    parser = argparse.ArgumentParser(description="Perform"
                    " the Lineage Analysis derived from Wagner et al. (2018)"
                    " for a given input.")
    parser.add_argument("-N", "--num_shufflings", metavar='N'
            , type=int, nargs=1
            , help="The number of random shufflings"
                                    " (of cell-type assignments) to be made"
                                    " (default {}).".format(num_shufflings))
    parser.add_argument("-l", "--csv_url", metavar='csv_data_url'
                            , type=str, nargs=1
                            , help="The csv url containing the data"
                                    " (default '{}').".format(csv_data_url))
    parser.add_argument("-f", "--csv_file", metavar='csv_data_file'
                        , type=str, nargs=1
                        , help="The csv file containing the data"
                                " (default '{}').".format(csv_data_file))
    parser.add_argument("-M"
                , "--csv_metric_values_matrix_output_file"
                , metavar='csv_metric_values_matrix_output_file'
                , type=str, nargs=1
                , help="The output csv file where the"
                        " metric values per state pair"
                        " will be written (note that the folder must exist)"
                        " (default '{}').".format(
                                        csv_metric_values_matrix_output_file))
    parser.add_argument("-Z"
            , "--csv_lineage_coupling_scores_matrix_output_file"
            , metavar='csv_lineage_coupling_scores_matrix_output_file'
            , type=str, nargs=1
            , help="The output csv file where the lineage coupling scores"
                    "  will be written (note that the folder must exist)"
                    " (default '{}').".format(
                            csv_lineage_coupling_scores_matrix_output_file))
    parser.add_argument("-C"
            , "--csv_lineage_coupling_scores_correlation_matrix_output_file"
            , metavar='csv_lineage_coupling_scores'
                                            '_correlation_matrix_output_file'
            , type=str, nargs=1
            , help="The output csv file where the lineage coupling"
                " correlation scores will be written"
                " (note that the folder must exist)"
                " (default '{}').".format(
                    csv_lineage_coupling_scores_correlation_matrix_output_file)
            )
    parser.add_argument("-S"
            , "--csv_num_shared_clones_matrix_output_file"
            , metavar='csv_num_shared_clones_matrix_output_file'
            , type=str, nargs=1
            , help="The output csv file where the numbers of shared clones"
                " per cluster pair will be written"
                " (note that the folder must exist)"
                " (default '{}').".format(
                    csv_num_shared_clones_matrix_output_file)
            )
    parser.add_argument("-T"
            , "--csv_num_cells_of_shared_clones_matrix_output_file"
            , metavar='csv_num_cells_of_shared_clones_matrix_output_file'
            , type=str, nargs=1
            , help="The output csv file where the numbers of cells belonging"
                " to shared clones per cluster pair will be written"
                " (note that the folder must exist)"
                " (default '{}').".format(
                    csv_num_cells_of_shared_clones_matrix_output_file)
            )
    parser.add_argument("-U"
            , "--lineage_coupling_scores_clustermap_output_file"
            , metavar='lineage_coupling_scores_clustermap_output_file'
            , type=str, nargs=1
            , help="The output pdf file where the lineage coupling scores"
                    "  will be written (note that the folder must exist)"
                    " (default '{}').".format(
                    lineage_coupling_scores_clustermap_output_file)
            )
    parser.add_argument("-V"
            , "--lineage_coupling_scores_correlation_clustermap_output_file"
            , metavar='lineage_coupling_scores'
                                        '_correlation_clustermap_output_file'
            , type=str, nargs=1
            , help="The output csv file where the lineage coupling"
                " correlation scores will be written"
                " (note that the folder must exist)"
                " (default '{}').".format(
                    lineage_coupling_scores_correlation_clustermap_output_file)
            )
    parser.add_argument("-s"
            , "--min_num_shared_clones"
            , metavar='min_num_shared_clones'
            , type=int, nargs=1
            , help="The minimum number of shared clones a cluster pair"
                    " must have in order to count its"
                    " observed number of shared clones as > 0"
                    " (default '{}').".format(min_num_shared_clones))
    parser.add_argument("-c"
            , "--min_clonesize"
            , metavar='min_clonesize'
            , type=int, nargs=1
            , help="The minimum number of cells a clone must have in order to"
                    " be taken into account for the analysis"
                    " (default '{}').".format(min_clonesize))
    parser.add_argument("-u", "--clustermap_scale_min", metavar='scale_min'
                    , type=float, nargs=1
                    , help="The minimum number of the scale used"
                            " to plot the clustermap (default {}).".format(
                                                                        vmin))
    parser.add_argument("-v", "--clustermap_scale_max", metavar='scale_max'
                    , type=float, nargs=1
                    , help="The maximum number of the scale used"
                            " to plot the clustermap (default {}).".format(
                                                                        vmax))
    args = parser.parse_args()
    if(args.num_shufflings is not None):
        num_shufflings = args.num_shufflings[0]
    # The url has priority over the path
    if(args.csv_url is not None):
        csv_data_url = args.csv_url[0]
        csv_data_file = None
    elif(args.csv_file is not None):
        csv_data_url = None
        csv_data_file = args.csv_file[0]
    if(args.csv_metric_values_matrix_output_file is not None):
        csv_metric_values_matrix_output_file = (
                    args.csv_metric_values_matrix_output_file[0]
                )
    if(args.csv_lineage_coupling_scores_matrix_output_file is not None):
        csv_lineage_coupling_scores_matrix_output_file = (
                    args.csv_lineage_coupling_scores_matrix_output_file[0])
    if(args.csv_lineage_coupling_scores_correlation_matrix_output_file
                                                            is not None):
        csv_lineage_coupling_scores_correlation_matrix_output_file = (
                args.csv_lineage_coupling_scores_correlation_matrix_output_file[0]
                )
    if(args.csv_num_shared_clones_matrix_output_file is not None):
        csv_num_shared_clones_matrix_output_file = (
                    args.csv_num_shared_clones_matrix_output_file[0]
                )
    if(args.csv_num_cells_of_shared_clones_matrix_output_file is not None):
        csv_num_cells_of_shared_clones_matrix_output_file = (
                    args.csv_num_cells_of_shared_clones_matrix_output_file[0]
                )
    if(args.lineage_coupling_scores_clustermap_output_file is not None):
        lineage_coupling_scores_clustermap_output_file = (
                    args.lineage_coupling_scores_clustermap_output_file[0])
    if(args.lineage_coupling_scores_correlation_clustermap_output_file
                                                            is not None):
        lineage_coupling_scores_correlation_clustermap_output_file = (
                args.lineage_coupling_scores_correlation_clustermap_output_file[0]
                )
    if(args.min_num_shared_clones is not None):
        min_num_shared_clones = args.min_num_shared_clones[0]
    if(args.min_clonesize is not None):
        min_clonesize = args.min_clonesize[0]
    if(args.clustermap_scale_min is not None):
        vmin = args.clustermap_scale_min[0]
    if(args.clustermap_scale_max is not None):
        vmax = args.clustermap_scale_max[0]
    return (num_shufflings
            , csv_data_url, csv_data_file
            , csv_metric_values_matrix_output_file
            , csv_lineage_coupling_scores_matrix_output_file
            , csv_lineage_coupling_scores_correlation_matrix_output_file
            , min_num_shared_clones, min_clonesize, vmin, vmax
            , csv_num_shared_clones_matrix_output_file
            , csv_num_cells_of_shared_clones_matrix_output_file
            , lineage_coupling_scores_clustermap_output_file
            , lineage_coupling_scores_correlation_clustermap_output_file)

def extract_rows(csv_reader, csv_data_clone_cells):
    """
    Extract relevant row values from a csv.DictReader() object.

    Obs:
        In some datasets there may be a column, 'ident_name'
        which stores the names of the cell states.
        That column is used for parts that should output something
        about the cell states (e.g. metric values in the csv).
        However, in all datasets used for the analyses in the paper,
        there wasn't such column, so the column 'ident' itself is used
        as the column name.
    """
    row_tmp = {}
    for i, row in enumerate(csv_reader):
        # Store the register in a new dedicated list.
        # Discard all registers that had an empty value
        # for any of the fields
        # Also, for some reason, if cloneID is converted into a string
        # rather than a integer, then other, unrelated
        # floating computations become a little bit unstable 
        # (e.g. the metric values begin to differ in
        # their latter decimals when computed with
        # the exact same parameters)
        # row_tmp["cloneID"] = int(row["cloneID"])
        row_tmp["cloneID"] = str(row["cloneID"])
        if(row_tmp["cloneID"] == ""):
            print("Register number {}.".format(i))
            print("cloneID is empty."
                    .format(row["cloneID"])
                    )
            continue
        row_tmp["ident_name"] = str(row["ident"])
        if(row_tmp["ident_name"] == ""):
            print("Register number {}.".format(i))
            print("ident_name is empty.\n"
                    .format(row["ident"])
                    )
            continue
        row_tmp["ident"] = str(row["ident"])
        if(row_tmp["ident"] == ""):
            print("Register number {}.".format(i))
            print("ident is empty.\n"
                .format(row["ident"])
                )
            continue
        csv_data_clone_cells.append(row_tmp.copy())
    return csv_data_clone_cells

def parse_csv_data(csv_data_url, csv_data_file):
    """
    Return the relevant parsed data from a selected source.

    The relevant data are constituted by:
        1) Data of cells that belong to a clone,
        without including the name of the cluster.
        2) The names of each cluster.
    """
    csv_data_clone_cells = []
    # Before calling this function, if the selected source was
    # a csv file, then csv_data_url was set to None. Otherwise,
    # the url source has priority.
    if(csv_data_url is not None):
        response = requests.get(csv_data_url)
        with closing(requests.get(csv_data_url, stream=True)) as req:
            csv_reader = csv.DictReader(codecs.iterdecode(
                                                req.iter_lines(), 'utf-8'))
            csv_data_clone_cells = extract_rows(csv_reader
                                                    , csv_data_clone_cells)
    else:
        with open(csv_data_file, newline='') as csvfile:
            csv_reader = csv.DictReader(csvfile)
            csv_data_clone_cells = extract_rows(csv_reader
                                                    , csv_data_clone_cells)
            csvfile.close()
    return csv_data_clone_cells

def compute_clones_distribution(csv_data_clone_cells):
    """Return a counter of cells per clone."""
    clones = collections.Counter()
    for row in csv_data_clone_cells:
        cloneID = row["cloneID"]
        clones[cloneID] += 1
    return clones

def filter_clones(csv_data_clone_cells, clones, min_clonesize):
    """Filter the cells of clones of size lesser than min_clonesize."""
    print("Filtering clones with less than {} cells.".format(min_clonesize))
    print("Originally there are {} rows and {} clones.".format(
                                    len(csv_data_clone_cells), len(clones)))
    clones_to_delete = set()
    rows_to_delete = set()
    num_deleted_rows = 0
    for i, row in enumerate(csv_data_clone_cells):
        if(clones[row["cloneID"]] < min_clonesize):
            num_deleted_rows += 1
            rows_to_delete.add(i)
            clones_to_delete.add(row["cloneID"])

    for row_num in sorted(rows_to_delete, reverse=True):
        del csv_data_clone_cells[row_num]

    for cloneID in clones_to_delete:
        del clones[cloneID]
    
    print("{} deleted rows of {} clones.".format(num_deleted_rows, len(clones_to_delete)))
    print("Total examined rows: {}.".format(i+1))
    # print(clones)
    return csv_data_clone_cells, clones

def compute_clusters_distribution(csv_data_clone_cells):
    """Return a counter of cells per cluster, and their names."""
    clusters = collections.Counter()
    clusters_names = {}
    for row in csv_data_clone_cells:
        clusterID = row["ident"]
        cluster_name = row["ident_name"]
        clusters[clusterID] += 1
        clusters_names[clusterID] = cluster_name
    return clusters, clusters_names
    
def generate_cluster_pairs(clusters_list):
    """
    Return a dictionary-valued dictionary, indexed by the cluster pairs.

    The indices are 2-tuples, whose elements are the clusterIDs.
    Despite this, each pair is conceptually considered to be
    an unordered pair, and hereby appearing only once.
    This dictionary will later store all relevant data
    of each pair of clusters:
        - The metric values for the clusters
        (according to the definition given
        in the Methods section of the paper),
        resulting from the data in the dataset.
        - The parameters of the distribution of metric values
        resulting from doing the shuffling of cell cluster assignments.
        - The z-score of the first with respect to the second.
        - Other data such as num shared clones or num cells of shared
        clones between the cluster pair.
    """
    cluster_pairs = {}
    for cluster_pair in itertools.combinations_with_replacement(
                                                    clusters_list, 2):
        cluster_pairs[cluster_pair] = {}
        cluster_pairs[cluster_pair]["num_shared_clones"] = 0
        cluster_pairs[cluster_pair]["num_cells_in_shared_clones"] = 0
        cluster_pairs[cluster_pair]["metric_value"] = 0.0
    return cluster_pairs

def compute_num_cells_per_clusters_per_clone(csv_data_clone_cells):
    """
    Return a multi-level counter with different cell counts per level.

    Each count is stored in a field named "total".
    What is counted on each level:
        1) Total number of cells.
        2) Total number of cells per clone.
        3) Total number of cells per cluster per clone.
    """
    num_cells = collections.Counter()
    num_cells["total"] = 0
    for cell in csv_data_clone_cells:
        clone = cell["cloneID"]
        cluster = cell["ident"]
        num_cells["total"] += 1
        if(clone not in num_cells):
            num_cells[clone] = collections.Counter()
            num_cells[clone]["total"] = 0
        num_cells[clone]["total"] += 1
        if(cluster not in num_cells[clone]):
            num_cells[clone][cluster] = collections.Counter()
            num_cells[clone][cluster]["total"] = 0
        num_cells[clone][cluster]["total"] += 1
    return num_cells

def compute_num_shared_clones_all_pairs(num_cells, cluster_pairs):
    """
    Return the number of shared clones of each cluster pair.

    (According to the definition given in the Methods section of the
    paper).
    It is stored in a dictionary field, indexed by the cluster pair.
    """
    for clone in num_cells:
        if(clone == "total"):
            continue
        clusters_clone = [key for key in num_cells[clone] if key!="total"]
        clusters_pairs_clone = itertools.combinations(clusters, r=2)
    for cluster_pair_clone in clusters_pairs_clone:
        if(cluster_pair_clone not in cluster_pairs):
            cluster_pair_clone = (cluster_pair_clone[1]
                                        , cluster_pair_clone[0])
        cluster_pairs[cluster_pair_clone]["num_shared_clones"] += 1
    for cluster_clone in clusters_clone:
        if(num_cells[clone][cluster_clone] != 0
                    and num_cells[clone][cluster_clone]["total"] >= 2):
            cluster_pair_clone = (cluster_clone, cluster_clone)
            cluster_pairs[cluster_pair_clone]["num_shared_clones"] += 1
    return cluster_pairs

def compute_cluster_pairs_values(num_cells, cluster_pairs):
    """
    Return different values for each cluster pair.

    The values are (for each cluster pair):
    - Number of shared clones
    - Number of cells of shared clones and in the cluster pair
    - Metric value
    They are stored in a dictionary field, indexed by the cluster pair.
    """
    for clone in num_cells:
        # Skip the field that stores the total number of cells
        if(clone == "total"):
            continue
        clusters_clone = [key for key in num_cells[clone] if key!="total"]
        # Generate pairs of different clusters of the cells in each clone
        cluster_pairs_clone = itertools.combinations(clusters_clone, r=2)
        # And count them/add their cells/add their metric values
        for cluster_pair_clone in cluster_pairs_clone:
            # Extra step needed because at the moment of generating
            # the cluster pair it might have been in the opposite order
            if(cluster_pair_clone not in cluster_pairs):
                cluster_pair_clone = (cluster_pair_clone[1]
                                            , cluster_pair_clone[0])
            # Count this clone as shared between the cluster pair
            cluster_pairs[cluster_pair_clone]["num_shared_clones"] += 1
            # Sum the cells of this shared clone that belonged to the
            # cluster pair
            num_cells_clone_cluster_pair = [None]*2
            num_cells_clone_cluster_pair[0] = num_cells[clone][cluster_pair_clone[0]]
            num_cells_clone_cluster_pair[0] = (
                    0 if num_cells_clone_cluster_pair[0]==0
                    else num_cells_clone_cluster_pair[0]["total"])
            num_cells_clone_cluster_pair[1] = num_cells[clone][
                                                        cluster_pair_clone[1]]
            num_cells_clone_cluster_pair[1] = (
                    0 if num_cells_clone_cluster_pair[1]==0
                    else num_cells_clone_cluster_pair[1]["total"])
            num_cells_clone_cluster_pair = sum(num_cells_clone_cluster_pair)
            cluster_pairs[cluster_pair_clone]["num_cells_in_shared_clones"] += (
                    num_cells_clone_cluster_pair)
            # Sum the % of this shared clone that belonged to the
            # cluster pair
            cluster_pairs[cluster_pair_clone]["metric_value"] += (
                    float(num_cells_clone_cluster_pair)
                    / float(num_cells[clone]["total"]))
        # Do the same as above, but for the pairs
        # composed by the same cluster
        for cluster_clone in clusters_clone:
            # There must be at least 2 cells of a clone in the cluster
            # for that clone to be considered as shared by
            # the cluster with itself
            if(num_cells[clone][cluster_clone] != 0
                    and num_cells[clone][cluster_clone]["total"] >= 2):
                cluster_pair_clone = (cluster_clone, cluster_clone)
                cluster_pairs[cluster_pair_clone]["num_shared_clones"] += 1
                num_cells_clone_cluster_clone = (
                                    num_cells[clone][cluster_clone]["total"])
                cluster_pairs[cluster_pair_clone][
                        "num_cells_in_shared_clones"] += (
                                                num_cells_clone_cluster_clone)
                cluster_pairs[cluster_pair_clone]["metric_value"] += (
                        float(num_cells_clone_cluster_clone)
                        / float(num_cells[clone]["total"]))
    return cluster_pairs

def filter_cluster_pairs_values(cluster_pairs, min_value
                        , value_to_filter_on="num_shared_clones"
                        , value_to_filter="metric_value"):
    """
    Filter values of the cluster that don't reach a minimum thershold.

    The value is specified in value_field, and the values of
    those cluster pairs whose value is less than min_value
    are set to 0.
    """
    for cluster_pair in cluster_pairs:
        if(cluster_pairs[cluster_pair][value_to_filter_on] < min_value):
            cluster_pairs[cluster_pair][value_to_filter] = 0
    return cluster_pairs

def random_sample_with_frequency(freq_table):
    """
    Return a random sample of a given distribution.

    The elements are determined by freq_table, which is expected to be
    a collections.Counter() object
    (and hereby store counts of elements).
    """
    elements = list(freq_table.elements())
    total_elements = len(elements)
    sample = tuple(random.sample(elements, k=total_elements))
    return list(sample)

def print_matrix_onto_csv(array_2d, clusters_names_array_order
    , csv_matrix_output_file):
    """
    Write a matrix of values for each cluster pair into a csv file.
    """
    with open(csv_matrix_output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        # Write the column indexes in the first row
        writer.writerow(["ident_name"] + clusters_names_array_order)
        # Write each row, beginning by its corresponding index
        for cluster_name, row in zip(
                                clusters_names_array_order, array_2d):
            writer.writerow([cluster_name] + list(row))
        csvfile.close()

def write_matrix(title, fieldname, cluster_pairs, clusters, clusters_names
    , csv_matrix_output_file):
    """
    Write a matrix of values of each cluster pair into a csv file.

    Also, return the matrix in a 2D list form,
    and a list with the cluster names.
    Parameters:
        title: string that contains the title of the matrix to be
        outputed. It's just for outputing a message in the console.
        fieldname: in the dictionary cluster_pairs, it's the
        name of the field whose values are desired to be written.
        cluster_pairs: dictionary indexed by the cluster pairs,
        among whose attributes should be fieldname.
        clusters: a collections.Counter() object with the
        frequency distribution of the clusters.
        clusters_names: a list that stores, for each clusterID,
        its corresponding name.
        csv_matrix_output_file: string with the file path where
        the matrix is desired to be written.
    """
    # First, transform the data
    # to a structure accepted by the used libraries.

    # Create a 2-D array form of the values in cluster_pairs,
    # and a list that stores the clusters in the order
    # in which they occur in that array (along the columns or rows)
    matrix, clusters_list_array_order = (
                dict_to_array(cluster_pairs, len(clusters)
                                            , fieldname)
            )
    # print(title + ":")
    # print_2D_float_list(matrix)
    # Create a list that stores the cluster names, in the order
    # in which they occur in the 2D array (along the columns or rows)
    clusters_names_array_order = [
            clusters_names[cluster] for cluster in clusters_list_array_order]
    print_matrix_onto_csv(matrix, clusters_names_array_order
                            , csv_matrix_output_file)

    return matrix, clusters_names_array_order

def compute_distribution_parameters(number_list):
    """Return the mean and standard deviation of a number list."""
    return (float(statistics.mean(number_list))
            , float(statistics.stdev(number_list)))

def compute_distributions(cluster_pairs, num_shufflings
    , csv_data_clone_cells, clusters):
    """
    Return distribution parameters for each cluster pair.

    Do the cell cluster assignment shufflings,
    simulate the resulting assignments, and for each pair of clusters
    compute and store in a list
    the metric values resulting from each shuffling.
    Then compute and store the means and standard deviations
    of the resulting distributions.
    """
    # List that will store the simulated data in each shuffling
    csv_data_clone_cells_local_shuffle = copy.deepcopy(csv_data_clone_cells)
    # Dictionary that will store in lists the metric values
    # obtained over the shufflings for each cluster pair
    cluster_pairs_shufflings = generate_cluster_pairs(clusters.keys())
    for cluster_pair in cluster_pairs_shufflings:
        cluster_pairs_shufflings[cluster_pair]["metric_value"] = []

    for shuffling in range(0, num_shufflings):
        # Perform the random shuffling of cluster assigments
        clusters_assignment = random_sample_with_frequency(clusters)
        # Simulate cell cluster assignments
        for i, csv_data_clone_cell_local in enumerate(
                                        csv_data_clone_cells_local_shuffle):
            csv_data_clone_cell_local["ident"] = clusters_assignment[i]
        # Compute the metric value for each cluster pair
        # observed in the current sample.
        num_cells_per_clusters_per_clone = compute_num_cells_per_clusters_per_clone(
                                        csv_data_clone_cells_local_shuffle)
        # Dictionary, local to each shuffling,
        # that will store the metric values for each cluster pair
        cluster_pairs_local_shuffle = generate_cluster_pairs(clusters.keys())
        cluster_pairs_local_shuffle = compute_cluster_pairs_values(
                num_cells_per_clusters_per_clone, cluster_pairs_local_shuffle)
        # Append the current shuffling's results to the lists
        for cluster_pair in cluster_pairs_local_shuffle:
            cluster_pairs_shufflings[
                    cluster_pair]["metric_value"].append(
                                        cluster_pairs_local_shuffle[
                                                cluster_pair]["metric_value"])

    # Compute, for each cluster pair, the parameters
    # of the metric values distributions resulting from the simulations
    for cluster_pair in cluster_pairs_shufflings:
        (cluster_pairs[cluster_pair]["avg_metric_value"]
                , cluster_pairs[cluster_pair]["stdv_metric_value"]) \
                = compute_distribution_parameters(
                                        cluster_pairs_shufflings[
                                                cluster_pair]["metric_value"])
    return cluster_pairs

def compute_z_scores(cluster_pairs):
    """
    Return the z-scores for each pair of clusters.

    The z-scores are of the metric values of the clusters
    (according to the definition
    given in the Methods section of the paper),
    resulting from the data in the dataset,
    with respect to their corresponding distributions of metric values
    resulting from doing the cell cluster assignments shufflings.
    Obs:
        There is a special case when a pair of clusters
        obtains the same metric value in all shufflings
        (and hereby its distribution's stdev=0
        and the z-score is not defined).
        A "N/A" value is assigned to those clusters in this case
        instead of a z-score value.
        In practice, however, this case is unlikely, and even more
        when N is large.
    """
    cluster_pairs_infinite_z_score = []
    for cluster_pair in cluster_pairs:
        deviation = (float(cluster_pairs[cluster_pair]["metric_value"])
                        -cluster_pairs[cluster_pair]["avg_metric_value"])
        if(cluster_pairs[cluster_pair]["stdv_metric_value"] == 0):
            print("Pair {} had a stdev = 0, with a deviation = {}"
                                            .format(cluster_pair, deviation))
            print(cluster_pairs[cluster_pair])
            if(deviation == 0):
               cluster_pairs[
                    cluster_pair]["metric_value_z_score"] = 0
            else:
                cluster_pairs_infinite_z_score.append(cluster_pair)
        else:
            cluster_pairs[
                    cluster_pair]["metric_value_z_score"] = (
                                    deviation
                                    / cluster_pairs[
                                            cluster_pair]["stdv_metric_value"]
                            )
    for cluster_pair in cluster_pairs_infinite_z_score:
        cluster_pairs[cluster_pair]["metric_value_z_score"] = "N/A"
    return cluster_pairs

def print_2D_float_list(l):
    """Print a 2D list of floats in a row-wise format."""
    for row in l:
        for col in row:
            print("{:8.1f}".format(col), end=" ")
        print("")

def dict_to_array(cluster_pairs, num_clusters, key_name):
    """
    Return a 2D list with a value for each cluster pair, and their IDs.

    The value for each cluster pair is obtained from the key
    named key_name.
    The 2D list to be returned is conceptually considered to be
    a 2D array, with the rows and columns corresponding to
    the first and second clusters of the cluster pairs,
    and viceversa (to fill the whole matrix).
    The cell value will be the corresponding item stored in the field
    key_name.
    The list of row/column IDs is in insertion order, which is the same
    that the order in which they appear in a cluster pair
    when the latter are iterated over.
    Obs:
        The latter means that the order of the cluster IDs (idents)
        is not maintained, but it depends on
        when it appears within a cluster pair,
        and when that cluster pair appears
        in the iteration of cluster pairs.
    """
    # Ordered Dictionary that will store the indexes in the array
    # corresponding to each of the clusters
    array_indexes = collections.OrderedDict()
    max_array_index = -1
    # Initialize the array to a (num_clusters x num_clusters) shape
    array_2D = np.array([
                                [0.0 for i in range(num_clusters)]
                            for j in range(num_clusters)])
    for cluster_pair in cluster_pairs:
        for cluster in cluster_pair:
            if(cluster not in array_indexes):
                max_array_index += 1
                # print("array_indexes[{}] = {}"
                #         .format(cluster, max_array_index))
                # As array_indexes is an Ordered Dict, the order in
                # which the indexes were inserted (and the dict keys
                # created) will correspond to the indexes themselves
                array_indexes[cluster] = max_array_index
        # Fill the matrix with the corresponding value
        array_2D[
                array_indexes[cluster_pair[0]]][
                array_indexes[cluster_pair[1]]] = cluster_pairs[cluster_pair][
                                                                    key_name]
        # In both orders (to fill the whole matrix)
        array_2D[
                array_indexes[cluster_pair[1]]][
                array_indexes[cluster_pair[0]]] = cluster_pairs[cluster_pair][
                                                                    key_name]
    # The keys in insertion order of array_indexes are the
    # clusterIDs in the order that they will appear in the new 2D array
    keys_in_insertion_order = [key for key in array_indexes]
    return array_2D, keys_in_insertion_order

def plot_clustermap(title, output_file, array_2D, labels_array_order
                                                                , vmin, vmax):
    """
    Plot clustermap of the values in array_2D.

    An optional step of hiding the upper triangular matrix
    can be taken.
    """
    # print(title)
    # print_2D_float_list(array_2D)
    mask = np.zeros_like(array_2D)
    # Set the next assigment to False for plotting
    # the whole matrix, and to True for only the lower triangular matrix
    # with the diagonal
    mask[np.triu_indices_from(mask)] = False
    mask[np.diag_indices_from(mask)] = False
    with sns.axes_style("white"):
        ax = sns.clustermap(array_2D
                    , xticklabels = labels_array_order
                    , yticklabels = labels_array_order
                    , method = "average"
                    , mask=mask, vmin=vmin, vmax=vmax
                    # , annot=True
                    # , square=True
                    ,  cmap="bwr")
        plt.setp(ax.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    plt.savefig(output_file)
    plt.show()

def output_lineage_coupling_computations(
    cluster_pairs, clusters, clusters_names
    , csv_lineage_coupling_scores_matrix_output_file
    , lineage_coupling_scores_clustermap_output_file
    , vmin, vmax
    , csv_lineage_coupling_scores_correlation_matrix_output_file
    , lineage_coupling_scores_correlation_clustermap_output_file):
    """
    Output computations related to clusters' lineage couplings.

    The computations of the lineage coupling scores,
    outputed in the form of:
        - a matrix onto a csv
        - and a clustermap in a plot
    """
    # First, transform data to a structure accepted by seaborn methods.

    # The 2D array form of the z-scores,
    # and a list that stores the order of the clusters
    # in which they occur in that array (along the columns or rows)
    z_scores_array, clusters_list_array_order = (
                dict_to_array(cluster_pairs, len(clusters)
                        , "metric_value_z_score")
            )
    # print("Lineage Coupling Matrix:")
    # print_2D_float_list(z_scores_array)
    clusters_names_array_order = [
            clusters_names[cluster] for cluster in clusters_list_array_order]
    print_matrix_onto_csv(z_scores_array, clusters_names_array_order
                , csv_lineage_coupling_scores_matrix_output_file)
    plot_clustermap("Lineage Coupling Scores"
                    , lineage_coupling_scores_clustermap_output_file
                    , z_scores_array, clusters_names_array_order, vmin, vmax)
    z_scores_array_corr = np.corrcoef(z_scores_array)
    print_matrix_onto_csv(z_scores_array_corr, clusters_names_array_order
                , csv_lineage_coupling_scores_correlation_matrix_output_file)
    plot_clustermap("Lineage Coupling Scores Correlation"
                    , lineage_coupling_scores_correlation_clustermap_output_file
                    , z_scores_array_corr, clusters_names_array_order, -1.0, 1.0)

if __name__ == "__main__":
    # Obtain the parameters of the program
    (num_shufflings, csv_data_url, csv_data_file
                , csv_metric_values_matrix_output_file
                , csv_lineage_coupling_scores_matrix_output_file
                , csv_lineage_coupling_scores_correlation_matrix_output_file
                , min_num_shared_clones, min_clonesize , vmin, vmax
                , csv_num_shared_clones_matrix_output_file
                , csv_num_cells_of_shared_clones_matrix_output_file
                , lineage_coupling_scores_clustermap_output_file
                , lineage_coupling_scores_correlation_clustermap_output_file
            ) = parse_args()
    # The url is given priority as the data source over the local file
    csv_data_source_str = (csv_data_url
                                if (csv_data_url is not None)
                                else csv_data_file)
    print("Run 'python lineage_coupling_analysis.py -h' for displaying usage.")
    # Iterable that stores the parsed data
    csv_data_clone_cells = parse_csv_data(csv_data_url, csv_data_file)
    # Counter of cells per clone
    clones = compute_clones_distribution(csv_data_clone_cells)
    
    # Remove cells that are of clones of size lesser than the minimum
    csv_data_clone_cells, clones = filter_clones(csv_data_clone_cells, clones, min_clonesize)

    # Counter of cells per cluster, and names of each cluster
    clusters, clusters_names = compute_clusters_distribution(
                                                        csv_data_clone_cells)
    # Dictionary-valued dictionary indexed by the cluster pairs.
    cluster_pairs = generate_cluster_pairs(clusters.keys())
    # Multi-level counter with different cell counts per level
    num_cells_per_clusters_per_clone = compute_num_cells_per_clusters_per_clone(csv_data_clone_cells)
    cluster_pairs = compute_cluster_pairs_values(
                                            num_cells_per_clusters_per_clone, cluster_pairs)

    cluster_pairs = filter_cluster_pairs_values(cluster_pairs
                                , min_value=min_num_shared_clones 
                                , value_to_filter_on="num_shared_clones"
                                , value_to_filter="num_shared_clones")
    cluster_pairs = filter_cluster_pairs_values(cluster_pairs
                                , min_value=min_num_shared_clones 
                                , value_to_filter_on="num_shared_clones"
                                , value_to_filter="num_cells_in_shared_clones")
    cluster_pairs = filter_cluster_pairs_values(cluster_pairs
                                , min_value=min_num_shared_clones 
                                , value_to_filter_on="num_shared_clones"
                                , value_to_filter="metric_value")
    print("Writing numbers of shared clones per pair into csv file '{}'"
            .format(csv_num_shared_clones_matrix_output_file))
    write_matrix("Matrix of numbers of shared clones per pair"
                , "num_shared_clones", cluster_pairs
                , clusters, clusters_names
                , csv_num_shared_clones_matrix_output_file)
    print("Writing numbers of cells of shared clones in pairs"
            " into csv file '{}'"
            .format(csv_num_cells_of_shared_clones_matrix_output_file))
    write_matrix("Matrix of numbers of cells of shared clones in pairs"
                , "num_cells_in_shared_clones", cluster_pairs
                , clusters, clusters_names
                , csv_num_cells_of_shared_clones_matrix_output_file)
    print("Writing metric values per pair into csv file '{}'"
            .format(csv_metric_values_matrix_output_file))
    write_matrix("Matrix of metric values per pair"
                , "metric_value", cluster_pairs
                , clusters, clusters_names
                , csv_metric_values_matrix_output_file)
    print("Computing lineage coupling scores for data in '{}'"
        ", with {} shufflings, clustermap_minimum {}, and clustermap_maximum {}"
        .format(csv_data_source_str, num_shufflings, vmin, vmax))
    cluster_pairs = compute_distributions(cluster_pairs, num_shufflings
                        , csv_data_clone_cells, clusters)
    cluster_pairs = compute_z_scores(cluster_pairs)    
    # pp.pprint(cluster_pairs)
    output_lineage_coupling_computations(
        cluster_pairs, clusters, clusters_names
        , csv_lineage_coupling_scores_matrix_output_file
        , lineage_coupling_scores_clustermap_output_file
        , vmin, vmax
        , csv_lineage_coupling_scores_correlation_matrix_output_file
        , lineage_coupling_scores_correlation_clustermap_output_file)