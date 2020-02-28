#!/usr/bin/env python3

"""Takes kmer a set of kmer count files and calculates the distance between them.
"""

import argparse
import sys
import time
import os
import re
from scipy.cluster.hierarchy import dendrogram, linkage
import seaborn as sns
import matplotlib.pyplot as plt
import collections as c
import numpy as np

def read_kmer_array(kmer_count_file):

    lbc=25
    ubc=75

    kmer_array = c.defaultdict(int)
    ksa = []

    with open(kmer_count_file, "r") as f:
        argc = f.readline()
        header = f.readline()

        for line in f.readlines():
            kmer, count, _ = re.split("[\t|\n]", line)
            if "N" in kmer:
                continue
            # kmer_array[kmer] = int(count)
            ksa.append([kmer, count])

    ksa_s = sorted(ksa, key=lambda x: x[1], reverse=True)
    lb = int((len(ksa_s)/100)*lbc)
    ub = int((len(ksa_s)/100)*ubc)

    for kmer, count in ksa_s[lb:ub+1]:
        kmer_array[kmer] = int(count)

    return kmer_array

def ed(m1, m2):

	a = []
	for f, c in m1.items():
		if f in m2:
			a.append([c, m2[f]])
		else:
			a.append([c, 0])
	return distance_calc(a)

def distance_calc(dataset):
	""" Calculates the Euclidian Distance of the target from 0 """
	ED = np.sqrt(np.sum([float(i[0]-i[1])**2 for i in dataset]))
	return ED

def kmer_array_Dn(kmer_master_dict):

	ED_dict = c.defaultdict(dict)

	for i1, m1 in kmer_master_dict.items():
		for i2, m2 in kmer_master_dict.items():
			ED = ed(m1, m2)
			ED_dict[i1][i2]=ED
	return ED_dict

def heatmap(ED_dict):

    ht_array = []

    for key, value in ED_dict.items():
        b=[v for k, v in value.items()]
        ht_array.append(np.array(b))

    print(np.array(ht_array))

    # fig, ax = plt.subplots()
    # im = ax.imshow(ht_array)
    #
    # ax.set_xticks(np.arange(len(ht_array)))
    # ax.set_yticks(np.arange(len(ht_array)))
    #
    # tl = [k for k, v in ED_dict.items()]
    #
    # ax.set_xticklabels(tl)
    # ax.set_yticklabels(tl)
    #
    # plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
    # rotation_mode="anchor")
    #
    # # for i in range(len(tl)):
    # #     for j in range(len(tl)):
    # #         text = ax.text(j, i, ht_array[i, j],
    # #                        ha="center", va="center", color="w")
    #
    plt.tight_layout()

    # sns.clustermap(ht_array)


    plt.style.use('ggplot')

    Z = linkage(ht_array, 'ward')

    plt.title('Hierarchical Clustering Dendrogram')
    plt.xlabel('ncov19 accession (clinical isoalte)')
    plt.ylabel('distance (Ward)')
    dendrogram(Z, labels=[key for key, value in ED_dict.items()], leaf_rotation=90)

    plt.show()

def is_file(filename):
    """Checks if a path is a file"""

    if not os.path.isfile(filename):

        msg = "{0} is not a file".format(filename)
        raise argparse.ArgumentTypeError(msg)

    else:
        return os.path.abspath(os.path.realpath(os.path.expanduser(filename)))

def is_dir(direname):
    """Checks if a path is a directory"""

    if not os.path.isdir(direname):

        msg = "{0} is not a directory".format(direname)
        raise argparse.ArgumentTypeError(msg)

    else:
        return os.path.abspath(os.path.realpath(os.path.expanduser(direname)))

def parse_args(argv):

    parser = argparse.ArgumentParser()

    parser.add_argument('script_path', action='store', help=argparse.SUPPRESS)
    parser.add_argument('kmer_array', type=is_file, nargs='*', action='store', help='Kmer array files.')

    args = parser.parse_args(argv)

    return args

def main(argv):

    args = parse_args(argv)

    kmer_master_dict = c.defaultdict()

    for kmer_count_file in args.kmer_array:
        kmer_array = read_kmer_array(kmer_count_file)
        fname = os.path.basename(kmer_count_file).split(".")[0]
        kmer_master_dict[fname] = kmer_array

    ED_dict = kmer_array_Dn(kmer_master_dict)

    heatmap(ED_dict)

if __name__ == "__main__":
    if len(sys.argv) < 2 or "-h" in sys.argv:
        print(__doc__)
        main([sys.argv[0], "-h"])
        sys.exit(1)
    main(sys.argv)
