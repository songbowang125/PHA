import datetime
import os
import random
import numpy as np
import pysam
import multiprocessing
import argparse
import sys



def parse_arguments(arguments=sys.argv[1:]):
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', dest="input_folder", type=os.path.abspath, required=True, help='Absolute path to input')
    parser.add_argument('-o', dest="output_folder", type=os.path.abspath, required=True, help='Absolute path to output')

    return parser.parse_args(arguments)


if __name__ == '__main__':

    options = parse_arguments()

    if not os.path.exists(options.output_folder):
        os.mkdir(options.output_folder)

    for sample in os.listdir(options.input_folder):

        sample_final_folder = os.path.join(options.output_folder, sample)
        if not os.path.exists(sample_final_folder):
            os.mkdir(sample_final_folder)

        sample_folder = os.path.join(options.input_folder, sample)

        os.system("cp {}/res_ragtag/ragtag.scaffold.agp {}".format(sample_folder, sample_final_folder))
        # os.system("cp {}/aSAT* {}".format(sample_folder, sample_final_folder))
        os.system("cp {}/aSAT.monomer.bed {}".format(sample_folder, sample_final_folder))

        os.system("cp {}/aSAT.hor.stv.bed {}".format(sample_folder, sample_final_folder))
        os.system("cp {}/aSAT.hor.stv_stats.table {}".format(sample_folder, sample_final_folder))
        os.system("cp {}/aSAT.hor.stv_stats.tsv {}".format(sample_folder, sample_final_folder))

        os.system("cp {}/oSAT* {}".format(sample_folder, sample_final_folder))
        os.system("cp {}/res_fine_aSAT/stats.cent_regions.txt {}".format(sample_folder, sample_final_folder))

        os.system("cp {}/rpmask.out {}".format(sample_folder, sample_final_folder))
