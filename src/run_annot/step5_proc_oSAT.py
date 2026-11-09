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

    return parser.parse_args(arguments)


def run_proc_oSAT(sample_folder):

    sample_other_SAT_folder = os.path.join(sample_folder, "res_other_SAT")


    script_out_path = os.path.join(sample_other_SAT_folder, "run_proc_oSAT.sh")

    sample_other_res_tmp = os.path.join(sample_folder, "oSAT.tmp.bed")
    if os.path.exists(sample_other_res_tmp):
        os.remove(sample_other_res_tmp)

    sample_other_res = os.path.join(sample_folder, "oSAT.bed")
    if os.path.exists(sample_other_res):
        os.remove(sample_other_res)

    sample_rpmask_res = os.path.join(sample_folder, "rpmask.out")
    if os.path.exists(sample_rpmask_res):
        os.remove(sample_rpmask_res)


    with open(script_out_path, "w") as fout:
        fout.write("cd {} \n".format(sample_other_SAT_folder))

        for chrom in os.listdir(sample_other_SAT_folder):

            if "chr" not in chrom:
                continue

            hSat23_res = os.path.join(sample_other_SAT_folder, "{}/{}_RagTag.peri_cent.HSat2and3_Regions.bed".format(chrom, chrom))
            if os.path.exists(hSat23_res):
                fout.write("awk '{print $1\"\t\"$2\"\t\"$3\"\t\"$4\"\t\"$6}' " + hSat23_res + " >> " + sample_other_res_tmp + "\n")

            rpmask_res = os.path.join(sample_other_SAT_folder, "{}/{}_RagTag.peri_cent.fa.out".format(chrom, chrom))
            if os.path.exists(rpmask_res):
                fout.write("grep -w \"SAR\" " +  rpmask_res + " | awk '{print $5\"\t\"$6\"\t\"$7\"\tHSat1A\"\"\t\"$9}' " +  " >> " + sample_other_res_tmp + "\n")
                fout.write("grep -w \"HSATI\" " +  rpmask_res + " | awk '{print $5\"\t\"$6\"\t\"$7\"\tHSat1B\"\"\t\"$9}' " +  " >> " + sample_other_res_tmp + "\n")
                fout.write("grep -E \"BSAT|Beta|LSAU\" " +  rpmask_res + " | awk '{print $5\"\t\"$6\"\t\"$7\"\tBSat\"\"\t\"$9}' " +  " >> " + sample_other_res_tmp + "\n")
                fout.write("grep \"GSAT\" " +  rpmask_res + " | awk '{print $5\"\t\"$6\"\t\"$7\"\tGSat\"\"\t\"$9}' " +  " >> " + sample_other_res_tmp + "\n")

                fout.write("cat {} >> {} \n".format(rpmask_res, sample_rpmask_res))
            else:
                print("Warning: failed rpmask at {}".format(rpmask_res))

            rpmask_tmp1 = os.path.join(sample_other_SAT_folder, "{}/{}_RagTag.peri_cent.fa.masked".format(chrom, chrom))
            rpmask_tmp2 = os.path.join(sample_other_SAT_folder, "{}/{}_RagTag.peri_cent.fa.cat.gz".format(chrom, chrom))

            if os.path.exists(rpmask_tmp1):
                os.remove(rpmask_tmp1)
            if os.path.exists(rpmask_tmp2):
                os.remove(rpmask_tmp2)

        fout.write("bedtools sort -i {} > {} \n".format(sample_other_res_tmp, sample_other_res))

    os.system("bash {}".format(script_out_path))
    os.remove(sample_other_res_tmp)

if __name__ == '__main__':

    options = parse_arguments()

    cnt = 0
    for sample_folder in os.listdir(options.input_folder):

        cnt += 1
        sample_folder = os.path.join(options.input_folder, sample_folder)

        print(sample_folder)

        run_proc_oSAT(sample_folder)

