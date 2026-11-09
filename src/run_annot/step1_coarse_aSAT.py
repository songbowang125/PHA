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
    parser.add_argument('-s', dest="code_path", type=os.path.abspath, required=True, help='Absolute path to codes')

    parser.add_argument('-t', dest="thread", type=int, required=True, help='Thread number')
    parser.add_argument('-n', dest='nodes', metavar="IN", nargs="+", type=str, help='Cluster nodes')
    parser.add_argument('-c', dest='cluster', type=str, required=True, choices=['pbs', 'slurm'], help='Cluster ID')

    return parser.parse_args(arguments)



def run_coarse_aSAT(in_fasta, sample_output_folder, code_path, thread, nodes, cluster):
    if not os.path.exists(sample_output_folder):
        os.mkdir(sample_output_folder)

    nhmmer_tbl_out = os.path.join(sample_output_folder,  "coarse_aSAT.tbl")
    nhmmer_bed_out = os.path.join(sample_output_folder, "coarse_aSAT.bed")
    nhmmer_bed_10k_out = os.path.join(sample_output_folder, "coarse_aSAT.10k.bed")

    hmm_model = os.path.join(code_path, "hmm_models/DF_ALR_DF000000029.hmm")

    if cluster == "pbs":
        script_out_path = os.path.join(sample_output_folder, "run_coarse_aSAT.pbs")
    else:
        script_out_path = os.path.join(sample_output_folder, "run_coarse_aSAT.slurm")

    with open(script_out_path, "w") as fout:
        if cluster == "pbs":

            fout.write("#!/bin/bash\n")
            fout.write("#PBS -N {}\n".format("run_coarse_aSAT"))
            fout.write("#PBS -e {}\n".format(os.path.join(sample_output_folder, "run_coarse_aSAT.elog")))
            fout.write("#PBS -o {}\n".format(os.path.join(sample_output_folder, "run_coarse_aSAT.olog")))
            fout.write("#PBS -q {}\n".format("batch"))
            fout.write("#PBS -l nodes={}:ppn={}\n".format(random.choice(nodes), thread))
            fout.write("##PBS -V\n")
        else:
            fout.write("#!/bin/bash\n")
            fout.write("#SBATCH --job-name={}\n".format("run_coarse_aSAT"))
            fout.write("#SBATCH --partition=CPUqueue_fcp1\n")
            fout.write("#SBATCH --account=songbo\n")
            fout.write("#SBATCH --qos=songbo_1_qos\n")
            fout.write("#SBATCH --nodes=1\n")
            fout.write("#SBATCH --ntasks-per-node=1\n")
            fout.write("#SBATCH --cpus-per-task={}\n".format(thread))
            fout.write("#SBATCH --mem={}G\n".format(int(500 / (64 / thread))))
            fout.write("#SBATCH -o {}\n".format(os.path.join(sample_output_folder, "run_coarse_aSAT.olog")))
            fout.write("#SBATCH -e {}\n".format(os.path.join(sample_output_folder, "run_coarse_aSAT.elog")))
            fout.write("#SBATCH -w {}\n".format(random.choice(nodes)))

        fout.write("start_time=$(date +%s)" + "\n")

        # # run nhmmer for coarse aSAT search
        cmd_str = "nhmmer -o /dev/null --notextw --noali --dna --cpu {} --tblout {} {} {}".format(thread, nhmmer_tbl_out, hmm_model, in_fasta)
        fout.write(cmd_str + "\n")

        # # convert tbl to bed
        cmd_str = "awk 'BEGIN{OFS=\"\\t\"}{start = ($7 < $8 ? $7 : $8);end   = ($7 < $8 ? $8 : $7); print $1, start, end, $3}' " + nhmmer_tbl_out + " | bedtools sort | bedtools merge > " + nhmmer_bed_out
        fout.write(cmd_str + "\n")

        cmd_str = ("bedtools merge -d 10000 -i {} > {}".format(nhmmer_bed_out, nhmmer_bed_10k_out))
        fout.write(cmd_str + "\n")

        fout.write("end_time=$(date +%s)" + "\n")
        fout.write("elapsed_seconds=$((end_time - start_time))" + "\n")
        fout.write("elapsed_minutes=$((elapsed_seconds / 60))" + "\n")
        fout.write("echo \"Total $elapsed_seconds Seconds, $elapsed_minutes Minutes\"" + "\n")

    # os.system("qsub {}".format(script_out_path))
    # os.system("bash {}".format(pbs_out_path))
    if cluster == "pbs":
        os.system("qsub {}".format(script_out_path))
    else:
        os.system("sbatch {}".format(script_out_path))

if __name__ == '__main__':



    options = parse_arguments()

    for sample_folder in os.listdir(options.input_folder):

        sample_coarse_aSAT_folder = os.path.join(options.input_folder, "{}/res_coarse_aSAT".format(sample_folder))

        sample_ragtag_fa = os.path.join(options.input_folder, "{}/res_ragtag/ragtag.scaffold.fasta".format(sample_folder))

        run_coarse_aSAT(sample_ragtag_fa, sample_coarse_aSAT_folder, options.code_path, options.thread, options.nodes, options.cluster)

        # print(sample_folder)
        # break

