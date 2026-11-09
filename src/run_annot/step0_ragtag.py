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
    parser.add_argument('-o', dest="output_folder", type=os.path.abspath, required=True, help='Absolute path to output')
    parser.add_argument('-i', dest="input_folder", type=os.path.abspath, required=True, help='Absolute path to input')
    parser.add_argument('-r', dest="ref_path", type=os.path.abspath, required=True, help='Absolute path to reference')

    parser.add_argument('-t', dest="thread", type=int, required=True, help='Thread number')
    parser.add_argument('-n', dest='nodes', metavar="IN", nargs="+", type=str, help='Cluster nodes')
    parser.add_argument('-c', dest='cluster', type=str, required=True, choices=['pbs', 'slurm'], help='Cluster ID')

    return parser.parse_args(arguments)

def run_ragtag_for_one_fa(ref_path, sample_name, sample_path, output_folder, thread, nodes, cluster):

    sample_output_folder = os.path.join(output_folder, "{}/res_ragtag".format(sample_name))

    os.system("mkdir -p {}\n".format(sample_output_folder))

    if cluster == "pbs":
        script_out_path = os.path.join(sample_output_folder, "run_ragtag.pbs")
    else:
        script_out_path = os.path.join(sample_output_folder, "run_ragtag.slurm")

    with open(script_out_path, "w") as fout:

        if cluster == "pbs":
            fout.write("#!/bin/bash\n")
            fout.write("#PBS -N {}\n".format("run_ragtag"))
            fout.write("#PBS -e {}\n".format(os.path.join(sample_output_folder, "run_ragtag.elog")))
            fout.write("#PBS -o {}\n".format(os.path.join(sample_output_folder, "run_ragtag.olog")))
            fout.write("#PBS -q {}\n".format("batch"))
            fout.write("#PBS -l nodes={}:ppn={}\n".format(random.choice(nodes), thread))
            fout.write("##PBS -V\n")
        else:
            fout.write("#!/bin/bash\n")
            fout.write("#SBATCH --job-name={}\n".format("run_ragtag"))
            fout.write("#SBATCH --partition=CPUqueue_fcp1\n")
            fout.write("#SBATCH --account=songbo\n")
            fout.write("#SBATCH --qos=songbo_1_qos\n")
            fout.write("#SBATCH --nodes=1\n")
            fout.write("#SBATCH --ntasks-per-node=1\n")
            fout.write("#SBATCH --cpus-per-task={}\n".format(thread))
            fout.write("#SBATCH --mem={}G\n".format(int(500 / thread)))
            fout.write("#SBATCH -o {}\n".format(os.path.join(sample_output_folder, "run_ragtag.olog")))
            fout.write("#SBATCH -e {}\n".format(os.path.join(sample_output_folder, "run_ragtag.elog")))
            fout.write("#SBATCH -w {}\n".format(random.choice(nodes)))


        fout.write("start_time=$(date +%s)" + "\n")

        fout.write("ragtag.py scaffold -u -o {} -t {} {} {}\n".format(sample_output_folder, thread, ref_path, sample_path))
        fout.write("samtools faidx {}\n".format(os.path.join(sample_output_folder, "ragtag.scaffold.fasta")))

        fout.write("end_time=$(date +%s)" + "\n")
        fout.write("elapsed_seconds=$((end_time - start_time))" + "\n")
        fout.write("elapsed_minutes=$((elapsed_seconds / 60))" + "\n")
        fout.write("echo \"Total $elapsed_seconds Seconds, $elapsed_minutes Minutes\"" + "\n")

    if not os.path.exists(os.path.join(sample_output_folder, "ragtag.scaffold.fasta.fai")):
        if cluster == "pbs":
            os.system("qsub {}".format(script_out_path))
        else:
            os.system("sbatch {}".format(script_out_path))


if __name__ == '__main__':

    options = parse_arguments()

    if not os.path.exists(options.output_folder):
        os.mkdir(options.output_folder)

    for item in os.listdir(options.input_folder):
        if os.path.isdir(os.path.join(options.input_folder, item)):
            for sample in os.listdir(os.path.join(options.input_folder, item)):
                if not (sample.endswith("fasta") or sample.endswith("fasta.gz") or sample.endswith(".fa") or sample.endswith(".fa.gz")):
                    continue

                run_ragtag_for_one_fa(options.ref_path, sample, os.path.join(options.input_folder, "{}/{}".format(item, sample)), options.output_folder, options.thread, options.nodes, options.cluster)

        else:

            if not (item.endswith("fasta") or item.endswith("fasta.gz") or item.endswith(".fa") or item.endswith(".fa.gz")):
                continue

            run_ragtag_for_one_fa(options.ref_path, item, os.path.join(options.input_folder, item), options.output_folder, options.thread, options.nodes, options.cluster)

