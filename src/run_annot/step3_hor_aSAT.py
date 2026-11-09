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


def run_stv_hor_aSAT(sample_folder, code_path, thread, nodes, cluster):

    sample_fine_aSAT_folder = os.path.join(sample_folder, "res_fine_aSAT")

    if cluster == "pbs":
        script_out_path = os.path.join(sample_fine_aSAT_folder, "run_hor_aSAT.pbs")
    else:
        script_out_path = os.path.join(sample_fine_aSAT_folder, "run_hor_aSAT.slurm")

    with open(script_out_path, "w") as fout:
        if cluster == "pbs":
            fout.write("#!/bin/bash\n")
            fout.write("#PBS -N {}\n".format("run_hor_aSAT.pbs"))
            fout.write("#PBS -e {}\n".format(os.path.join(sample_fine_aSAT_folder, "run_hor_aSAT.elog")))
            fout.write("#PBS -o {}\n".format(os.path.join(sample_fine_aSAT_folder, "run_hor_aSAT.olog")))
            fout.write("#PBS -q {}\n".format("batch"))
            fout.write("#PBS -l nodes={}:ppn={}\n".format(random.choice(nodes), thread))
            fout.write("##PBS -V\n")
        else:
            fout.write("#!/bin/bash\n")
            fout.write("#SBATCH --job-name={}\n".format("run_hor_aSAT.pbs"))
            fout.write("#SBATCH --partition=CPUqueue_fcp1\n")
            fout.write("#SBATCH --account=songbo\n")
            fout.write("#SBATCH --qos=songbo_1_qos\n")
            fout.write("#SBATCH --nodes=1\n")
            fout.write("#SBATCH --ntasks-per-node=1\n")
            fout.write("#SBATCH --cpus-per-task={}\n".format(thread))
            fout.write("#SBATCH --mem={}G\n".format(int(500 / (64 / thread))))
            fout.write("#SBATCH -o {}\n".format(os.path.join(sample_fine_aSAT_folder, "run_hor_aSAT.olog")))
            fout.write("#SBATCH -e {}\n".format(os.path.join(sample_fine_aSAT_folder, "run_hor_aSAT.elog")))
            fout.write("#SBATCH -w {}\n".format(random.choice(nodes)))

        fout.write("start_time=$(date +%s)" + "\n")

        aSAT_monomer_path = os.path.join(sample_folder, "aSAT.monomer.bed")
        aSAT_hor_live_hor_bed = os.path.join(sample_folder, "aSAT.hor.AS_liveHORs.bed")
        aSAT_hor_stv_raw_bed = os.path.join(sample_folder, "aSAT.hor.stv_raw.bed")
        aSAT_hor_stv_colored_bed = os.path.join(sample_folder, "aSAT.hor.stv_colored.bed")
        aSAT_hor_stv_bed = os.path.join(sample_folder, "aSAT.hor.stv.bed")
        aSAT_hor_stv_stats_tsv = os.path.join(sample_folder, "aSAT.hor.stv_stats.tsv")

        fout.write("cat {}/chr*_RagTag.fine_aSAT_HOR.tbl.final.bed | bedtools sort > {} \n".format(sample_fine_aSAT_folder, aSAT_monomer_path) + "\n")

        # # run stv
        stv_script_path = os.path.join(code_path, "stv/scripts")

        fout.write("python {}/live_HORs_filter.py {} > {} \n".format(stv_script_path, aSAT_monomer_path, aSAT_hor_live_hor_bed))
        fout.write("python {}/mon2stv.py {} > {} \n".format(stv_script_path, aSAT_hor_live_hor_bed, aSAT_hor_stv_raw_bed))
        fout.write("python {}/bed2stat.py {} > {} \n".format(stv_script_path, aSAT_hor_stv_raw_bed, aSAT_hor_stv_stats_tsv))
        fout.write("python {}/stats2table.py {} \n".format(stv_script_path, aSAT_hor_stv_stats_tsv))
        fout.write("python {}/coloring.py {} > {} \n".format(stv_script_path, aSAT_hor_stv_raw_bed, aSAT_hor_stv_colored_bed))
        fout.write("python {}/numbering.py {} > {} \n".format(stv_script_path, aSAT_hor_stv_colored_bed, aSAT_hor_stv_bed))

        fout.write("rm {} {}\n".format(aSAT_hor_stv_colored_bed, aSAT_hor_live_hor_bed))

        # # deal with SF
        aSAT_sf_path = os.path.join(sample_folder, "aSAT.sf.bed")
        fout.write("cat {}/chr*_RagTag.fine_aSAT_SF.tbl.final.bed | bedtools sort > {} \n".format(sample_fine_aSAT_folder, aSAT_sf_path) + "\n")

        fout.write("end_time=$(date +%s)" + "\n")
        fout.write("elapsed_seconds=$((end_time - start_time))" + "\n")
        fout.write("elapsed_minutes=$((elapsed_seconds / 60))" + "\n")
        fout.write("echo \"Total $elapsed_seconds Seconds, $elapsed_minutes Minutes\"" + "\n")

    if cluster == "pbs":
        os.system("qsub {}".format(script_out_path))
    else:
        os.system("sbatch {}".format(script_out_path))

if __name__ == '__main__':

    options = parse_arguments()

    for sample_folder in os.listdir(options.input_folder):


        sample_folder = os.path.join(options.input_folder, sample_folder)

        run_stv_hor_aSAT(sample_folder,  options.code_path, options.thread, options.nodes, options.cluster)

        print(sample_folder)

