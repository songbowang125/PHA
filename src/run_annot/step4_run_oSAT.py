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


def run_other_SAT(sample_folder, code_path, thread, nodes, cluster):

    sample_fine_aSAT_folder = os.path.join(sample_folder, "res_fine_aSAT")
    sample_other_SAT_folder = os.path.join(sample_folder, "res_other_SAT")

    if not os.path.exists(sample_other_SAT_folder):
        os.mkdir(sample_other_SAT_folder)

    for peri_cent_fa in os.listdir(sample_fine_aSAT_folder):
        if not peri_cent_fa.endswith("_RagTag.peri_cent.fa"):
            continue

        chrom = peri_cent_fa.split("_")[0]

        chrom_other_SAT_folder = os.path.join(sample_other_SAT_folder, chrom)
        if not os.path.exists(chrom_other_SAT_folder):
            os.mkdir(chrom_other_SAT_folder)

        peri_cent_fa = os.path.join(sample_fine_aSAT_folder, peri_cent_fa)

        if cluster == "pbs":
            script_out_path = os.path.join(chrom_other_SAT_folder, "run_other_SAT.{}.pbs".format(chrom))
        else:
            script_out_path = os.path.join(chrom_other_SAT_folder, "run_other_SAT.{}.slurm".format(chrom))

        finish_out = os.path.join(chrom_other_SAT_folder, "run_other_SAT.finish_{}.log".format(chrom))

        with open(script_out_path, "w") as fout:
            if cluster == "pbs":
                fout.write("#!/bin/bash\n")
                fout.write("#PBS -N {}\n".format("run_other_SAT.{}.pbs".format(chrom)))
                fout.write("#PBS -e {}\n".format(os.path.join(chrom_other_SAT_folder, "run_other_SAT.{}.elog".format(chrom))))
                fout.write("#PBS -o {}\n".format(os.path.join(chrom_other_SAT_folder, "run_other_SAT.{}.olog".format(chrom))))
                fout.write("#PBS -q {}\n".format("batch"))
                fout.write("#PBS -l nodes={}:ppn={}\n".format(random.choice(nodes), thread))
                fout.write("##PBS -V\n")
            else:
                fout.write("#!/bin/bash\n")
                fout.write("#SBATCH --job-name={}\n".format("run_other_SAT.{}.pbs".format(chrom)))
                fout.write("#SBATCH --partition=CPUqueue_fcp1\n")
                fout.write("#SBATCH --account=songbo\n")
                fout.write("#SBATCH --qos=songbo_1_qos\n")
                fout.write("#SBATCH --nodes=1\n")
                fout.write("#SBATCH --ntasks-per-node=1\n")
                fout.write("#SBATCH --cpus-per-task={}\n".format(thread))
                fout.write("#SBATCH --mem={}G\n".format(int(500 / (64 / thread))))
                fout.write("#SBATCH -o {}\n".format(os.path.join(chrom_other_SAT_folder, "run_other_SAT.{}.olog".format(chrom))))
                fout.write("#SBATCH -e {}\n".format(os.path.join(chrom_other_SAT_folder, "run_other_SAT.{}.elog".format(chrom))))
                fout.write("#SBATCH -w {}\n".format(random.choice(nodes)))

            fout.write("start_time=$(date +%s)" + "\n")

            fout.write("cd {} \n".format(chrom_other_SAT_folder))

            fout.write("cp {}/HSat23/Assembly_HSat2and3_v3.pl {}/HSat23/HSat2_kmers.txt {}/HSat23/HSat3_kmers.txt {}\n".format(code_path, code_path, code_path, chrom_other_SAT_folder))

            fout.write("perl {}/Assembly_HSat2and3_v3.pl {} {} > /dev/null \n".format(chrom_other_SAT_folder, peri_cent_fa, chrom_other_SAT_folder))
            fout.write("rm {}/Assembly_HSat2and3_v3.pl {}/HSat2_kmers.txt {}/HSat3_kmers.txt \n".format(chrom_other_SAT_folder, chrom_other_SAT_folder, chrom_other_SAT_folder))

            fout.write("RepeatMasker -engine rmblast -species human -dir {} -s -pa {} {} > /dev/null \n".format(chrom_other_SAT_folder, int(thread / 4), peri_cent_fa))

            fout.write("end_time=$(date +%s)" + "\n")
            fout.write("elapsed_seconds=$((end_time - start_time))" + "\n")
            fout.write("elapsed_minutes=$((elapsed_seconds / 60))" + "\n")

            fout.write("echo \"Total $elapsed_seconds Seconds, $elapsed_minutes Minutes\"" + "\n")

            fout.write("echo done > {}".format(finish_out) + "\n")

        if not os.path.exists(finish_out):
            if cluster == "pbs":
                os.system("qsub {}".format(script_out_path))
            else:
                os.system("sbatch {}".format(script_out_path))

        # print(peri_cent_fa)
        # break

if __name__ == '__main__':

    options = parse_arguments()

    cnt = 0
    for sample_folder in os.listdir(options.input_folder):

        cnt += 1
        sample_folder = os.path.join(options.input_folder, sample_folder)

        run_other_SAT(sample_folder,  options.code_path, options.thread, options.nodes, options.cluster)

        print(sample_folder)

        # break
