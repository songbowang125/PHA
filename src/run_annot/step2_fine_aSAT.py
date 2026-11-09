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


def find_cent_regions(sample_coarse_aSAT_folder):
    avg_peri_cent_length = {'chr1': 20917949, 'chr10': 6377898, 'chr11': 6954504, 'chr12': 4679634, 'chr13': 9645842, 'chr14': 11966278, 'chr15': 12216493, 'chr16': 14059081, 'chr17': 5660880, 'chr18': 5834728, 'chr19': 9264118, 'chr2': 5313196, 'chr20': 7268969, 'chr21': 8972917, 'chr22': 17307563, 'chr3': 5920141, 'chr4': 4920277, 'chr5': 5675134, 'chr6': 4634109, 'chr7': 9618651, 'chr8': 4125243, 'chr9': 25440684, 'chrX': 3910857, 'chrY': 862561}
    chm13v2_cent_region = {'chr1': {'start': 121796048, 'end': 126300487, 'length': 4504439}, 'chr10': {'start': 39633793, 'end': 41664589, 'length': 2030796}, 'chr11': {'start': 51035789, 'end': 54450838, 'length': 3415049}, 'chr12': {'start': 34620838, 'end': 37202490, 'length': 2581652}, 'chr13': {'start': 15547593, 'end': 17498291, 'length': 1950698}, 'chr14': {'start': 10092112, 'end': 12708411, 'length': 2616299}, 'chr15': {'start': 16678794, 'end': 17694466, 'length': 1015672}, 'chr16': {'start': 35848286, 'end': 37829521, 'length': 1981235}, 'chr17': {'start': 23892419, 'end': 27486939, 'length': 3594520}, 'chr18': {'start': 15965699, 'end': 20933550, 'length': 4967851}, 'chr19': {'start': 25817676, 'end': 29768171, 'length': 3950495}, 'chr2': {'start': 92333543, 'end': 94673023, 'length': 2339480}, 'chr20': {'start': 26925852, 'end': 29099655, 'length': 2173803}, 'chr21': {'start': 10962853, 'end': 11306205, 'length': 343352}, 'chr22': {'start': 12788180, 'end': 15711065, 'length': 2922885}, 'chr3': {'start': 91738002, 'end': 96415026, 'length': 4677024}, 'chr4': {'start': 49705154, 'end': 55199795, 'length': 5494641}, 'chr5': {'start': 47039134, 'end': 49596625, 'length': 2557491}, 'chr6': {'start': 58286706, 'end': 61058390, 'length': 2771684}, 'chr7': {'start': 60414372, 'end': 63714499, 'length': 3300127}, 'chr8': {'start': 44215832, 'end': 46325080, 'length': 2109248}, 'chr9': {'start': 44951775, 'end': 47582595, 'length': 2630820}, 'chrX': {'start': 57820107, 'end': 60927025, 'length': 3106918}, 'chrY': {'start': 10565750, 'end': 10883085, 'length': 317335}}
    chm13v2_chrom_length = {'chr1': 248387328, 'chr10': 134758134, 'chr11': 135127769, 'chr12': 133324548, 'chr13': 113566686, 'chr14': 101161492, 'chr15': 99753195, 'chr16': 96330374, 'chr17': 84276897, 'chr18': 80542538, 'chr19': 61707364, 'chr2': 242696752, 'chr20': 66210255, 'chr21': 45090682, 'chr22': 51324926, 'chr3': 201105948, 'chr4': 193574945, 'chr5': 182045439, 'chr6': 172126628, 'chr7': 160567428, 'chr8': 146259331, 'chr9': 150617247, 'chrX': 154259566, 'chrY': 62460029}
    chm13v2_chroms = chm13v2_chrom_length.keys()

    # # STEP: find the core aSAT region using the results of coarse aSAT
    coarse_aSAT_file = os.path.join(sample_coarse_aSAT_folder, "coarse_aSAT.10k.bed")
    longer_than_100k_cord_dict = {}
    with open(coarse_aSAT_file) as fin:
        for line in fin:
            line_split = line.strip().split("\t")
            aSAT_chrom, aSAT_start, aSAT_end, aSAT_length = line_split[0].replace("_RagTag", ""), int(line_split[1]), int(line_split[2]), int(line_split[2]) - int(line_split[1])

            if aSAT_chrom not in chm13v2_chroms:
                continue

            if aSAT_length >= 100000 * 0.8:
                if aSAT_chrom not in longer_than_100k_cord_dict:
                    longer_than_100k_cord_dict[aSAT_chrom] = []

                longer_than_100k_cord_dict[aSAT_chrom].append(aSAT_start)
                longer_than_100k_cord_dict[aSAT_chrom].append(aSAT_end)

    final_core_aSAT_regions = {}
    final_core_cent_regions = {}
    final_peri_cent_regions = {}

    for chrom in longer_than_100k_cord_dict:
        core_aSAT_start, core_aSAT_end = min(longer_than_100k_cord_dict[chrom]), max(longer_than_100k_cord_dict[chrom])

        core_cent_start, core_cent_end = max(1, core_aSAT_start - 1000000), core_aSAT_end + 1000000
        peri_cent_start, peri_cent_end = max(1, core_aSAT_start - avg_peri_cent_length[chrom]), core_aSAT_end + avg_peri_cent_length[chrom]

        final_core_aSAT_regions[chrom] = [str(core_aSAT_start), str(core_aSAT_end), str(core_aSAT_end - core_aSAT_start)]
        final_core_cent_regions[chrom] = [str(core_cent_start), str(core_cent_end), str(core_cent_end - core_cent_start)]
        final_peri_cent_regions[chrom] = [str(peri_cent_start), str(peri_cent_end), str(peri_cent_end - peri_cent_start)]

    return final_core_aSAT_regions, final_core_cent_regions, final_peri_cent_regions

def run_fine_aSAT(sample_folder, code_path, thread, nodes, cluster):

    sample_coarse_aSAT_folder = os.path.join(sample_folder, "res_coarse_aSAT")
    sample_fine_aSAT_folder = os.path.join(sample_folder, "res_fine_aSAT")
    sample_ragtag_fa = os.path.join(sample_folder, "res_ragtag/ragtag.scaffold.fasta")

    if not os.path.exists(sample_fine_aSAT_folder):
        os.mkdir(sample_fine_aSAT_folder)

    # # find core cent region
    final_core_aSAT_regions, final_core_cent_regions, final_peri_cent_regions = find_cent_regions(sample_coarse_aSAT_folder)

    with open(os.path.join(sample_fine_aSAT_folder, "stats.cent_regions.txt"), "w") as fout:
        for chrom in final_core_aSAT_regions:
            fout.write("{}_RagTag\t{}\t{}\t{}\n".format(chrom, "\t".join(final_core_aSAT_regions[chrom]), "\t".join(final_core_cent_regions[chrom]), "\t".join(final_peri_cent_regions[chrom])))

    # # run fine aSAT search
    for chrom in final_core_cent_regions:

        hmm_hor_path = os.path.join(code_path, "hmm_models/{}.hmm".format(chrom))
        hmm_sf_path = os.path.join(code_path, "hmm_models/AS-SFs-hmmer3.0.290621.hmm")
        awk_file = os.path.join(code_path, "hmmertblout2bed.awk")
        overlap_py_file = os.path.join(code_path, "overlap_filter.py")

        if cluster == "pbs":
            script_out_path = os.path.join(sample_fine_aSAT_folder, "run_fine_aSAT.{}.pbs".format(chrom))
        else:
            script_out_path = os.path.join(sample_fine_aSAT_folder, "run_fine_aSAT.{}.slurm".format(chrom))

        with open(script_out_path, "w") as fout:
            if cluster == "pbs":
                fout.write("#!/bin/bash\n")
                fout.write("#PBS -N {}\n".format("run_fine_aSAT.{}.pbs".format(chrom)))
                fout.write("#PBS -e {}\n".format(os.path.join(sample_fine_aSAT_folder, "run_fine_aSAT.{}.elog".format(chrom))))
                fout.write("#PBS -o {}\n".format(os.path.join(sample_fine_aSAT_folder, "run_fine_aSAT.{}.olog".format(chrom))))
                fout.write("#PBS -q {}\n".format("batch"))
                fout.write("#PBS -l nodes={}:ppn={}\n".format(random.choice(nodes), thread))
                fout.write("##PBS -V\n")
            else:
                fout.write("#!/bin/bash\n")
                fout.write("#SBATCH --job-name={}\n".format("run_fine_aSAT.{}.pbs".format(chrom)))
                fout.write("#SBATCH --partition=CPUqueue_fcp1\n")
                fout.write("#SBATCH --account=songbo\n")
                fout.write("#SBATCH --qos=songbo_1_qos\n")
                fout.write("#SBATCH --nodes=1\n")
                fout.write("#SBATCH --ntasks-per-node=1\n")
                fout.write("#SBATCH --cpus-per-task={}\n".format(thread))
                fout.write("#SBATCH --mem={}G\n".format(int(500 / (64 / thread))))
                fout.write("#SBATCH -o {}\n".format(os.path.join(sample_fine_aSAT_folder, "run_fine_aSAT.{}.olog".format(chrom))))
                fout.write("#SBATCH -e {}\n".format(os.path.join(sample_fine_aSAT_folder, "run_fine_aSAT.{}.elog".format(chrom))))
                fout.write("#SBATCH -w {}\n".format(random.choice(nodes)))

            fout.write("start_time=$(date +%s)" + "\n")

            # # fetch seq
            chrom_core_cent_path = os.path.join(sample_fine_aSAT_folder, "{}_RagTag.core_cent.fa".format(chrom))
            chrom_peri_cent_path = os.path.join(sample_fine_aSAT_folder, "{}_RagTag.peri_cent.fa".format(chrom))
            fout.write("samtools faidx {} {}_RagTag:{}-{} > {}\n".format(sample_ragtag_fa, chrom, final_core_cent_regions[chrom][0], final_core_cent_regions[chrom][1], chrom_core_cent_path))
            fout.write("samtools faidx {} {}_RagTag:{}-{} > {}\n".format(sample_ragtag_fa, chrom, final_peri_cent_regions[chrom][0], final_peri_cent_regions[chrom][1], chrom_peri_cent_path))

            # # run nhmmer
            hor_tbl_out = os.path.join(sample_fine_aSAT_folder, "{}_RagTag.fine_aSAT_HOR.tbl".format(chrom))
            sf_tbl_out = os.path.join(sample_fine_aSAT_folder, "{}_RagTag.fine_aSAT_SF.tbl".format(chrom))
            fout.write("nhmmer -o /dev/null --notextw --noali --dna --cpu {} --tblout {} {} {}\n".format(thread, hor_tbl_out, hmm_hor_path, chrom_core_cent_path))
            fout.write("nhmmer -o /dev/null --notextw --noali --dna --cpu {} --tblout {} {} {}\n".format(thread, sf_tbl_out, hmm_sf_path, chrom_core_cent_path))

            # # post process of nhmmer res
            for tbl_out in [hor_tbl_out, sf_tbl_out]:
                tbl_out_awk = tbl_out + ".awk"
                tbl_out_awk_sort = tbl_out + ".awk" + ".sort"
                tbl_out_awk_sort_bedmap = tbl_out + ".awk" + ".sort" + ".bedmap"
                tbl_out_awk_sort_bedmap_awk = tbl_out + ".awk" + ".sort" + ".bedmap" + ".awk"
                tbl_out_awk_sort_bedmap_awk_overlap = tbl_out + ".awk" + ".sort" + ".bedmap" + ".awk" + ".overlap"
                tbl_out_final_bed = tbl_out + ".final.bed"

                fout.write("awk -v th=0.7 -f {} {} > {}\n".format(awk_file, tbl_out, tbl_out_awk))
                fout.write("sort -k 1.4,1 -k 2,2n {} > {}\n".format(tbl_out_awk, tbl_out_awk_sort))
                fout.write("bedmap --max-element --fraction-either 0.1 {} > {}\n".format(tbl_out_awk_sort, tbl_out_awk_sort_bedmap))
                fout.write("awk \"{if(!(\$0 in a)){a[\$0]; print}}\" " + tbl_out_awk_sort_bedmap + " > " + tbl_out_awk_sort_bedmap_awk + "\n")
                fout.write("python {} {} > {}\n".format(overlap_py_file, tbl_out_awk_sort_bedmap_awk, tbl_out_awk_sort_bedmap_awk_overlap))

                fout.write("rm {} {} {}\n".format(tbl_out_awk_sort, tbl_out_awk_sort_bedmap, tbl_out_awk_sort_bedmap_awk))
                fout.write("mv {} {}\n".format(tbl_out_awk_sort_bedmap_awk_overlap, tbl_out_final_bed))

            fout.write("end_time=$(date +%s)" + "\n")
            fout.write("elapsed_seconds=$((end_time - start_time))" + "\n")
            fout.write("elapsed_minutes=$((elapsed_seconds / 60))" + "\n")
            fout.write("echo \"Total $elapsed_seconds Seconds, $elapsed_minutes Minutes\"" + "\n")

        if not (os.path.exists(hor_tbl_out + ".final.bed") and os.path.exists(sf_tbl_out + ".final.bed")):
            if cluster == "pbs":
                os.system("qsub {}".format(script_out_path))
            else:
                os.system("sbatch {}".format(script_out_path))

if __name__ == '__main__':

    options = parse_arguments()

    cnt = 0
    for sample_folder in os.listdir(options.input_folder):

        cnt += 1
        sample_folder = os.path.join(options.input_folder, sample_folder)

        run_fine_aSAT(sample_folder,  options.code_path, options.thread, options.nodes, options.cluster)

        print(sample_folder)
