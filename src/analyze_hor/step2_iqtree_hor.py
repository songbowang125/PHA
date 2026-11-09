import os
import pandas as pd
from collections import Counter



def load_hor(hor_bed):

    hor_list = []

    with open(hor_bed) as fin:
        for line in fin:

            line_split = line.strip().split()

            chrom, hor, orient = line_split[0], line_split[3], line_split[4]

            hor_list.append("{}".format(hor))

    return hor_list


def generate_hor_symbol_dict(hor_list):

    available_ids = "ABCDEFGHIJKLMNOPQRSTUV0123456789"
    available_ids_cnt = len(available_ids)

    hor_symbol_dict = {}
    symbol_hor_dict = {}

    Counter(hor_list).most_common()

    id = 0
    for hor, count in Counter(hor_list).most_common():

        if id >= available_ids_cnt - 1:
            id = available_ids_cnt - 1


        hor_symbol_dict[hor] = available_ids[id]
        symbol_hor_dict[available_ids[id]] = hor

        id += 1


    return hor_symbol_dict, symbol_hor_dict

def generate_data_and_run_mafft_iqtree(sample_hor_dict, output_folder, spec_chrom, sample_cluster_dict, add_ape_as_root, thread):

    all_hor_list = []
    for sample in sample_hor_dict:
        all_hor_list.extend(sample_hor_dict[sample])

    hor_symbol_dict, symbol_hor_dict = generate_hor_symbol_dict(all_hor_list)

    output_fa_path = os.path.join(output_folder, "{}.hor.fa".format(spec_chrom))
    output_mafft_path = os.path.join(output_folder, "{}.hor.mafft.fa".format(spec_chrom))

    with open(output_fa_path, "w") as fout:
        for sample in sample_hor_dict:
            fout.write(">{}_{}\n".format(sample_cluster_dict[sample], sample))
            fout.write("{}\n".format("".join([hor_symbol_dict[hor] for hor in sample_hor_dict[sample]])))


    with open(os.path.join(output_folder, "{}.hor_dict.txt".format(spec_chrom)), "w") as fout:
        for hor in hor_symbol_dict:
            fout.write("{}\t{}\t{}\n".format(hor, hor_symbol_dict[hor], all_hor_list.count(hor)))

    # # run mafft
    if run_machine == "slurm":

        script_file = os.path.join(output_folder, "{}.run_iqtree.slurm".format(spec_chrom))
        with open(script_file, "w") as fout:
            fout.write("#!/bin/bash\n")
            fout.write("#SBATCH --job-name={}\n".format("{}.run_iqtree".format(spec_chrom)))
            fout.write("#SBATCH --account=songbo\n")
            # fout.write("#SBATCH --partition=CPUqueue_fcp1\n")
            # fout.write("#SBATCH --qos=songbo_1_qos\n")
            fout.write("#SBATCH --partition=nodes_fcp1\n")
            fout.write("#SBATCH --qos=songbo_4_qos\n")

            fout.write("#SBATCH --nodes=1\n")
            fout.write("#SBATCH --ntasks-per-node=1\n")
            fout.write("#SBATCH --cpus-per-task={}\n".format(thread))
            fout.write("#SBATCH --mem={}G\n".format(memory))
            fout.write("#SBATCH -o {}\n".format(os.path.join(output_folder, "{}.run_iqtree.olog".format(spec_chrom))))
            fout.write("#SBATCH -e {}\n".format(os.path.join(output_folder, "{}.run_iqtree.elog".format(spec_chrom))))

            fout.write("mafft --thread {} --text {} > {}\n".format(thread, output_fa_path, output_mafft_path))

            if add_ape_as_root:
                fout.write("/home/data/user/songbo/workspace/tools/iqtree-3.1.3-Linux/bin/iqtree3 -nt {} --seqtype MORPH -B 1000 -s {} -o 8_GCA_028858775.3_NHGRI_mPanTro3-v2.1_pri_genomic.fna.gz\n".format(thread, output_mafft_path))
            else:
                fout.write("/home/data/user/songbo/workspace/tools/iqtree-3.1.3-Linux/bin/iqtree3 -nt {} --seqtype MORPH -B 1000 -s {}\n".format(thread, output_mafft_path))

        os.system("sbatch {}".format(script_file))

    else:
        os.system("mafft --thread {} --text {} > {}".format(thread, output_fa_path, output_mafft_path))

        if add_ape_as_root:
            os.system("/mnt/d/tools/iqtree-3.1.3-Linux/bin/iqtree3 -nt {} --seqtype MORPH -B 1000 -s {} -o 8_GCA_028858775.3_NHGRI_mPanTro3-v2.1_pri_genomic.fna.gz".format(thread, output_mafft_path))
            # os.system("/mnt/d/tools/iqtree-3.1.3-Linux/bin/iqtree3 -nt {} --seqtype MORPH -B 1000 -s {} -o Out_mPanTro3_1".format(thread, output_mafft_path))
        else:
            os.system("/mnt/d/tools/iqtree-3.1.3-Linux/bin/iqtree3 -nt {} --seqtype MORPH -B 1000 -s {}".format(thread, output_mafft_path))


def select_sample(feature_file, select_num):

    df = pd.read_csv(feature_file, sep="\t")

    df = df[~df[spec_cluster].astype(str).str.contains("-1", regex=False)].copy()

    if spec_cluster in ["cluster_3", "cluster_4"]:
        df = df[df[spec_cluster].astype(str).str.contains("Q3_Q1", regex=False)].copy()

    if cpc2_sample_only:
        df = df[df["Sample"].astype(str).str.contains("CPCERZ", regex=False)].copy()

    df = df[df["Sample"] != "T2T-CHM13.fa.gz"].copy()

    global_mean = df["aSat_length"].mean()

    selected_list = []
    selected_samples = []
    sample_cluster_dict = {}

    for cluster, group in df.groupby(spec_cluster):
        if select_num == "all":
            selected = group

        else:
            if select_mode == "random":
                # # random select
                selected = group.sample(n=min(select_num, len(group)), random_state=10)
            elif select_mode == "length":
                # # select by cluster mean
                cluster_mean = group["aSat_length"].mean()
                group["distance_to_cluster_mean"] = (group["aSat_length"] - cluster_mean).abs()
                selected = group.nsmallest(select_num, "distance_to_cluster_mean")
            else:
                print("No this select mode", select_mode)
                exit()

        for _, row in selected.iterrows():
            selected_list.append([row["Sample"], row["detail_popu"],  row["aSat_length"], row[spec_cluster]])
            selected_samples.append(row["Sample"])

            sample_cluster_dict[row["Sample"]] = row[spec_cluster]

    return selected_samples, sample_cluster_dict




if __name__ == '__main__':


    add_ape_as_root = False
    cpc2_sample_only = False

    sample_num_each_cluster = "all"

    select_mode = "random"

    run_machine = "slurm"

    spec_cluster = "cluster_3"

    if run_machine == "slurm":
        thread = 48
        memory = 375

        input_folder = "/home/data/user/songbo/workspace/cpc2_centro/res_all_annot"

        cluster_folder = "/home/data/user/songbo/workspace/cpc2_centro/principal_hor_analysis_outlier"

    else:
        thread = 24

        input_folder = "/mnt/e/workspace/cpc2_centro/res_all_annot"

        cluster_folder = "/mnt/e/workspace/cpc2_centro/res_all_post/principal_hor_analysis_outlier"

    if cpc2_sample_only:

        output_folder = os.path.join(cluster_folder, "{}-hor-{}-n{}-cpc".format(spec_cluster, select_mode, sample_num_each_cluster))
    else:
        output_folder = os.path.join(cluster_folder, "{}-hor-{}-n{}".format(spec_cluster, select_mode, sample_num_each_cluster))

    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    chrom_order = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr16", "chr17", "chr19", "chr20", "chr21"]


    chrom_inv_sample_dict = {"chr1": [],
                             "chr2": ["CPCERZ200144.2.fa.gz", "CPCERZ200272.2.fa.gz", "CPCERZ200294.1.fa.gz", "CPCERZ200353.1.fa.gz", "CPCERZ200359.2.fa.gz", "CPCERZ200369.1.fa.gz"],
                             "chr3": ["NA19043_hap1_hprc_r2_v1.0.1.fa.gz"],
                             "chr4": [], "chr5": [], "chr6": [],
                             "chr7": ["CPCERZ200318.2.fa.gz", "NA18879_hap1_hprc_r2_v1.0.1.fa.gz", "NA19185_hap1_hprc_r2_v1.0.1.fa.gz", "NA19338_hap2_hprc_r2_v1.0.1.fa.gz", "NA20355.vrk-ps-sseq.asm-hap2.fasta", "NA21309_pat_hprc_r2_v1.0.1.fa.gz"],
                             "chr8": [],
                             "chr9": ["CPCERZ200139.2.fa.gz", "CPCERZ200179.1.fa.gz", "CPCERZ200185.1.fa.gz", "CPCERZ200217.2.fa.gz", "CPCERZ200267.1.fa.gz", "CPCERZ200332.1.fa.gz", "HG00864.vrk-ps-sseq.asm-hap2.fasta", "HG01890.vrk-ps-sseq.asm-hap2.fasta", "HG02011.vrk-ps-sseq.asm-hap1.fasta", "HG02630_pat_hprc_r2_v1.0.1.fa.gz"],
                             "chr10": [], "chr11": [], "chr12": [], "chr13": [], "chr14": ["CPCERZ200146.2.fa.gz", "CPCERZ200156.2.fa.gz", "CPCERZ200256.2.fa.gz", "CPCERZ200295.2.fa.gz"],
                             "chr15": ["CPCERZ200206.1.fa.gz", "HG00235_hap2_hprc_r2_v1.0.1.fa.gz"], "chr16": ["CPCERZ200258.1.fa.gz"],
                             "chr17": [], "chr18": [], "chr19": [], "chr20": [], "chr21": [],
                             "chr22": ["CPCERZ200082.1.fa.gz", "CPCERZ200453.1.fa.gz", "CPCERZ200458.2.fa.gz", "HG02282.vrk-ps-sseq.asm-hap2.fasta", "HG02953.vrk-ps-sseq.asm-hap1.fasta", "NA18944_hap2_hprc_r2_v1.0.1.fa.gz", "NA19159_hap1_hprc_r2_v1.0.1.fa.gz", "NA19347.vrk-ps-sseq.asm-hap2.fasta"],
                             "chrX": ["CPCERZ200045.1.fa.gz"], "chrY": ["CPCERZ200084.1.fa.gz", "CPCERZ200233.1.fa.gz", "CPCERZ200271.1.fa.gz", "CPCERZ200353.1.fa.gz", "CPCERZ200401.1.fa.gz", "CPCERZ200434.2.fa.gz"]}


    for spec_chrom in chrom_order:
        os.system("rm -f {}/{}.*".format(output_folder, spec_chrom))

        if spec_cluster == "cluster_1":
            selected_samples, sample_cluster_dict = select_sample(os.path.join(cluster_folder, "{}.features.c1.txt".format(spec_chrom)), sample_num_each_cluster)
        elif spec_cluster == "cluster_2":
            selected_samples, sample_cluster_dict = select_sample(os.path.join(cluster_folder, "{}.features.c2.txt".format(spec_chrom)), sample_num_each_cluster)
        elif spec_cluster in ["cluster_3", "cluster_4"]:
            selected_samples, sample_cluster_dict = select_sample(os.path.join(cluster_folder, "{}.features.c3_c4_c5_c6.txt".format(spec_chrom)), sample_num_each_cluster)
        else:
            print("No this spec cluster", spec_cluster)
            exit()

        if add_ape_as_root:
            selected_samples.extend(["GCA_028858775.3_NHGRI_mPanTro3-v2.1_pri_genomic.fna.gz"])
            sample_cluster_dict["GCA_028858775.3_NHGRI_mPanTro3-v2.1_pri_genomic.fna.gz"] = 8

        chrom_exclude_samples = chrom_inv_sample_dict[spec_chrom]

        chrom_input_folder = os.path.join(input_folder, spec_chrom)

        sample_hor_dict = {}
        for hor_file in os.listdir(chrom_input_folder):
            if not hor_file.endswith("#stv_hor.bed"):
                continue

            sample_name = hor_file.replace("#stv_hor.bed", "")

            if sample_name not in selected_samples:
                continue

            if sample_name not in chrom_exclude_samples:
                sample_hor_dict[sample_name] = load_hor(os.path.join(chrom_input_folder, hor_file))

                print(spec_chrom, sample_name)

            else:
                chrom_exclude_samples.remove(sample_name)

        # if len(chrom_exclude_samples) != 0:
        #     print("ERROR: there is still exclude samples", chrom_exclude_samples)
        #     exit()

        generate_data_and_run_mafft_iqtree(sample_hor_dict, output_folder, spec_chrom, sample_cluster_dict, add_ape_as_root, thread)


