import os

import numpy as np
import pandas as pd
import pysam


def generate_data_and_run_mafft_iqtree(sample_seq_dict, output_folder, spec_chrom, sample_cluster_dict, add_ape_as_root, thread):


    output_fa_path = os.path.join(output_folder, "{}.flank.fa".format(spec_chrom))
    output_date_path = os.path.join(output_folder, "{}.flank.date.txt".format(spec_chrom))

    output_mafft_path = os.path.join(output_folder, "{}.flank.mafft.fa".format(spec_chrom))

    with open(output_fa_path, "w") as fout:

        for sample in sample_seq_dict:
            fout.write(">{}_{}\n".format(sample_cluster_dict[sample], sample))
            fout.write("{}\n".format(sample_seq_dict[sample]))

    with open(output_date_path, "w") as fout_date:
        for sample in sample_seq_dict:
            if sample != "GCA_028858775.3_NHGRI_mPanTro3-v2.1_pri_genomic.fna.gz":
                fout_date.write("8_GCA_028858775.3_NHGRI_mPanTro3-v2.1_pri_genomic.fna.gz,{}_{} -6.2\n".format(sample_cluster_dict[sample], sample))

                break

    # # run mafft
    if run_machine == "slurm":
        script_file = os.path.join(output_folder, "{}.flank_iqtree.slurm".format(spec_chrom))
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

            fout.write("mafft --thread {} --auto {} > {}\n".format(thread, output_fa_path, output_mafft_path))

            if add_ape_as_root:

                fout.write("/home/data/user/songbo/workspace/tools/iqtree-3.1.3-Linux/bin/iqtree3 -nt {} -B 1000 -s {} -o 8_GCA_028858775.3_NHGRI_mPanTro3-v2.1_pri_genomic.fna.gz\n".format(thread, output_mafft_path))
                fout.write("/home/data/user/songbo/workspace/tools/iqtree-3.1.3-Linux/bin/iqtree3 -nt {} -s {} -keep-ident -te {}.treefile --date {} --date-tip 0 --prefix {}.treefile.lsd -o 8_GCA_028858775.3_NHGRI_mPanTro3-v2.1_pri_genomic.fna.gz\n".format(thread, output_mafft_path, output_mafft_path, output_date_path, output_mafft_path))

            else:
                fout.write("/home/data/user/songbo/workspace/tools/iqtree-3.1.3-Linux/bin/iqtree3 -nt {} -B 1000 -s {}\n".format(thread, output_mafft_path))
        os.system("sbatch {}".format(script_file))

    else:
        os.system("mafft --thread {} --auto {} > {}".format(thread, output_fa_path, output_mafft_path))

        if add_ape_as_root:

            os.system("/mnt/d/tools/iqtree-3.1.3-Linux/bin/iqtree3 -nt {} -B 1000 -s {} -o 8_GCA_028858775.3_NHGRI_mPanTro3-v2.1_pri_genomic.fna.gz".format(thread, output_mafft_path))
            os.system("/mnt/d/tools/iqtree-3.1.3-Linux/bin/iqtree3 -nt {} -s {} -keep-ident -te {}.treefile --date {} --date-tip 0 --prefix {}.treefile.lsd -o 8_GCA_028858775.3_NHGRI_mPanTro3-v2.1_pri_genomic.fna.gz".format(thread, output_mafft_path, output_mafft_path, output_date_path, output_mafft_path))

        else:
            os.system("/mnt/d/tools/iqtree-3.1.3-Linux/bin/iqtree3 -nt {} -B 1000 -s {}".format(thread, output_mafft_path))


def select_sample(feature_file, select_num):

    df = pd.read_csv(feature_file, sep="\t")

    df = df[~df[spec_cluster].astype(str).str.contains("-1", regex=False)].copy()

    if spec_cluster in ["cluster_3", "cluster_4"]:
        df = df[df[spec_cluster].astype(str).str.contains("Q3_Q1", regex=False)].copy()

    if cpc2_sample_only:
        df = df[df["Sample"].astype(str).str.contains("CPCERZ", regex=False)].copy()

    # df = df[df[spec_cluster].astype(str).str.contains("2_", regex=False)].copy()

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

    add_ape_as_root = True

    cpc2_sample_only = False

    sample_num_each_cluster = "all"

    select_mode = "random"

    run_machine = "slurm"

    arm = "p"

    flank = "50k"

    spec_cluster = "cluster_2"

    if run_machine == "slurm":
        thread = 48
        memory = 375

        if flank == "all":
            input_folder = "/home/data/user/songbo/workspace/cpc2_centro/res_all_flanking/{}-arm".format(arm)
        else:
            input_folder = "/home/data/user/songbo/workspace/cpc2_centro/res_all_flanking_{}/{}-arm".format(flank, arm)

        cluster_folder = "/home/data/user/songbo/workspace/cpc2_centro/principal_hor_analysis_outlier"
        work_folder = "/home/data/user/songbo/workspace/cpc2_centro/"

    else:
        thread = 12

        if flank == "all":
            input_folder = "/mnt/g/res_all_flanking/{}-arm".format(arm)
        else:
            input_folder = "/mnt/e/workspace/cpc2_centro/res_all_flanking_{}/{}-arm".format(flank, arm)

        cluster_folder = "/mnt/e/workspace/cpc2_centro/res_all_post/principal_hor_analysis_outlier"
        work_folder = "/mnt/e/workspace/cpc2_centro/res_all_post/"

    if add_ape_as_root:
        if cpc2_sample_only:
            output_folder = os.path.join(cluster_folder, "{}-flank_{}_{}-{}-n{}-cpc-ape".format(spec_cluster, arm, flank, select_mode, sample_num_each_cluster))
        else:
            output_folder = os.path.join(cluster_folder, "{}-flank_{}_{}-{}-n{}-ape".format(spec_cluster, arm, flank, select_mode, sample_num_each_cluster))

    else:
        if cpc2_sample_only:
            output_folder = os.path.join(cluster_folder, "{}-flank_{}_{}-{}-n{}-cpc".format(spec_cluster, arm, flank, select_mode, sample_num_each_cluster))
        else:
            output_folder = os.path.join(cluster_folder, "{}-flank_{}_{}-{}-n{}".format(spec_cluster, arm, flank, select_mode, sample_num_each_cluster))

    if not os.path.exists(output_folder):
        os.mkdir(output_folder)


    # chrom_order = ["chr1", "chr2", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr16", "chr17", "chr19", "chr20", "chr21"]
    # chrom_order = ["chr1", "chr19", "chr5", "chr21", "chr13", "chr10", "chr20", "chr16", "chr17",  "chr2", "chr4", "chr6", "chr7", "chr8", "chr9", "chr11", "chr12", ]
    chrom_order = ["chr21", "chr13", "chr19", "chr5", "chr1", "chr10", "chr20", "chr16", "chr17",  "chr2", "chr4", "chr6", "chr7", "chr8", "chr9", "chr11", "chr12"]

    #
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

        chrom_flanking_stats_file = os.path.join(work_folder, "flanking_{}_{}_{}/{}.max_align.stats.txt".format(arm, flank, spec_cluster, spec_chrom))
        print(chrom_flanking_stats_file)
        cluster_flanking_stats_dict = {}
        with open(chrom_flanking_stats_file) as fin:
            for line in fin:
                if "median" in line:
                    continue
                line_split = line.strip().split()

                # _, cluster, min, max, median, q1, q3, iqr, lower_bound, upper_bound, _, _ = line.strip().split()

                cluster = line_split[1]
                q1 = line_split[5]
                q3 = line_split[6]
                # lower = line_split[8]
                # upper = line_split[9]

                cluster_flanking_stats_dict[cluster] = [float(q1), float(q3)]
                # cluster_flanking_stats_dict[cluster] = [float(lower), float(upper)]


        os.system("rm -f {}/{}.*".format(output_folder, spec_chrom))

        if spec_cluster == "cluster_1":
            selected_samples, sample_cluster_dict = select_sample(os.path.join(cluster_folder, "{}.features.c1.txt".format(spec_chrom)), sample_num_each_cluster)
        elif spec_cluster == "cluster_2":
            selected_samples, sample_cluster_dict = select_sample(os.path.join(cluster_folder, "{}.features.c2.txt".format(spec_chrom)), sample_num_each_cluster)
        elif spec_cluster in ["cluster_3", "cluster_4"]:
            selected_samples, sample_cluster_dict = select_sample(os.path.join(cluster_folder, "{}.features.c3_c4_c5_c6.txt".format(spec_chrom)), sample_num_each_cluster)

        if add_ape_as_root:
            selected_samples.extend(["GCA_028858775.3_NHGRI_mPanTro3-v2.1_pri_genomic.fna.gz"])
            sample_cluster_dict["GCA_028858775.3_NHGRI_mPanTro3-v2.1_pri_genomic.fna.gz"] = 8

        chrom_exclude_samples = chrom_inv_sample_dict[spec_chrom]

        chrom_input_folder = os.path.join(input_folder, spec_chrom)

        sample_seq_dict = {}
        for sample_fasta in os.listdir(chrom_input_folder):
            if not sample_fasta.endswith(".flanking.fa"):
                continue

            sample_name = sample_fasta.replace(".flanking.fa", "")

            if sample_name not in selected_samples:
                continue

            sample_cluster = str(sample_cluster_dict[sample_name])

            # # filter the low mapped sample
            if sample_name != "GCA_028858775.3_NHGRI_mPanTro3-v2.1_pri_genomic.fna.gz":
                sample_paf = os.path.join(chrom_input_folder, "{}.paf".format(sample_name))
                max_align_ratio = 0
                with open(sample_paf) as fin:
                    for line in fin:
                        line_split = line.strip().split()

                        target_length, align_length = int(line_split[6]), int(line_split[9])

                        align_ratio = align_length / target_length

                        if align_ratio > max_align_ratio:
                            max_align_ratio = align_ratio
                            # print(align_length, target_length, align_ratio)

                # if max_align_ratio < 0.95:
                #     continue
                if not (cluster_flanking_stats_dict[sample_cluster][0] <= max_align_ratio <= cluster_flanking_stats_dict[sample_cluster][1]):
                    continue

            # # filter the large inv samples
            if sample_name not in chrom_exclude_samples:
                with pysam.FastaFile(os.path.join(chrom_input_folder, sample_fasta)) as fin:
                    sample_seq_dict[sample_name] = fin.fetch(fin.references[0])
            else:
                chrom_exclude_samples.remove(sample_name)

        generate_data_and_run_mafft_iqtree(sample_seq_dict, output_folder, spec_chrom, sample_cluster_dict, add_ape_as_root, thread)
