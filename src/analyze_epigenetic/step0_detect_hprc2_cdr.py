import gzip
import os
import re
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import ranksums
from scipy.stats import mannwhitneyu
from scipy.stats import ks_2samp
from scipy.stats import mannwhitneyu
from scipy.stats import permutation_test
import multiprocessing
import json
from scipy.stats import kruskal

plt.rcParams["svg.fonttype"] = "none"


def reverse_hor(in_hor):
    name = in_hor.split(".")[0]
    combine = in_hor.split(".")[1]

    # # reverse combine
    ids = re.split(r'[-_]', combine)
    seps = re.findall(r'[-_]', combine)

    ids.reverse()
    seps.reverse()

    return name + "." + ids[0] + ''.join(sep + id_ for sep, id_ in zip(seps, ids[1:]))


def find_candidate_cent_region_on_contig(hprc2_samples, output_folder):

    ragtag_folder = "/mnt/e/workspace/cpc2_centro/res_all_annot/chrall"
    annot_folder = "/mnt/e/workspace/cpc2_centro/res_all_annot/"

    cent_on_contig_dict = {}

    for sample in hprc2_samples:

        cent_on_contig_dict[sample] = {}

        # # load ragtag res
        sample_ragtag_file = os.path.join(ragtag_folder, "{}#ragtag.scaffold.agp".format(sample))
        sample_ragtag_res = {}
        with open(sample_ragtag_file) as fin:
            for line in fin:
                if "agp" in line or "AGP" in line or "scaffold" in line:
                    continue

                line_split = line.strip().split("\t")

                chrom, ref_start, ref_end, contig_name, contig_start, contig_end, contig_orient = line_split[0], int(line_split[1]), int(line_split[2]), line_split[5].split("#")[-1], int(line_split[6]), int(line_split[7]), line_split[8]

                if chrom not in sample_ragtag_res:
                    sample_ragtag_res[chrom] = {}
                if contig_name not in sample_ragtag_res[chrom]:
                    sample_ragtag_res[chrom][contig_name] = []
                sample_ragtag_res[chrom][contig_name] = [ref_start, ref_end, contig_start, contig_end, contig_orient]

        # #
        # sample_cent_file = os.path.join(ragtag_folder, "{}/stats.cent_regions.txt".format(sample))

        with (open(os.path.join(output_folder, "{}.contig_chrom.txt".format(sample)), "w") as fout,  open(os.path.join(output_folder, "{}.hor_regions.bed".format(sample)), "w") as hor_fout):
            for chrom in chrom_order:

                chrom_with_tag = chrom + "_RagTag"
                if chrom_with_tag not in sample_ragtag_res:
                    continue

                # # load asat hor
                # core_cent_start, core_cent_end, core_cent_length  = int(line_split[1]), int(line_split[2]), int(line_split[3])

                # # load hor region
                sample_chrom_hor_file = os.path.join(annot_folder, "{}/{}#stv_hor_osat_te.pha.bed".format(chrom, sample))
                if not os.path.exists(sample_chrom_hor_file):
                    continue

                # hor_cent_start, hor_cent_end = np.inf, -np.inf
                with open(sample_chrom_hor_file) as fin:
                    fin_lines = fin.readlines()
                    if len(fin_lines) == 0:
                        continue

                    first_line_split = fin_lines[0].strip().split()
                    last_line_split = fin_lines[-1].strip().split()

                    hor_shift = int(first_line_split[0].split(":")[1].split("-")[0])
                    hor_cent_start = int(first_line_split[1]) + hor_shift
                    hor_cent_end = int(last_line_split[2]) + hor_shift

                    hor_fout.write("{}\t{}\t{}\n".format(chrom, hor_cent_start, hor_cent_end))

                hor_cent_length = hor_cent_end - hor_cent_start + 1

                for contig_name in sample_ragtag_res[chrom_with_tag]:

                    ref_start, ref_end, contig_start, contig_end, contig_orient = sample_ragtag_res[chrom_with_tag][contig_name]

                    if ref_start <= hor_cent_start < hor_cent_end <= ref_end:
                        distance_to_ref_start = hor_cent_start - ref_start
                        if contig_orient == "+":
                            cut_contig_start = contig_start + distance_to_ref_start
                            cut_contig_end = cut_contig_start + hor_cent_length

                        else:
                            cut_contig_end = contig_end - distance_to_ref_start
                            cut_contig_start = cut_contig_end - hor_cent_length

                        cent_on_contig_dict[sample][contig_name] = [chrom, cut_contig_start, cut_contig_end]

                        fout.write("{}\t{}\t{}\t{}\t{}\n".format(contig_name, chrom, ref_start, ref_end, contig_orient))


    return cent_on_contig_dict


def extract_methy_for_cent(sample, sample_contig_info, methy_folder, output_folder, bin_size=5000):

    sample_name, hap_name = sample.split("_")[0], "hap1" if "pat" in sample or "hap1" in sample else "hap2"
    print(sample, sample_name, hap_name)
    # # initial a 5k bin
    sample_contig_bin_dict = {}
    for contig_name in sample_contig_info:
        contig_cent_length = sample_contig_info[contig_name][2] - sample_contig_info[contig_name][1]

        contig_bin_num = int(contig_cent_length / bin_size) + 1
        sample_contig_bin_dict[contig_name] = [[] for i in range(contig_bin_num)]


    # # traverse methy file
    sample_methy_file = os.path.join(os.path.join(methy_folder, "{}.methylation.PacBio.{}.methylc.gz".format(sample_name, hap_name)))

    if not os.path.exists(sample_methy_file):
        return -1

    sample_raw_out_file = os.path.join(output_folder, "{}.cent_methy.bed".format(sample))
    sample_bin_out_file = os.path.join(output_folder, "{}.cent_methy.bin.bed".format(sample, bin_size))

    # with gzip.open(sample_methy_file, "rt") as fin, open(sample_raw_out_file, "w") as fout:
    with gzip.open(sample_methy_file, "rt") as fin:

        for line in fin:
            line_split = line.strip().split()

            contig_name, start, end, methy_val = line_split[0].split("#")[-1], int(line_split[1]), int(line_split[2]), float(line_split[4])

            if contig_name not in sample_contig_info:
                continue
            if sample_contig_info[contig_name][1] < start < end < sample_contig_info[contig_name][2]:
                # fout.write("{}\t{}\t{}\t{}\n".format(contig_name, start, end, methy_val))
                bin_index = int((start - sample_contig_info[contig_name][1]) / bin_size)
                sample_contig_bin_dict[contig_name][bin_index].append(methy_val)

    # print(sample_contig_bin_dict)

    # # average and output each bin
    with open(sample_bin_out_file, "w") as fout:
        for contig_name in sample_contig_bin_dict:
            contig_chrom = sample_contig_info[contig_name][0]
            contig_start = sample_contig_info[contig_name][1]

            for bin_index in range(len(sample_contig_bin_dict[contig_name])):
                bin_vals = sample_contig_bin_dict[contig_name][bin_index]
                if len(bin_vals) == 0:
                    continue

                avg_bin_vals = round(np.average(bin_vals) * 100, 2)

                bin_start = contig_start + bin_index * bin_size
                bin_end = bin_start + bin_size

                fout.write("{}\t{}\t{}\t{}\n".format(contig_name, bin_start, bin_end, avg_bin_vals))


def run_cdr_finder(bin_bed, output_folder):

    # # default params from cdf-finder
    thr_height_perc_valley = 0.34
    thr_prom_perc_valley = 0.3
    #
    # # cenmap para
    # thr_height_perc_valley = 0.2
    # thr_prom_perc_valley = 0.1

    # # test para
    # thr_height_perc_valley = 0.5
    # thr_prom_perc_valley = 0.5

    sample = bin_bed.replace(".cent_methy.bin.clean.bed", "")

    sample_bin_bed = os.path.join(output_folder, "{}.cent_methy.bin.clean.bed".format(sample))
    sample_cdr_bed = os.path.join(output_folder, "{}.cdr.bed".format(sample))

    os.system("python CDR-finder/cdr_finder.py -i {} -o {} --thr_height_perc_valley {} --thr_prom_perc_valley {} --bp_edge 500000 --extend_edges_std -1".format(sample_bin_bed, sample_cdr_bed, thr_height_perc_valley, thr_prom_perc_valley))


def merge_res(output_folder, thresh):

    with (open(os.path.join(output_folder, "all.merged_res.l_{}.stats.txt".format(thresh)), "w") as stats_fout,
          open(os.path.join(output_folder, "all.merged_res.l_{}.details.txt".format(thresh)), "w") as details_fout):

        cdr_bed_cnt = 0
        for cdr_bed in os.listdir(output_folder):
            if ".cdr.bed" not in cdr_bed:
                continue

            cdr_bed_cnt += 1
            # if cdr_bed_cnt > 10:
            #     break

            sample = cdr_bed.replace(".cdr.bed", "")
            print(sample)

            # if sample != "HG005_mat_hprc_r2_v1.0.1.fa.gz":
            #     continue

            sample_cdr_bed = os.path.join(output_folder, "{}.cdr.bed".format(sample))
            if not os.path.exists(sample_cdr_bed):
                continue

            # # load bin methylation
            sample_bin_bed = os.path.join(output_folder, "{}.cent_methy.bin.clean.bed".format(sample))
            sample_bin_dict = {}
            with open(sample_bin_bed) as fin:
                for line in fin:
                    line_split = line.strip().split()
                    bin_chrom, bin_start, bin_end, bin_val = line_split[0], int(line_split[1]), int(line_split[2]), float(line_split[3])
                    if bin_chrom not in sample_bin_dict:
                        sample_bin_dict[bin_chrom] = []
                    sample_bin_dict[bin_chrom].append([bin_start, bin_end, bin_val])

            # # load sample cent info
            sample_cent_file = os.path.join(output_folder, "{}.hor_regions.bed".format(sample))
            sample_cent_dict = {}
            with open(sample_cent_file) as fin:
                for line in fin:
                    line_split = line.strip().split()
                    cent_chrom, cent_start, cent_end =  line_split[0], int(line_split[1]), int(line_split[2])
                    if cent_chrom not in sample_cent_dict:
                        sample_cent_dict[cent_chrom] = []
                    sample_cent_dict[cent_chrom] = [cent_start, cent_end]

            # # load sample TE info
            sample_TE_dict = {}
            for chrom in ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]:
                sample_TE_dict[chrom + "_RagTag"] = []

                sample_TE_file = os.path.join(annot_folder, "{}/{}#rpmask.bed".format(chrom, sample))
                if not os.path.exists(sample_TE_file):
                    continue

                with open(sample_TE_file) as fin:
                    for line in fin:
                        line_split = line.strip().split()
                        TE_chrom, TE_shift, TE_start, TE_end, TE_type = line_split[0].split(":")[0], int(line_split[0].split(":")[1].split("-")[0]), int(line_split[1]), int(line_split[2]), line_split[3].split("/")[0]
                        # if TE_type in ["Simple_repeat", "Satellite", "Low_complexity"]:
                        #     continue
                        if TE_type not in ["LINE", "SINE", "LTR", "DNA"]:
                            continue

                        sample_TE_dict[TE_chrom].append([TE_start + TE_shift, TE_end + TE_shift, TE_type])


            # # load cdr
            df = pd.read_csv(sample_cdr_bed, sep="\t", header=None, usecols=[0, 1, 2], names=["chr", "start", "end"])
            df = df.sort_values(["chr", "start", "end"]).reset_index(drop=True)

            details_info = []
            stats_cdr_length_dict = {}
            stats_asat_length_dict = {}
            stats_dip_cnt_dict = {}
            stats_methy_val_dict = {}
            stats_normal_pos_dict = {}
            stats_TE_dict = {}

            for chrom, group in df.groupby("chr", sort=False):
                current_start = None
                current_end = None

                for start, end in zip(group["start"], group["end"]):

                    if current_start is None:
                        current_start = start
                        current_end = end

                    elif start - current_end <= thresh:
                        current_end = max(current_end, end)

                    else:
                        # # calculate average methylation values
                        methy_vals = []
                        for bin_start, bin_end, bin_val in sample_bin_dict[chrom]:
                            if current_start <= bin_start <= bin_end <= current_end:
                                methy_vals.append(bin_val)

                        # # calculate the normalized pos of cdr
                        normalized_pos = round((current_start - sample_cent_dict[chrom][0]) / (sample_cent_dict[chrom][1] - sample_cent_dict[chrom][0]), 2)

                        # # calculate the overlapped TE length
                        TE_info = 0
                        for TE_start, TE_end, TE_type in sample_TE_dict[chrom + "_RagTag"]:
                            if current_start < TE_end and TE_start < current_end:
                                # TE_info += 1
                                TE_info += (TE_end - TE_start + 1)

                        # # append to details info
                        details_info.append([chrom, current_start, current_end, round(np.average(methy_vals), 2) , normalized_pos, TE_info])

                        # # append to stats info
                        if chrom not in stats_cdr_length_dict:
                            stats_cdr_length_dict[chrom] = 0
                            stats_dip_cnt_dict[chrom] = 0
                            stats_methy_val_dict[chrom] = []
                            stats_normal_pos_dict[chrom] = []
                            stats_asat_length_dict[chrom] = -1
                            stats_TE_dict[chrom] = 0
                        stats_cdr_length_dict[chrom] += (current_end - current_start)
                        stats_dip_cnt_dict[chrom] += 1
                        stats_methy_val_dict[chrom].extend(methy_vals)
                        stats_normal_pos_dict[chrom].append(normalized_pos)
                        stats_asat_length_dict[chrom] = (sample_cent_dict[chrom][1] - sample_cent_dict[chrom][0] + 1)
                        stats_TE_dict[chrom] += TE_info

                        current_start = start
                        current_end = end

                if current_start is not None:
                    # # calculate average methylation values
                    methy_vals = []
                    for bin_start, bin_end, bin_val in sample_bin_dict[chrom]:
                        if current_start <= bin_start <= bin_end <= current_end:
                            methy_vals.append(bin_val)

                    # # calculate the normalized pos of cdr
                    normalized_pos = round((current_start - sample_cent_dict[chrom][0]) / (sample_cent_dict[chrom][1] - sample_cent_dict[chrom][0]), 2)

                    # # calculate the overlapped TE length
                    TE_info = 0
                    for TE_start, TE_end, TE_type in sample_TE_dict[chrom + "_RagTag"]:
                        if current_start < TE_end and TE_start < current_end:
                            TE_info += 1

                    # # append to details info
                    details_info.append([chrom, current_start, current_end, round(np.average(methy_vals), 2), normalized_pos, TE_info])

                    # # append to stats info
                    if chrom not in stats_cdr_length_dict:
                        stats_cdr_length_dict[chrom] = 0
                        stats_dip_cnt_dict[chrom] = 0
                        stats_methy_val_dict[chrom] = []
                        stats_normal_pos_dict[chrom] = []
                        stats_asat_length_dict[chrom] = -1
                        stats_TE_dict[chrom] = 0
                    stats_cdr_length_dict[chrom] += (current_end - current_start)
                    stats_dip_cnt_dict[chrom] += 1
                    stats_methy_val_dict[chrom].extend(methy_vals)
                    stats_normal_pos_dict[chrom].append(normalized_pos)
                    stats_asat_length_dict[chrom] = (sample_cent_dict[chrom][1] - sample_cent_dict[chrom][0] + 1)
                    stats_TE_dict[chrom] += TE_info

            # # output for each sample
            # merged_bin_file = os.path.join(output_folder, cdr_bed.replace(".bed", ".details_info.l_{}.bed".format(thresh)))
            # print(merged_bin_file)
            # with open(merged_bin_file, "w") as sample_fout:
            #     for chrom, start, end in details_info:
            #         sample_fout.write("{}\t{}\t{}\n".format(chrom, start, end))

            # # output details info
            for chrom, start, end, val, pos, TE in details_info:
                details_fout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(sample, chrom, start, end, val, pos, TE))

            # # output stats info
            for chrom in stats_cdr_length_dict:
                stats_fout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(sample, chrom.replace("_RagTag", ""), stats_cdr_length_dict[chrom], stats_dip_cnt_dict[chrom], round(np.average(stats_methy_val_dict[chrom]), 2), round(np.average(stats_normal_pos_dict[chrom]), 2), stats_asat_length_dict[chrom], stats_TE_dict[chrom]))


            # break




def compare_CDR_among_cluster(output_folder, pha_input_folder):

    cdr_res = "/mnt/e/workspace/cpc2_centro/res_all_post/res_hprc2_cdr_hor/all.merged_res.l_0.stats.txt"

    compare_mode = "length" # # choose from length, cnt, methy, pos, te

    output_folder = os.path.join(output_folder, "cdr_compare_{}".format(compare_mode))
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    # chrom_order = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"]
    # chrom_order = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]

    with open(os.path.join(output_folder, "res.txt"), "w") as fout, open(os.path.join(output_folder, "res.stats.txt"), "w") as stats_fout:
        for spec_chrom in chrom_order:
            print(spec_chrom)

            chrom_centhap_res = os.path.join(pha_input_folder, "{}.features.c2.txt".format(spec_chrom))
            cluster_CDR_dict = {}

            # # load cluste res
            sample_cluster_dict = {}
            with open(chrom_centhap_res) as fin:
                for line in fin:
                    if "Sample" in line:
                        continue

                    line_split = line.strip().split()
                    sample, cluster = line_split[0], line_split[-1]

                    if cluster == "-1":
                        continue
                    sample_cluster_dict[sample] = cluster

                    if cluster not in cluster_CDR_dict:
                        cluster_CDR_dict[cluster] = []

            # # load CDR res
            with open(cdr_res) as fin:
                for line in fin:
                    line_split = line.strip().split()

                    sample, chrom = line_split[0], line_split[1]
                    if compare_mode == "length":
                        cdr_info = int(line_split[2])
                    elif compare_mode == "cnt":
                        cdr_info = int(line_split[3])
                    elif compare_mode == "methy":
                        cdr_info = float(line_split[4])
                    elif compare_mode == "pos":
                        cdr_info = float(line_split[5])
                    elif compare_mode == "te":
                        cdr_info = float(line_split[7])

                    else:
                        print("No this mode", compare_mode)
                        exit()

                    if chrom != spec_chrom:
                        continue

                    if sample not in sample_cluster_dict:
                        continue

                    sample_cluster = sample_cluster_dict[sample]

                    cluster_CDR_dict[sample_cluster].append(cdr_info)

            all_clusters = sorted(list(cluster_CDR_dict.keys()))

            if "-1" in all_clusters:
                all_clusters.remove("-1")

            fout.write(spec_chrom + "\n")

            try:
                groups = {cluster: np.asarray(cluster_CDR_dict[cluster], dtype=float) for cluster in all_clusters}
                _, kruskal_p = kruskal(*groups.values())
                kruskal_p = round(kruskal_p, 3) if kruskal_p > 0.01 else format(kruskal_p, '.1e')
            except:
                kruskal_p = "."
            stats_fout.write("{}, Kruskal-Wallis test, P={}\n".format(spec_chrom, kruskal_p))

            for cluster in all_clusters:
                fout.write("{}\t{}\t{}".format(spec_chrom, cluster, len(cluster_CDR_dict[cluster])))
                stats_fout.write("{}\t{}\t{}".format(spec_chrom, cluster, len(cluster_CDR_dict[cluster])))

                if len(cluster_CDR_dict[cluster]) != 0:
                    q1 = np.percentile(cluster_CDR_dict[cluster], 25)
                    median = np.percentile(cluster_CDR_dict[cluster], 50)
                    q3 = np.percentile(cluster_CDR_dict[cluster], 75)
                    iqr = round(q3 - q1, 3)
                    lower_bound = round(q1 - 1.5 * iqr, 3)
                    upper_bound = round(q3 + 1.5 * iqr, 3)
                else:
                    q1 = "."
                    median = "."
                    q3 = "."
                    iqr = "."
                    lower_bound = "."
                    upper_bound = "."

                p_value_list = []

                for target_cluster in all_clusters:
                    if len(cluster_CDR_dict[target_cluster]) == 0:
                        p_value_list.append(".")
                    else:
                        p_value = mannwhitneyu(cluster_CDR_dict[cluster], cluster_CDR_dict[target_cluster], alternative="two-sided").pvalue
                        p_value_list.append(str(round(p_value, 3) if p_value > 0.01 else format(p_value, '.3e')))

                    if p_value < 0.001:
                        pass_flag = "***"
                    elif p_value < 0.01:
                        pass_flag = "**"
                    elif p_value < 0.05:
                        pass_flag = "*"
                    else:
                        pass_flag = "ns"

                    # ks_p = round(ks_p, 3) if ks_p > 0.01 else format(ks_p, '.3e')
                    fout.write("\t{}".format(pass_flag))

                for cdr_info in cluster_CDR_dict[cluster]:
                    fout.write("\t{}".format(cdr_info))
                fout.write("\n")

                stats_fout.write("\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(median, q1, q3, iqr, lower_bound, upper_bound, "\t".join(p_value_list)))



            df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in cluster_CDR_dict.items()]))

            df_long = df.melt(var_name="cluster", value_name="value")

            if len(cluster_CDR_dict.keys()) == 2:
                fig, ax = plt.subplots(figsize=(5, 7))
            elif len(cluster_CDR_dict.keys()) == 3:
                fig, ax = plt.subplots(figsize=(6, 7))
            elif len(cluster_CDR_dict.keys()) == 4:
                fig, ax = plt.subplots(figsize=(7, 7))
            else:
                fig, ax = plt.subplots(figsize=(8, 7))

            # lower = df_long.groupby('cluster', observed=True)['value'].transform('quantile', 0.01)
            # upper = df_long.groupby('cluster', observed=True)['value'].transform('quantile', 0.99)
            # df_plot = df_long[(df_long['value'] >= lower) &(df_long['value'] <= upper)].copy()

            sns.violinplot(data=df_long,  saturation=1, x="cluster",y="value", inner="box", palette=cluster_colors, cut=0, )
            # sns.swarmplot(data=df_long, x="cluster", y="value", size=1, color="#fe8d24")

            ax.set_xlabel("Cluster", fontsize=21)
            if compare_mode == "length":
                ax.set_ylabel("CDR length (bp)", fontsize=21)
            elif compare_mode == "pos":
                ax.set_ylabel("Normalized CDR position", fontsize=21)

            # 调整横轴 tick label 字体大小
            ax.tick_params(axis="x", labelsize=21)

            # 如果需要，也可以调整纵轴 tick label 字体大小
            ax.tick_params(axis="y", labelsize=21)

            plt.tight_layout()
            plt.savefig(os.path.join(output_folder, "{}.cdr_compare.png".format(spec_chrom)))
            plt.savefig(os.path.join(output_folder, "{}.cdr_compare.svg".format(spec_chrom)))
            plt.close()


def overlap_cdr_with_hor(cdr_file):

    annot_folder = "/mnt/e/workspace/cpc2_centro/res_all_annot/"

    # # load cdr results
    with open(cdr_file) as fin, open(cdr_file.replace(".txt", ".add_hor.txt"), "w") as fout:

        # # traverse each cdr region
        cnt = 0
        for line in fin:
            cnt += 1
            # if cnt > 5000:
            #     break

            line_split = line.strip().split()

            sample, chrom, start, end = line_split[0], line_split[1], int(line_split[2]), int(line_split[3])
            # if sample != "HG00097_hap1_hprc_r2_v1.0.1.fa.gz" or chrom != "chr3_RagTag":
            #     continue
            print(sample, chrom)

            # # traverse hor
            hor_annot_file = os.path.join(annot_folder, "{}/{}#stv_hor.bed".format(chrom.replace("_RagTag", ""), sample))

            chrom_selected_hors = c1_principle_hor_dict[chrom.replace("_RagTag", "")]
            hor_cnt_dict = {hor: 0 for hor in chrom_selected_hors + ["Other"]}

            total_overlap_length = 0
            with open(hor_annot_file) as hor_fin:
                for hor_line in hor_fin:
                    line_split = hor_line.strip().split()

                    shift, target_start, target_end, hor_type = int(line_split[0].split(":")[1].split("-")[0]), int(line_split[1]), int(line_split[2]), line_split[3]

                    hor_type = re.sub(r'\{\d+\}', '{n}', hor_type)  # # for chr1, chr5 chr19
                    target_start += shift
                    target_end += shift

                    # # overlap
                    if start < target_end and target_start < end:
                        overlap_start, overlap_end = max(start, target_start), min(end, target_end)
                        overlap_length = overlap_end - overlap_start + 1
                        total_overlap_length += overlap_length

                        if hor_type in chrom_selected_hors:
                            hor_cnt_dict[hor_type] += 1
                        elif reverse_hor(hor_type) in chrom_selected_hors:
                            hor_cnt_dict[reverse_hor(hor_type)] += 1
                        else:
                            hor_cnt_dict["Other"] += 1

            hor_freq_dict = {hor: 0 for hor in chrom_selected_hors + ["Other"]}

            hor_total_cnt = sum(hor_cnt_dict.values())

            # if hor_total_cnt != 0 and total_overlap_length / (end - start + 1) > 0.5:
            if hor_total_cnt != 0:

                for hor in hor_cnt_dict:
                    hor_freq_dict[hor] = round(hor_cnt_dict[hor] / hor_total_cnt, 2)


            fout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(sample, chrom, start, end, ",".join([str(hor_freq_dict[hor]) for hor in hor_freq_dict]), ",".join([str(hor_cnt_dict[hor]) for hor in hor_cnt_dict]), total_overlap_length, (end - start + 1), round(total_overlap_length / (end - start + 1), 2)))


def clean_bin_bed_with_hor(output_folder):

    for bin_bed in os.listdir(output_folder):

        if ".cent_methy.bin.bed" not in bin_bed:
            continue
        print(bin_bed)
        sample = bin_bed.replace(".cent_methy.bin.bed", "")

        # # load hor region
        hor_region_file = os.path.join(output_folder, "{}.hor_regions.bed".format(sample))
        hor_region_dict = {}
        with open(hor_region_file) as fin:
            for line in fin:
                line_split = line.strip().split()
                chrom, start, end = line_split[0], int(line_split[1]), int(line_split[2])
                if chrom not in hor_region_dict:
                    hor_region_dict[chrom] = []
                hor_region_dict[chrom].append([start, end])

        # # load contig chrom convert file
        contig_chrom_file = os.path.join(output_folder, "{}.contig_chrom.txt".format(sample))
        contig_chrom_dict = {}
        with open(contig_chrom_file) as fin:
            for line in fin:
                line_split = line.strip().split()
                contig, chrom, contig_start_on_chrom = line_split[0], line_split[1], int(line_split[2]) - 1
                contig_chrom_dict[contig] = [chrom, contig_start_on_chrom]

        # # traverse each bin
        with open(os.path.join(output_folder, bin_bed)) as fin,   open(os.path.join(output_folder, bin_bed.replace(".bin.bed", ".bin.clean.bed")), "w") as fout:

            for line in fin:
                line_split = line.strip().split()

                contig, start, end, val = line_split[0], int(line_split[1]), int(line_split[2]), line_split[3]

                chrom, shift = contig_chrom_dict[contig][0], contig_chrom_dict[contig][1]
                start += shift
                end += shift

                for hor_start, hor_end in hor_region_dict[chrom]:

                    # # if overlap
                    if start < hor_end and hor_start < end:
                        fout.write("{}\t{}\t{}\t{}\n".format(chrom, start, end, val))
                        break


if __name__ == '__main__':

    input_folder = "/mnt/f/data_hprc2_methy_pacbio"
    annot_folder = "/mnt/e/workspace/cpc2_centro/res_all_annot"
    pha_input_folder = "/mnt/e/workspace/cpc2_centro/res_all_post/principal_hor_analysis_outlier"

    output_folder = "/mnt/e/workspace/cpc2_centro/res_all_post/res_hprc2_cdr_hor"
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    chrom_order = ["chr1", "chr2", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr16", "chr17", "chr19", "chr20", "chr21"]

    c1_principle_hor_dict = {"chr1": ["S1C1/5/19H1L.1-5(_6/4_5){n}-6", "S1C1/5/19H1L.1-3(_6/4_5){n}-6", "S1C1/5/19H1L.5_6/4", "S1C1/5/19H1L.5-6", "S1C1/5/19H1L.1-5_6/4"], "chr2": ["S2C2H1L.1-4", "S2C2H1L.1-2", "S2C2H1L.1_4", "S2C2H1L.1/3_4"],
                               "chr3": ["S01/1C3H1L.1-17", "S01/1C3H1L.1-4_12/5_6-17", "S01/1C3H1L.1-6_5/7_8-17", "S01/1C3H1L.1-4_5/7_8-17"], "chr4": ["S2C4H1L.1-19", "S2C4H1L.1-2_7-19", "S2C4H1L.1-9_3-4_12-19", "S2C4H1L.1-18", "S2C4H1L.1-3_16-19"],
                               "chr5": ["S1C1/5/19H1L.1-5(_6/4_5){n}_2/6", "S1C1/5/19H1L.1_2/4_5_2/6", "S1C1/5/19H1L.1-2_5(_6/4_5){n}_2/6", "S1C1/5/19H1L.1-5_6/4"], "chr6": ["S01C6H1L.1-18", "S01C6H1L.1-12_16-18"],
                               "chr7": ["S1C7H1L.1-6", "S1C7H1L.1-5_2-4_3-6"], "chr8": ["S2C8H1L.1-11", "S2C8H1L.1-7", "S2C8H1L.1-3_4/7_8-11"],
                               "chr9": ["S2C9H1L.1-7", "S2C9H1L.1-2_6-7", "S2C9H1L.7", "S2C9H1L.1-3_7"], "chr10": ["S1C10H1L.1-8", "S1C10H1L.1-2_5-8", "S1C10H1L.1-2_5_6/8"],
                               "chr11": ["S3C11H1L.1-5", "S3C11H1L.1"], "chr12": ["S1C12H1L.1-8", "S1C12H1L.1_6-8"],
                               "chr13": ["S2C13/21H1L.1-11", "S2C13/21H1L.8-11", "S2C13/21H1L.1-7", "S2C13/21H1L.1-2_9-11"], "chr14": ["S2C14/22H1L.1-8"],
                               "chr15": ["S2C15H1L.1-11", "S2C15H1L.1-4", "S2C15H1L.1", "S2C15H1L.7"], "chr16": ["S1C16H1L.1-10", "S1C16H1L.2-9"],
                               "chr17": ["S3C17H1L.1-16", "S3C17H1L.1-10_15_15-16", "S3C17H1L.1-13_15_15-16", "S3C17H1L.1-13_15-16"], "chr18": ["S2C18H1L.1-12", "S2C18H1L.1-2_3/7_8-12", "S2C18H1L.1-2_7_4-12"],
                               "chr19": ["S1C1/5/19H1L.5(_6/4_5){n}-6", "S1C1/5/19H1L.5-6", "S1C1/5/19H1L.1-3(_6/4_5){n}_2/6", "S1C1/5/19H1L.1", "S1C1/5/19H1L.1-5(_6/4_5){n}_2/6", "S1C1/5/19H1L.1_2/6", "S1C1/5/19H1L.1(_6/4_5){n}_2/6"], "chr20": ["S2C20H1L.1-16", "S2C20H1L.7-16"],
                               "chr21": ["S2C13/21H1L.1-11", "S2C13/21H1L.1-2_9-11", "S2C13/21H1L.1-5", "S2C13/21H1L.8-11"], "chr22": ["S2C14/22H1L.1-8", "S2C14/22H1L.7"],
                               "chrX": ["S3CXH1L.1-12"]}

    cluster_colors = {
        "0": "#547A95",  # 灰蓝
        '1': "#E5A06A",  # 柔和橙
        "2": "#8D80B8",  # 灰紫
        "3": "#6FA69A",  # 青绿
        "4": "#D4776A",  # 珊瑚红
        "5": "#C5A65A",  # 金色
    }

    # # STEP 1: load sample info
    with open("/mnt/e/workspace/cpc2_centro/res_all_post/stats.whole_core_cent.json") as file:
        whole_cent_dict = dict(json.load(file))

    hprc2_samples = []
    with open("../data/data_population_hgsvc3_hprc2_cpc2.txt") as fin:
        for line in fin:
            if "HPRC2" in line:
                hprc2_samples.append(line.strip().split()[1])

    cent_on_contig_dict = find_candidate_cent_region_on_contig(hprc2_samples, output_folder)

    # #
    # # STEP 2: extract cent methy and group by bin
    process_pool = multiprocessing.Pool(20)
    for sample in cent_on_contig_dict:
        sample_contig_info = cent_on_contig_dict[sample]
        extract_methy_for_cent(sample, sample_contig_info, input_folder, output_folder, 5000)
        # process_pool.apply_async(extract_methy_for_cent, (sample, sample_contig_info, input_folder, output_folder, 5000))
    process_pool.close()
    process_pool.join()

    # # STEP: remove methylation bins that don't overlap with hor
    clean_bin_bed_with_hor(output_folder)

    # # STEP 3: run cdr-finder
    process_pool = multiprocessing.Pool(10)
    for bin_bed in os.listdir(output_folder):
        if ".cent_methy.bin.clean.bed" not in bin_bed:
            continue
        process_pool.apply_async(run_cdr_finder, (bin_bed, output_folder))
    process_pool.close()
    process_pool.join()

    # STEP 4: merge res
    # merge_res(output_folder, thresh=150000)
    merge_res(output_folder, thresh=0)

    # # STEP 5: overlap cdr res with hor
    overlap_cdr_with_hor(os.path.join(output_folder, "all.merged_res.l_0.details.txt"))

    # # STEP 6: compare
    compare_CDR_among_cluster(output_folder, pha_input_folder)


