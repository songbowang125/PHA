import os
import json
from scipy.stats import pearsonr
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import ks_2samp, mannwhitneyu
import diptest
from scipy.signal import find_peaks
from scipy.stats import gaussian_kde
from scipy.stats import kruskal
from statsmodels.stats.multitest import multipletests
import scikit_posthocs as sp
from itertools import combinations
plt.rcParams["svg.fonttype"] = "none"


def significance_label(p):
    if p < 0.0001:
        return "****"
    elif p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return "ns"


def stats_aSat_length_from_monomer(merge_thresh="0k", spec_aSat=None):

    if merge_thresh not in ["0k", "10k"]:
        print("No this merge thresh", merge_thresh)
        exit()

    if spec_aSat not in [None, "active", "inactive"]:
        print("no this spec aSat", spec_aSat)
        exit()

    if spec_aSat is None:
        output_sat_length_json = open(os.path.join(output_folder, "stats.aSat_length_monomer.merge_{}.json".format(merge_thresh)), "w")
    else:
        output_sat_length_json = open(os.path.join(output_folder, "stats.aSat_length_monomer.{}.merge_{}.json".format(spec_aSat, merge_thresh)), "w")

    sample_sat_length_dict = {}      # # format Sat_length_collect_sample[sample][chrom] = length

    for spec_chrom in os.listdir(annot_res_folder):
        chrom_input_folder = os.path.join(annot_res_folder, spec_chrom)
        for sample_monomer_file in os.listdir(chrom_input_folder):
            if "#monomer.bed" not in sample_monomer_file:
                continue
            sample = sample_monomer_file.replace("#monomer.bed", "")

            if sample not in all_included_samples:
                continue

            if "{}_RagTag".format(spec_chrom) not in whole_core_cent[sample]:
                continue

            print(spec_chrom, sample)

            if sample not in sample_sat_length_dict:
                sample_sat_length_dict[sample] = {}

            if spec_chrom not in sample_sat_length_dict[sample]:
                sample_sat_length_dict[sample][spec_chrom] = 0

            if spec_aSat is None:
                sample_monomer_file = os.path.join(chrom_input_folder, sample_monomer_file)
                sample_monomer_file_merged = sample_monomer_file + "tmp_merge"

                if merge_thresh == "0k":
                    os.system("bedtools merge -d {} -i {} > {}".format(0, sample_monomer_file, sample_monomer_file_merged))
                else:
                    os.system("bedtools merge -d {} -i {} > {}".format(10000, sample_monomer_file, sample_monomer_file_merged))

            else:
                sample_monomer_file = os.path.join(chrom_input_folder, sample_monomer_file)
                sample_monomer_live_file = sample_monomer_file + "tmp_live"
                sample_monomer_file_merged = sample_monomer_file + "tmp_merge"

                with open(sample_monomer_file) as fin, open(sample_monomer_live_file, "w") as fout:
                    for line in fin:
                        line_split = line.strip().split()
                        start, end, monomer = int(line_split[1]), int(line_split[2]), line_split[3]

                        if spec_aSat == "active":
                            if "L." in monomer:
                                fout.write(line)
                        else:
                            if "L." not in monomer:
                                fout.write(line)

                if merge_thresh == "0k":
                    os.system("bedtools merge -d {} -i {} > {}".format(0, sample_monomer_live_file, sample_monomer_file_merged))
                else:
                    os.system("bedtools merge -d {} -i {} > {}".format(10000, sample_monomer_live_file, sample_monomer_file_merged))

                os.remove(sample_monomer_live_file)

            # # traverse file
            with open(sample_monomer_file_merged) as fin:
                for line in fin:
                    line_split = line.strip().split()
                    start, end = int(line_split[1]), int(line_split[2])

                    sample_sat_length_dict[sample][spec_chrom] += (end - start) + 1

            os.remove(sample_monomer_file_merged)

    json.dump(sample_sat_length_dict, output_sat_length_json, indent=4)
    output_sat_length_json.close()



def stats_oSAT_lengths_from_rpmask(merge_thresh):

    if merge_thresh not in ["0k", "10k"]:
        print("No this merge thresh", merge_thresh)
        exit()


    allowed_SAT_types = ["HSat1A", "HSat1B", "HSat2", "HSat3", "BSat", "GSat"]

    # # begin
    for Sat in allowed_SAT_types:
        output_sat_length_json = open(os.path.join(output_folder, "stats.{}_length_rpmask.merge_{}.json".format(Sat, merge_thresh)), "w")

        sample_sat_length_dict = {}  # # format Sat_length_collect_sample[sample][chrom] = length

        for spec_chrom in os.listdir(annot_res_folder):
            chrom_input_folder = os.path.join(annot_res_folder, spec_chrom)

            for sample_oSat_file in os.listdir(chrom_input_folder):
                if "#osat.bed" not in sample_oSat_file:
                    continue
                sample = sample_oSat_file.replace("#osat.bed", "")

                if sample not in all_included_samples:
                    continue

                if "{}_RagTag".format(spec_chrom) not in whole_peri_cent[sample]:
                    continue

                print(Sat, sample, spec_chrom)
                if sample not in sample_sat_length_dict:
                    sample_sat_length_dict[sample] = {}
                if spec_chrom not in sample_sat_length_dict[sample]:
                    sample_sat_length_dict[sample][spec_chrom] = 0

                sample_oSat_file = os.path.join(chrom_input_folder, sample_oSat_file)
                sample_oSat_file_merged = sample_oSat_file + "tmp_merge"

                if merge_thresh == "0k":
                    os.system("grep {} {} | bedtools merge -d {} > {}".format(Sat, sample_oSat_file, 0, sample_oSat_file_merged))
                else:
                    os.system("grep {} {} | bedtools merge -d {} > {}".format(Sat, sample_oSat_file, 10000, sample_oSat_file_merged))

                # # traverse file
                with open(sample_oSat_file_merged) as fin:
                    for line in fin:
                        line_split = line.strip().split()
                        start, end = int(line_split[1]), int(line_split[2])

                        sample_sat_length_dict[sample][spec_chrom] += (end - start) + 1

                os.remove(sample_oSat_file_merged)

        json.dump(sample_sat_length_dict, output_sat_length_json, indent=4)
        output_sat_length_json.close()



def check_length_corr_among_popu(output_folder, json_1, json_2, popu_mode):

    if popu_mode == "hyper":
        popu_dict = hyper_popu_dict
        popu_order = hyper_popu_order
    else:
        popu_dict = detail_popu_dict
        popu_order = detail_popu_order


    chrom_order = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]

    with open(os.path.join(output_folder, "stats.whole_core_cent.json"), 'r') as file:
        whole_core_cent = dict(json.load(file))

    with open(os.path.join(output_folder, json_1)) as file:
        length_dict_1 = dict(json.load(file))

    with open(os.path.join(output_folder, json_2)) as file:
        length_dict_2 = dict(json.load(file))

    for chrom in chrom_order:
        for spec_popu in popu_order:
            length_list_1 = []
            length_list_2 = []


            for sample in length_dict_1:
                if chrom + "_RagTag" not in whole_core_cent[sample]:
                    continue

                if chrom not in length_dict_1[sample] or chrom not in length_dict_2[sample]:
                    continue

                sample_popu = popu_dict[sample]
                if sample_popu != spec_popu:
                    continue

                length_list_1.append(length_dict_1[sample][chrom])
                length_list_2.append(length_dict_2[sample][chrom])

            r, p = pearsonr(length_list_1, length_list_2)

            print(chrom, spec_popu, round(r, 3), p, sep="\t")

def check_length_corr_among_cluster(output_folder, json_1, json_2, cluster_folder):



    chrom_order = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]

    with open(os.path.join(output_folder, "stats.whole_core_cent.json"), 'r') as file:
        whole_core_cent = dict(json.load(file))

    with open(os.path.join(output_folder, json_1)) as file:
        length_dict_1 = dict(json.load(file))

    with open(os.path.join(output_folder, json_2)) as file:
        length_dict_2 = dict(json.load(file))

    for chrom in chrom_order:

        # # load cluste res
        chrom_cluster_file = os.path.join(cluster_folder, "{}.features.c2.txt".format(chrom))

        collect_dict = {}
        sample_cluster_dict = {}
        with open(chrom_cluster_file) as fin:
            for line in fin:
                if "Sample" in line:
                    continue

                line_split = line.strip().split()
                sample, cluster = line_split[0], line_split[-1]

                sample_cluster_dict[sample] = cluster

                if cluster not in collect_dict:
                    collect_dict[cluster] = []
        cluster_order = sorted(list(collect_dict.keys()))

        for spec_cluster in cluster_order:
            length_list_1 = []
            length_list_2 = []

            for sample in length_dict_1:
                if chrom + "_RagTag" not in whole_core_cent[sample]:
                    continue

                if chrom not in length_dict_1[sample] or chrom not in length_dict_2[sample]:
                    continue

                if sample not in sample_cluster_dict:
                    continue

                sample_cluster = sample_cluster_dict[sample]
                if sample_cluster != spec_cluster:
                    continue

                length_list_1.append(length_dict_1[sample][chrom])
                length_list_2.append(length_dict_2[sample][chrom])

            r, p = pearsonr(length_list_1, length_list_2)

            print(chrom, spec_cluster, round(r, 3), p, sep="\t")

def genearte_SAT_length_matrix(input_folder, popu_mode, spec_aSat=None):
    if spec_aSat not in [None, "active", "inactive"]:
        print("no this spec aSat", spec_aSat)
        exit()

    if spec_aSat is None:
        plot_output_folder = os.path.join(input_folder, "Sat_length_{}_popu".format(popu_mode))
    else:
        plot_output_folder = os.path.join(input_folder, "Sat_length_{}_popu_{}".format(popu_mode, spec_aSat))
    if not os.path.exists(plot_output_folder):
        os.mkdir(plot_output_folder)
    chrom_order = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]

    if spec_aSat is None:
        sat_order = ["aSat", "HSat1A", "HSat1B", "HSat2", "HSat3", "BSat", "GSAT"]
    else:
        sat_order = ["aSat"]

    for sat in sat_order:
        # # load data
        if sat == "aSat":
            if spec_aSat is None:
                with open(os.path.join(input_folder, "stats.aSat_length_monomer.merge_0k.json")) as file:
                    sample_sat_length_dict = dict(json.load(file))
            else:
                with open(os.path.join(input_folder, "stats.aSat_length_monomer.{}.merge_0k.json".format(spec_aSat))) as file:
                    sample_sat_length_dict = dict(json.load(file))

        else:
            with open(os.path.join(input_folder, "stats.{}_length_rpmask.merge_0k.json".format(sat))) as file:
                sample_sat_length_dict = dict(json.load(file))

        with open(os.path.join(plot_output_folder, "res.{}_length.matrix.txt".format(sat)), "w") as fout:
            fout.write("Sample\t{}\n".format("\t".join(chrom_order)))

            for sample in sample_sat_length_dict:

                length_list = []
                for chrom in chrom_order:
                    if chrom not in sample_sat_length_dict[sample]:
                        length_list.append(".")
                    else:
                        length_list.append(str(sample_sat_length_dict[sample][chrom]))

                fout.write("{}\t{}\n".format(sample, "\t".join(length_list)))


def plot_SAT_length_among_popu(input_folder, popu_mode, spec_aSat=None):

    if spec_aSat not in [None, "active", "inactive"]:
        print("no this spec aSat", spec_aSat)
        exit()

    if spec_aSat is None:
        plot_output_folder = os.path.join(input_folder, "Sat_length_{}_popu".format(popu_mode))
    else:
        plot_output_folder = os.path.join(input_folder, "Sat_length_{}_popu_{}".format(popu_mode, spec_aSat))

    if not os.path.exists(plot_output_folder):
        os.mkdir(plot_output_folder)

    # # set info
    chrom_order = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]

    if spec_aSat is None:
        sat_order = ["aSat", "HSat1A", "HSat1B", "HSat2", "HSat3", "BSat", "GSAT"]
    else:
        sat_order = ["aSat"]


    if popu_mode == "hyper":
        popu_dict = hyper_popu_dict
        popu_order = hyper_popu_order

    else:
        popu_dict = detail_popu_dict
        popu_order = detail_popu_order

    for sat in sat_order:
        # # load data
        if sat == "aSat":
            if spec_aSat is None:
                with open(os.path.join(input_folder, "stats.aSat_length_monomer.merge_0k.json")) as file:
                    sample_sat_length_dict = dict(json.load(file))
            else:
                with open(os.path.join(input_folder, "stats.aSat_length_monomer.{}.merge_0k.json".format(spec_aSat))) as file:
                    sample_sat_length_dict = dict(json.load(file))

        else:
            with open(os.path.join(input_folder, "stats.{}_length_rpmask.merge_0k.json".format(sat))) as file:
                sample_sat_length_dict = dict(json.load(file))

        with open(os.path.join(plot_output_folder, "res.{}_length.txt".format(sat)), "w") as fout, open(os.path.join(plot_output_folder, "res.{}_length.data.txt".format(sat)), "w") as data_fout:
            all_p = []
            for chrom in chrom_order:
                fout.write("\n=====================================================")

                collect_dict = {popu: [] for popu in popu_order}

                for sample in sample_sat_length_dict:

                    if sample == "T2T-CHM13.fa.gz":
                        continue
                    if sample not in popu_dict or popu_dict[sample] == "NOUSE":
                        continue
                    if chrom not in sample_sat_length_dict[sample]:
                        continue

                    sample_sat_length = int(sample_sat_length_dict[sample][chrom])

                    collect_dict[popu_dict[sample]].append(sample_sat_length)

                # # output the raw data
                for popu in popu_order:
                    data_fout.write("{}\t{}\t{}\n".format(chrom, popu, "\t".join([str(val) for val in collect_dict[popu]])))

                # # the overall tests
                groups = { popu: np.asarray(lengths, dtype=float) for popu, lengths in collect_dict.items() }

                _, kruskal_p = kruskal(*groups.values())

                fout.write("\nkruskal P:\n")
                fout.write("{}, Kruskal-Wallis test, P={}\n".format(chrom,  round(kruskal_p, 3) if kruskal_p > 0.01 else format(kruskal_p, '.1e')))

                all_p.append(kruskal_p)

                if kruskal_p < 0.05:
                    # # using dunn test
                    # df = pd.DataFrame([ {"population": popu, "length": value} for popu, values in groups.items() for value in values ])
                    #
                    # dunn = sp.posthoc_dunn( df, val_col="length", group_col="population", p_adjust="fdr_bh" )
                    #
                    # fout.write("\nDunn pairwise comparisons (BH-adjusted P):\n")
                    # fout.write(dunn.round(5).to_csv(sep="\t"))


                    # # using mwu test
                    fout.write("\nMannwhitneyu pairwise comparisons:\n")

                    populations = popu_order
                    pairs = list(combinations(populations, 2))

                    for pop1, pop2 in pairs:
                        x = np.asarray(collect_dict[pop1], dtype=float)
                        y = np.asarray(collect_dict[pop2], dtype=float)

                        stat, muw_p = mannwhitneyu( x, y, alternative="two-sided", method="auto")
                        muw_p_label = significance_label(muw_p)
                        muw_p = round(muw_p, 3) if muw_p > 0.01 else format(muw_p, '.1e')

                        stat, ks_p = ks_2samp(x, y)
                        ks_p_label = significance_label(ks_p)
                        ks_p = round(ks_p, 3) if ks_p > 0.01 else format(ks_p, '.1e')

                        fout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chrom, pop1, pop2, muw_p, muw_p_label, ks_p, ks_p_label))

                fout.write("\nStats values:\n")

                for popu in popu_order:
                    if len(collect_dict[popu]) == 0:
                        continue

                    q1 = int(np.percentile(collect_dict[popu], 25))
                    median = int(np.percentile(collect_dict[popu], 50))
                    q3 = int(np.percentile(collect_dict[popu], 75))
                    iqr = q3 - q1
                    lower = int(q1 - 1.5 * iqr)
                    upper = int(q3 + 1.5 * iqr)
                    mean = int(np.average(collect_dict[popu]))

                    dip, dip_p = diptest.diptest(np.asarray(collect_dict[popu]))
                    dip_p_label = significance_label(dip_p)
                    # dip = round(dip, 3)
                    dip_p = round(dip_p, 3) if dip_p > 0.01 else format(dip_p, '.1e')


                    p_values = []
                    for target_popu in popu_order:
                        stat, muw_p = mannwhitneyu( collect_dict[popu], collect_dict[target_popu], alternative="two-sided", method="auto")
                        muw_p_label = significance_label(muw_p)
                        muw_p = round(muw_p, 3) if muw_p > 0.01 else format(muw_p, '.1e')

                        p_values.append(str(muw_p))


                    fout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chrom, popu, len(collect_dict[popu]), mean, median, lower, upper, dip_p, dip_p_label, "\t".join(p_values)))


                    #
                    # if popu_mode == "hyper":
                    #     length_diff = int(np.percentile(collect_dict[popu], 50)) - int(np.percentile(collect_dict["Other"], 50))
                    #
                    #     length_p = mannwhitneyu(collect_dict[popu], collect_dict["Other"], alternative="two-sided").pvalue
                    #     length_p_label = significance_label(length_p)
                    #     length_p = round(length_p, 3) if length_p > 0.01 else format(length_p, '.1e')
                    #
                    #     ks_length_p =  ks_2samp(collect_dict[popu], collect_dict["Other"]).pvalue
                    #     ks_length_p_label = significance_label(ks_length_p)
                    #     ks_length_p = round(ks_length_p, 3) if ks_length_p > 0.01 else format(ks_length_p, '.1e')
                    #
                    #
                    #     fout.write("{}\t{}\t{}\t{}\tMid:{}\tQ1:{}\tQ3:{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chrom, popu, len(collect_dict[popu]), mean, median, lower, upper, dip, dip_p, dip_p_label, length_diff, length_p, length_p_label, ks_length_p, ks_length_p_label))
                    #
                    # else:
                    #
                    #
                    #     fout.write("{}\t{}\t{}\t{}\tMid:{}\tQ1:{}\tQ3:{}\n".format(chrom, popu, len(collect_dict[popu]), mean, median, lower, upper))

                # # plot R
                all_values, all_groups = [], []

                for popu in popu_order:
                    all_values.extend(collect_dict[popu])
                    all_groups.extend([popu] * len(collect_dict[popu]))

                df = pd.DataFrame({"value": all_values, "group": all_groups})

                # # plot a boxplot
                fig, ax = plt.subplots(figsize=(3, 5))
                sns.boxplot(data=df, x="group", y="value", color="#5975a5", showfliers=False, ax=ax)
                sns.despine(ax=ax)

                ax.set_xlabel("")
                ax.set_ylabel("")
                ax.yaxis.get_offset_text().set_visible(False)
                plt.tight_layout()
                plt.savefig(os.path.join(plot_output_folder, "res.{}_length.{}.boxplot.png".format(sat, chrom)))
                plt.close()

                # # plot a violin plot
                # df["comparison"] = "All"

                fig, ax = plt.subplots(figsize=(3, 5))
                sns.violinplot(data=df, x="group", y="value", cut=0, linewidth=1, color="#5975a5", ax=ax)
                # sns.violinplot(data=df, x="comparison", y="value", hue="group", cut=0, linewidth=1, ax=ax, split=True)

                sns.despine(ax=ax)

                ax.set_xlabel("")
                ax.set_ylabel("")
                ax.yaxis.get_offset_text().set_visible(False)
                plt.tight_layout()
                plt.savefig(os.path.join(plot_output_folder, "res.{}_length.{}.violin.png".format(sat, chrom)))
                plt.close()

            # fdr_values = multipletests(all_p, method="fdr_bh")[1]
            # print(fdr_values)

def plot_SAT_length_among_cluster(input_folder, cluster_folder, spec_aSat=None):

    cluster_colors = {
        "0": "#547A95",  # 灰蓝
        '1': "#E5A06A",  # 柔和橙
        "2": "#8D80B8",  # 灰紫
        "3": "#6FA69A",  # 青绿
        "4": "#D4776A",  # 珊瑚红
        "5": "#C5A65A",  # 金色
    }

    if spec_aSat not in [None, "active", "inactive"]:
        print("no this spec aSat", spec_aSat)
        exit()

    if spec_aSat is None:
        plot_output_folder = os.path.join(input_folder, "Sat_length_cluster_2")
    else:
        plot_output_folder = os.path.join(input_folder, "Sat_length_cluster_2_{}".format(spec_aSat))

    if not os.path.exists(plot_output_folder):
        os.mkdir(plot_output_folder)

    # # set info
    chrom_order = ["chr1", "chr2", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr16", "chr17", "chr19", "chr20", "chr21"]

    if spec_aSat is None:
        sat_order = ["aSat", "HSat1A", "HSat1B", "HSat2", "HSat3", "BSat", "GSAT"]
    else:
        sat_order = ["aSat"]

    for sat in sat_order:
        # # load data
        if sat == "aSat":
            if spec_aSat is None:
                with open(os.path.join(input_folder, "stats.aSat_length_monomer.merge_0k.json")) as file:
                    sample_sat_length_dict = dict(json.load(file))
            else:
                with open(os.path.join(input_folder, "stats.aSat_length_monomer.{}.merge_0k.json".format(spec_aSat))) as file:
                    sample_sat_length_dict = dict(json.load(file))

        else:
            with open(os.path.join(input_folder, "stats.{}_length_rpmask.merge_0k.json".format(sat))) as file:
                sample_sat_length_dict = dict(json.load(file))

        with open(os.path.join(plot_output_folder, "res.{}_length.txt".format(sat)), "w") as fout, open(os.path.join(plot_output_folder, "res.{}_length.data.txt".format(sat)), "w") as data_fout:
            for chrom in chrom_order:

                # # load cluste res
                chrom_cluster_file = os.path.join(cluster_folder, "{}.features.c2.txt".format(chrom))

                collect_dict = {}
                cluster_dict = {}
                with open(chrom_cluster_file) as fin:
                    for line in fin:
                        if "Sample" in line:
                            continue

                        line_split = line.strip().split()
                        sample, cluster = line_split[0], line_split[-1]
                        if cluster == "-1":
                            continue

                        cluster_dict[sample] = cluster

                        if cluster not in collect_dict:
                            collect_dict[cluster] = []
                cluster_order = sorted(list(collect_dict.keys()))
                # ############
                collect_dict = {cluster: [] for cluster in cluster_order}

                for sample in sample_sat_length_dict:

                    if sample == "T2T-CHM13.fa.gz":
                        continue
                    if sample not in cluster_dict or cluster_dict[sample] == "-1":
                        continue
                    if chrom not in sample_sat_length_dict[sample]:
                        continue

                    sample_sat_length = int(sample_sat_length_dict[sample][chrom])

                    collect_dict[cluster_dict[sample]].append(sample_sat_length)

                # # output the raw data
                for cluster in cluster_order:
                    data_fout.write("{}\t{}\t{}\n".format(chrom, cluster, "\t".join([str(val) for val in collect_dict[cluster]])))

                # # the overall tests
                groups = {popu: np.asarray(lengths, dtype=float) for popu, lengths in collect_dict.items()}

                _, kruskal_p = kruskal(*groups.values())

                fout.write("{}, Kruskal-Wallis test, P={}\n".format(chrom, round(kruskal_p, 3) if kruskal_p > 0.01 else format(kruskal_p, '.1e')))

                # if kruskal_p < 0.05:
                #
                #     # # using mwu test
                #     fout.write("\nMannwhitneyu pairwise comparisons:\n")
                #
                #     clusters = cluster_order
                #     pairs = list(combinations(clusters, 2))
                #
                #     for pop1, pop2 in pairs:
                #         x = np.asarray(collect_dict[pop1], dtype=float)
                #         y = np.asarray(collect_dict[pop2], dtype=float)
                #
                #         stat, muw_p = mannwhitneyu(x, y, alternative="two-sided", method="auto")
                #         muw_p_label = significance_label(muw_p)
                #         muw_p = round(muw_p, 3) if muw_p > 0.01 else format(muw_p, '.1e')
                #
                #         stat, ks_p = ks_2samp(x, y)
                #         ks_p_label = significance_label(ks_p)
                #         ks_p = round(ks_p, 3) if ks_p > 0.01 else format(ks_p, '.1e')
                #
                #         fout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chrom, pop1, pop2, muw_p, muw_p_label, ks_p, ks_p_label))
                #
                # fout.write("\nStats values:\n")

                for cluster in cluster_order:
                    if len(collect_dict[cluster]) == 0:
                        continue

                    q1 = int(np.percentile(collect_dict[cluster], 25))
                    median = int(np.percentile(collect_dict[cluster], 50))
                    q3 = int(np.percentile(collect_dict[cluster], 75))
                    iqr = q3 - q1
                    lower = int(q1 - 1.5 * iqr)
                    upper = int(q3 + 1.5 * iqr)
                    mean = int(np.average(collect_dict[cluster]))

                    p_values = []
                    for target_popu in cluster_order:
                        stat, muw_p = mannwhitneyu(collect_dict[cluster], collect_dict[target_popu], alternative="two-sided", method="auto")
                        muw_p_label = significance_label(muw_p)
                        muw_p = round(muw_p, 3) if muw_p > 0.01 else format(muw_p, '.1e')

                        p_values.append(str(muw_p))

                    fout.write("{}\t{}\t{}\t{}\t{}\n".format(chrom, cluster, len(collect_dict[cluster]), median, "\t".join(p_values)))


                # # plot R
                all_values, all_groups = [], []

                for cluster in cluster_order:
                    all_values.extend(collect_dict[cluster])
                    all_groups.extend([cluster] * len(collect_dict[cluster]))

                df = pd.DataFrame({"value": all_values, "group": all_groups})

                # # plot a boxplot
                fig, ax = plt.subplots(figsize=(3, 5))
                sns.boxplot(data=df, x="group", y="value", palette=cluster_colors, showfliers=False, ax=ax,saturation=1,)
                sns.despine(ax=ax)

                ax.set_xlabel("Cluster", fontsize=16)
                ax.set_ylabel("αSat length (bp)", fontsize=16)
                # 调整横轴 tick label 字体大小
                ax.tick_params(axis="x", labelsize=16)

                # 如果需要，也可以调整纵轴 tick label 字体大小
                ax.tick_params(axis="y", labelsize=16)

                # ax.yaxis.get_offset_text().set_visible(False)
                plt.tight_layout()
                plt.savefig(os.path.join(plot_output_folder, "res.{}_length.{}.boxplot.png".format(sat, chrom)))
                plt.savefig(os.path.join(plot_output_folder, "res.{}_length.{}.boxplot.svg".format(sat, chrom)))

                plt.close()

                # # plot a violin plot
                # df["comparison"] = "All"

                # fig, ax = plt.subplots(figsize=(3, 5))
                # sns.violinplot(data=df, x="group", y="value", cut=0, linewidth=1, color="#5975a5", ax=ax)
                # # sns.violinplot(data=df, x="comparison", y="value", hue="group", cut=0, linewidth=1, ax=ax, split=True)
                #
                # sns.despine(ax=ax)
                #
                # ax.set_xlabel("")
                # ax.set_ylabel("")
                # ax.yaxis.get_offset_text().set_visible(False)
                # plt.tight_layout()
                # plt.savefig(os.path.join(plot_output_folder, "res.{}_length.{}.violin.png".format(sat, chrom)))
                # plt.close()

            # fdr_values = multipletests(all_p, method="fdr_bh")[1]
            # print(fdr_values)




if __name__ == '__main__':

    annot_res_folder = "/mnt/e/workspace/cpc2_centro/res_all_annot"
    output_folder = "/mnt/e/workspace/cpc2_centro/res_all_post"

    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    with open(os.path.join(output_folder, "stats.whole_peri_cent.json"), 'r') as file:
        whole_peri_cent = dict(json.load(file))
    with open(os.path.join(output_folder, "stats.whole_core_cent.json"), 'r') as file:
        whole_core_cent = dict(json.load(file))

    # # STEP: load population info
    all_included_samples = []
    detail_popu_dict = {}
    detail_popu_order = ["AFR", "EUR", "AMR", "SAS", "EAS", "CHN"]

    hyper_popu_dict = {}
    hyper_popu_order = ["Other", "CHN"]
    with open("../data/data_population_hgsvc3_hprc2_cpc2.txt") as fin:
        for line in fin:
            line_split = line.strip().split("\t")

            all_included_samples.append(line_split[1])
            detail_popu_dict[line_split[1]] = line_split[6]
            hyper_popu_dict[line_split[1]] = line_split[5]

    # # STEP: stats aSat length
    stats_aSat_length_from_monomer(merge_thresh="0k", spec_aSat=None)
    # stats_aSat_length_from_monomer(merge_thresh="0k", spec_aSat="active")
    # stats_aSat_length_from_monomer(merge_thresh="0k", spec_aSat="inactive")

    # # STEP: stats oSat length
    stats_oSAT_lengths_from_rpmask(merge_thresh="0k")

    # # STEP: plot
    plot_SAT_length_among_popu(output_folder, "detail", spec_aSat=None)
    # plot_SAT_length_among_popu(output_folder, "detail", spec_aSat="active")
    # plot_SAT_length_among_popu(output_folder, "detail", spec_aSat="inactive")

    # # STEP generate a data matrix for all chrom of all samples
    genearte_SAT_length_matrix(output_folder, "detail", spec_aSat=None)

    plot_SAT_length_among_cluster(output_folder, os.path.join(output_folder, "principal_hor_analysis_outlier"), spec_aSat=None)
    # plot_SAT_length_among_cluster(output_folder, os.path.join(output_folder, "principal_hor_analysis_outlier"), spec_aSat="inactive")

    # # STEP: check correlation
    # check_length_corr_among_popu(output_folder, "stats.aSat_length_monomer.inactive.merge_0k.json", "stats.aSat_length_monomer.active.merge_0k.json", popu_mode="detail")
    #
    # check_length_corr_among_cluster(output_folder, "stats.aSat_length_monomer.inactive.merge_0k.json", "stats.aSat_length_monomer.active.merge_0k.json", os.path.join(output_folder, "principal_hor_analysis_outlier"))


