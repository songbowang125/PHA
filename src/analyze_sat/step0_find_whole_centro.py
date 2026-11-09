import os
import json

import numpy as np


def find_whole_centro_for_one_sample(extension=5000000):
    chm13v2_chrom_length = {'chr1': 248387328, 'chr10': 134758134, 'chr11': 135127769, 'chr12': 133324548, 'chr13': 113566686, 'chr14': 101161492, 'chr15': 99753195, 'chr16': 96330374, 'chr17': 84276897, 'chr18': 80542538, 'chr19': 61707364, 'chr2': 242696752, 'chr20': 66210255, 'chr21': 45090682,
                            'chr22': 51324926, 'chr3': 201105948, 'chr4': 193574945, 'chr5': 182045439, 'chr6': 172126628, 'chr7': 160567428, 'chr8': 146259331, 'chr9': 150617247, 'chrX': 154259566, 'chrY': 62460029}

    # # load N file
    sample_nonN_dict = {}
    with open("../data/data_non-N_regions_hgsvc3_hprc2_cpc2.txt") as fin:

        for line in fin:

            project, sample_fasta, contig, start, end, length = line.strip().split()

            if sample_fasta not in sample_nonN_dict:
                sample_nonN_dict[sample_fasta] = {}

            if project == "hgsvc3":
                contig_split = contig.split("#")
                contig = "haplotype{}-{}_{}".format(contig_split[1], contig_split[2], contig_split[0])

            if contig not in sample_nonN_dict[sample_fasta]:
                sample_nonN_dict[sample_fasta][contig] = []
            sample_nonN_dict[sample_fasta][contig].append([int(start), int(end)])

    # # traverse each sample ragtag results
    ragtag_folder = os.path.join(work_folder, "chrall")

    for sample_ragtag_file in os.listdir(ragtag_folder):
        if not sample_ragtag_file.endswith("#ragtag.scaffold.agp"):
            continue

        sample = sample_ragtag_file.replace("#ragtag.scaffold.agp", "")
        if sample not in all_included_samples:
            continue
        # if "CPC"  in sample:
        #     continue
        print(sample)

        sample_ragtag_file = os.path.join(ragtag_folder, sample_ragtag_file)

        # # load mix and max sat cords
        sample_core_cent_region = {}    # # aSat with 1M extension
        sample_peri_cent_region = {}    # # aSat and oSat, with 1M extension

        for spec_chrom in chrom_order:
            chrom_input_folder = os.path.join(work_folder, spec_chrom)

            sample_aSat_file = os.path.join(chrom_input_folder, "{}#monomer.bed".format(sample))
            min_aSat_pos, max_aSat_pos = np.inf, -np.inf
            with open(sample_aSat_file) as fin:
                for line in fin:
                    line_split = line.strip().split()
                    contig_start, start, end = int(line_split[0].split(":")[1].split("-")[0]), int(line_split[1]), int(line_split[2])

                    start, end = start + contig_start, end + contig_start
                    if start < min_aSat_pos:
                        min_aSat_pos = start
                    if end > max_aSat_pos:
                        max_aSat_pos = end
            sample_core_cent_region[spec_chrom + "_RagTag"] = [max(1, min_aSat_pos - extension), max_aSat_pos + extension]

            # # load oSat
            sample_oSat_file = os.path.join(chrom_input_folder, "{}#osat.bed".format(sample))
            min_oSat_pos, max_oSat_pos = np.inf, -np.inf
            with open(sample_oSat_file) as fin:
                for line in fin:
                    line_split = line.strip().split()
                    contig_start, start, end = int(line_split[0].split(":")[1].split("-")[0]), int(line_split[1]), int(line_split[2])

                    start, end = start + contig_start, end + contig_start
                    if start < min_oSat_pos:
                        min_oSat_pos = start
                    if end > max_oSat_pos:
                        max_oSat_pos = end
            sample_peri_cent_region[spec_chrom + "_RagTag"] = [max(1, min_oSat_pos - extension), max_oSat_pos + extension]


        # # check
        candidate_contigs = []
        with open(sample_ragtag_file) as fin:

            for line in fin:
                line_split = line.strip().split()
                if len(line_split) != 9:
                    continue

                if line_split[6] == "scaffold":
                    continue

                candidate_contigs.append(line_split[5])

        for contig in sample_nonN_dict[sample]:
            if contig not in candidate_contigs:
                print("ERROR contig inconsistent", sample, contig)
                exit()

        # #
        chrom_include_contigs = {}
        with open(sample_ragtag_file) as fin:

            for line in fin:
                line_split = line.strip().split()
                if len(line_split) != 9:
                    continue

                if line_split[6] == "scaffold":
                    continue

                # # extract the mapped information for the contig
                chrom, start, end, contig = line_split[0], int(line_split[1]), int(line_split[2]), line_split[5]

                if chrom  not in sample_core_cent_region:
                    continue

                if chrom not in chrom_include_contigs:
                    chrom_include_contigs[chrom] = {"contig_cnt": 0, "contig_total_len": 0, "whole_peri_cent_flag": ".", "whole_core_cent_flag": "."}

                for nonN_start, nonN_end in sample_nonN_dict[sample][contig]:
                    chrom_include_contigs[chrom]["contig_cnt"] += 1
                    chrom_include_contigs[chrom]["contig_total_len"] += (nonN_end - nonN_start)

                    if nonN_start + start <= sample_core_cent_region[chrom][0] < sample_core_cent_region[chrom][1] <= nonN_end + start:
                        chrom_include_contigs[chrom]["whole_core_cent_flag"] = "Whole_Core_Centromere"

                    if nonN_start + start <= sample_peri_cent_region[chrom][0] < sample_peri_cent_region[chrom][1]  <= nonN_end + start:
                        chrom_include_contigs[chrom]["whole_peri_cent_flag"] = "Whole_Peri_Centromere"

        if extension == 5000000:
            sample_output_file = os.path.join(ragtag_folder, "{}#stats.whole_cent.5M.txt".format(sample))
        elif extension == 1000000:
            sample_output_file = os.path.join(ragtag_folder, "{}#stats.whole_cent.1M.txt".format(sample))
        elif extension == 10000:
            sample_output_file = os.path.join(ragtag_folder, "{}#stats.whole_cent.10K.txt".format(sample))
        elif extension == 0:
            sample_output_file = os.path.join(ragtag_folder, "{}#stats.whole_cent.0K.txt".format(sample))
        else:
            print("no this extension", extension)
            exit()

        with open(sample_output_file, "w") as fout:
            for chrom in chrom_include_contigs:
                chrom_no_ragtag = chrom.replace("_RagTag", "")

                if chrom_include_contigs[chrom]["contig_total_len"] >= chm13v2_chrom_length[chrom_no_ragtag] * 0.8:
                    chrom_include_contigs[chrom]["whole_chrom_flag"] = "Whole_chromosome"
                else:

                    chrom_include_contigs[chrom]["whole_chrom_flag"] = "."


                fout.write("{}\t{}\t{}\t{}\t{}\n".format(chrom, chrom_include_contigs[chrom]["contig_cnt"], chrom_include_contigs[chrom]["whole_chrom_flag"], chrom_include_contigs[chrom]["whole_core_cent_flag"], chrom_include_contigs[chrom]["whole_peri_cent_flag"]))

        # break

def stats_whole_centro_all_samples(out_path):

    if not os.path.exists(out_path):
        os.mkdir(out_path)

    ragtag_folder = os.path.join(work_folder, "chrall")

    with open(os.path.join(out_path, "stats.whole_core_cent.json"), "w") as whole_core_cent_fout, open(os.path.join(out_path, "stats.whole_peri_cent.json"), "w") as whole_peri_cent_fout:

        whole_core_cent = {}
        whole_peri_cent = {}

        for sample_stat_file in os.listdir(ragtag_folder):

            if not sample_stat_file.endswith("#stats.whole_cent.1M.txt"):
                continue
            sample = sample_stat_file.replace("#stats.whole_cent.1M.txt", "")
            sample_stat_file = os.path.join(ragtag_folder, sample_stat_file)

            if sample not in all_included_samples:
                continue

            whole_core_cent[sample] = []
            whole_peri_cent[sample] = []

            with open(sample_stat_file) as fin:
                for line in fin:
                    chrom = line.strip().split("\t")[0]
                    if "Whole_Core_Centromere" in line:
                        whole_core_cent[sample].append(chrom)
                    if "Whole_Peri_Centromere" in line:
                        whole_peri_cent[sample].append(chrom)

        json.dump(whole_core_cent, whole_core_cent_fout, indent=4)
        json.dump(whole_peri_cent, whole_peri_cent_fout, indent=4)



def merge_popu_and_whole_centro():

    with open(os.path.join(output_folder, "stats.whole_core_cent.json"), 'r') as file:
            whole_cent_dict = dict(json.load(file))

    with open("../data/data_population_hgsvc3_hprc2_cpc2.txt") as fin:
        for line in fin:
            line_split = line.strip().split("\t")

            project, sample, popu = line_split[0], line_split[1], line_split[6]

            sample_chroms = [chrom.replace("_RagTag", "") for chrom in whole_cent_dict[sample]]

            # sample_chroms = sorted(sample_chroms, key=lambda x: int(x[3:]) if x[3:].isdigit() else float('inf') if x == 'chrX' else float('inf'))
            print("{}\t{}\t{}\t{}".format(project, sample, popu, ",".join(sample_chroms)))

def stats_whole_centro_among_popu(cent_mode, popu_mode, CPC2_only=False):

    print(cent_mode, popu_mode, sep="\t")

    if cent_mode == "peri":
        with open(os.path.join(output_folder, "stats.whole_peri_cent.json"), 'r') as file:
            whole_cent_dict = dict(json.load(file))
    elif cent_mode == "core":
        with open(os.path.join(output_folder, "stats.whole_core_cent.json"), 'r') as file:
            whole_cent_dict = dict(json.load(file))
    else:
        print("No this cent mode", cent_mode)
        exit()


    if popu_mode == "hyper":
        popu_dict = hyper_popu_dict
        popu_order = hyper_popu_order
    elif popu_mode == "detail":
        popu_dict = detail_popu_dict
        popu_order = detail_popu_order
    else:
        print("No this popu_mode", popu_mode)
        exit()

    # # init
    pass_dict = {}
    all_dict = {}
    for chrom in chrom_order:
        pass_dict[chrom] = {}
        all_dict[chrom] = {}

        for popu in ["all"] + popu_order:
            pass_dict[chrom][popu] = 0
            all_dict[chrom][popu] = 0


    # # traverse
    for sample in whole_cent_dict:

        popu = popu_dict[sample]

        if CPC2_only:
            if popu == "CHN" and "CPCERZ2" not in sample:
                continue

        # # add to all dict
        for chrom in chrom_order:
            all_dict[chrom]["all"] += 1
            all_dict[chrom][popu] += 1

        # # add to pass dict
        for chrom_with_tag in whole_cent_dict[sample]:
            chrom = chrom_with_tag.replace("_RagTag", "")

            pass_dict[chrom]["all"] += 1
            pass_dict[chrom][popu] += 1

    # # output count
    print("chrom\tall\t{}".format("\t".join(popu_order)))
    for chrom in all_dict:
        print(chrom, end="\t")
        for popu in all_dict[chrom]:
            print(pass_dict[chrom][popu], end="\t")
        print()

    # # output freq
    print("chrom\tall\t{}".format("\t".join(popu_order)))
    for chrom in all_dict:
        print(chrom, end="\t")
        for popu in all_dict[chrom]:
            print(round(pass_dict[chrom][popu] / all_dict[chrom][popu], 3), end="\t")
        print()


    #
    # print(all_dict)
    #
    # print(pass_dict)


def stats_whole_centro_ratio_among_projects(work_path, centro_mode):


    if centro_mode == "core":
        with open(os.path.join(work_path, "stats.whole_core_cent.json"), 'r') as file:
            whole_cent = dict(json.load(file))
    elif centro_mode == "peri":
        with open(os.path.join(work_path, "stats.whole_peri_cent.json"), 'r') as file:
            whole_cent = dict(json.load(file))
    else:
        print("No this mode", centro_mode)
        exit()


    project_count = {"HGSVC3": 0, "HPRC2": 0, "CPC2": 0}
    project_completeness_count = {"HGSVC3": 0, "HPRC2": 0, "CPC2": 0}
    project_completeness_count_chroms = {"HGSVC3": {chrom: 0 for chrom in chrom_order}, "HPRC2": {chrom: 0 for chrom in chrom_order}, "CPC2": {chrom: 0 for chrom in chrom_order}}

    for sample in whole_cent:
        if "T2T-CHM13" in sample:
            continue

        if "vrk-ps-sseq" in sample:
            sample_project = "HGSVC3"
        elif "CPCERZ" in sample:
            sample_project = "CPC2"
        elif "hprc_r2" in sample or  "genbank" in sample:
            sample_project = "HPRC2"
        else:
            continue

        project_count[sample_project] += 1


        for chrom in whole_cent[sample]:
            # if "chrY" in chrom:
            #     continue

            project_completeness_count[sample_project] += 1
            project_completeness_count_chroms[sample_project][chrom.replace("_RagTag", "")] += 1


    print(project_count)
    print(project_completeness_count)

    project_completeness_ratio = {"HGSVC3": 0, "HPRC2": 0, "CPC2": 0}
    project_completeness_chroms_ratio = {"HGSVC3": {chrom: 0 for chrom in chrom_order}, "HPRC2": {chrom: 0 for chrom in chrom_order}, "CPC2": {chrom: 0 for chrom in chrom_order}}

    for project in project_completeness_ratio:
        project_completeness_ratio[project] = round(project_completeness_count[project] / (project_count[project] * 24), 3)

        for chrom in chrom_order:
            project_completeness_chroms_ratio[project][chrom] = round(project_completeness_count_chroms[project][chrom] / project_count[project], 3)

    print(project_completeness_ratio)
    # print(project_completeness_chroms_ratio)

    # # output

    print("count")
    for project in project_completeness_count_chroms:
        print(project, end="\t")
        for chrom in chrom_order:
            print(project_completeness_count_chroms[project][chrom], end="\t")

        print()

    print("ratio")
    for project in project_completeness_chroms_ratio:
        print(project, end="\t")
        for chrom in chrom_order:
            print(project_completeness_chroms_ratio[project][chrom], end="\t")

        print()



if __name__ == '__main__':

    work_folder = "/mnt/e/workspace/cpc2_centro/res_all_annot"
    output_folder = "/mnt/e/workspace/cpc2_centro/res_all_post"


    chrom_order = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]
    # chrom_order = [ "chr14"]


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


    # # STEP: the centromere completeness for each sample
    find_whole_centro_for_one_sample(extension=0)

    # # STEP: collect the result for each sample and output a json
    stats_whole_centro_all_samples(output_folder)

    # # STEP: compare among projects
    stats_whole_centro_ratio_among_projects(output_folder, "core")

    # # STEP: merge population and whole centro
    merge_popu_and_whole_centro()

    # # STEP: stats the whole ratio among chroms and among populations
    stats_whole_centro_among_popu("core", "detail", CPC2_only=False)

