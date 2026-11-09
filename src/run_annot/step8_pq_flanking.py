import os
import pysam
import multiprocessing

def reverse_complement_seq(seq):
    """
    given a seq, get its reverse complementary seq
    """
    inv_seq = ""
    for i in range(len(seq) - 1, -1, -1):
        bp = seq[i]

        if bp == 'A' or bp == "a":
            inv_bp = 'T'
        elif bp == 'T' or bp == "t":
            inv_bp = 'A'
        elif bp == 'C' or bp == "c":
            inv_bp = 'G'
        elif bp == 'G' or bp == "g":
            inv_bp = 'C'
        else:
            inv_bp = 'N'

        inv_seq += inv_bp

    return inv_seq


def run_minimap2(chrom, ref_fa, alt_fa, output_paf, output_fasta):

    os.system("minimap2 -t 1 -x asm20 --eqx -Y --secondary=no {} {} > {}".format(ref_fa, alt_fa, output_paf))

    max_align_length = -1
    max_align_contig, max_align_start_on_alt, max_align_end_on_alt, max_align_orient_on_alt = None, None, None, None
    max_align_start_on_ref, max_align_end_on_ref = None, None

    with open(output_paf) as fin:
        for line in fin:
            line_split = line.strip().split()

            alt_contig, alt_align_start, alt_align_end, alt_align_orient, alt_align_length = line_split[0], int(line_split[2]), int(line_split[3]), line_split[4], int(line_split[9])

            alt_contig_shift = int(alt_contig.split(":")[1].split("-")[0])

            if alt_align_length > max_align_length:
                max_align_contig, max_align_start_on_alt, max_align_end_on_alt, max_align_orient_on_alt = alt_contig, alt_align_start, alt_align_end, alt_align_orient
                max_align_start_on_ref, max_align_end_on_ref = alt_contig_shift + alt_align_start, alt_contig_shift + alt_align_end

                max_align_length = alt_align_length

    if max_align_start_on_alt is not None:
        print(max_align_contig, max_align_start_on_alt, max_align_end_on_alt, max_align_orient_on_alt)
        with open(output_fasta, "w") as fout, pysam.FastaFile(alt_fa) as fin:
            fout.write(">{}_RagTag:{}-{}-{}\n".format(chrom, max_align_start_on_ref, max_align_end_on_ref, max_align_orient_on_alt))

            if max_align_orient_on_alt == "+":
                fout.write(fin.fetch(max_align_contig, max_align_start_on_alt, max_align_end_on_alt) + "\n")
            else:
                fout.write(reverse_complement_seq(fin.fetch(max_align_contig, max_align_start_on_alt, max_align_end_on_alt)) + "\n")



if __name__ == '__main__':

    projects = ["cpc2", "ape", "hgsvc3", "hprc2", "apg1", "apr", "sea3k1", "tbn"]

    input_folder = "/home/data/user/songbo/workspace/cpc2_centro/"

    ref_folder = "/home/data/user/songbo/workspace/cpc2_centro/res_ref/T2T-CHM13.fa.gz"

    output_folder = "/home/data/user/songbo/workspace/cpc2_centro/res_all_flanking_50k"
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    chrom_order = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]

    # # load flanking bed
    chrom_flanking_dict = {}
    with open("chm13_v2.0_monomeric_p-qarm_only.50k.bed") as fin:
        for line in fin:
            line_split = line.strip().split()

            chrom, start, end, arm = line_split[0], int(line_split[1]), int(line_split[2]), line_split[3]

            if arm not in chrom_flanking_dict:
                chrom_flanking_dict[arm] = {}
            chrom_flanking_dict[arm][chrom] = [start, end]

    # # for ref
    for arm in chrom_flanking_dict:
        arm_output_folder = os.path.join(output_folder, arm)
        if not os.path.exists(arm_output_folder):
            os.mkdir(arm_output_folder)

        for chrom in chrom_flanking_dict[arm]:
            chrom_output_folder = os.path.join(arm_output_folder, chrom)
            if not os.path.exists(chrom_output_folder):
                os.mkdir(chrom_output_folder)

            chrom_peri_fasta = os.path.join(ref_folder, "res_fine_aSAT/{}_RagTag.peri_cent.fa".format(chrom))

            with pysam.FastaFile(chrom_peri_fasta) as fin, open(os.path.join(chrom_output_folder, "{}.flanking.fa".format("T2T-CHM13.fa.gz")), "w") as fout:
                chrom_contig_name = fin.references[0]

                chrom_contig_start = int(chrom_contig_name.split(":")[1].split("-")[0])

                chrom_flanking_start, chrom_flanking_end = chrom_flanking_dict[arm][chrom][0], chrom_flanking_dict[arm][chrom][1]

                if chrom_flanking_start < chrom_contig_start:
                    print("ERROR chrom contig start", arm, chrom, chrom_contig_start, chrom_flanking_dict[arm][chrom])
                    exit()

                fout.write(">{}_RagTag:{}-{}\n".format(chrom, chrom_flanking_start, chrom_flanking_end))
                fout.write(fin.fetch(chrom_contig_name, chrom_flanking_start - chrom_contig_start, chrom_flanking_end - chrom_contig_start) + "\n")

    # # for each file
    process_pool = multiprocessing.Pool(48)

    for arm in chrom_flanking_dict:
        arm_output_folder = os.path.join(output_folder, arm)

        for chrom in chrom_flanking_dict[arm]:
            chrom_output_folder = os.path.join(arm_output_folder, chrom)

            for project in projects:
                project_input_folder = os.path.join(input_folder, "res_{}".format(project))

                for sample in os.listdir(project_input_folder):
                    chrom_peri_fasta = os.path.join(project_input_folder, "{}/res_fine_aSAT/{}_RagTag.peri_cent.fa".format(sample, chrom))
                    chrom_ref_fasta = os.path.join(chrom_output_folder, "{}.flanking.fa".format("T2T-CHM13.fa.gz"))
                    chrom_paf_file = os.path.join(chrom_output_folder, "{}.paf".format(sample))
                    chrom_alt_fasta = os.path.join(chrom_output_folder, "{}.flanking.fa".format(sample))

                    process_pool.apply_async(run_minimap2, (chrom, chrom_ref_fasta, chrom_peri_fasta, chrom_paf_file, chrom_alt_fasta))
                    # run_minimap2(chrom, chrom_ref_fasta, chrom_peri_fasta, chrom_paf_file, chrom_alt_fasta)

    process_pool.close()
    process_pool.join()




