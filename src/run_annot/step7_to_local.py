import os
import json

if __name__ == '__main__':

    projects = ["cpc2", "ape", "hgsvc3", "hprc2", "ref"]

    for project in projects:

        input_folder = "/home/data/user/songbo/workspace/cpc2_centro/res_{}".format(project)

        output_folder = "/home/data/user/songbo/workspace/cpc2_centro/res_all_annot"

        if not os.path.exists(output_folder):
            os.mkdir(output_folder)


        chrom_order = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]

        for sample in os.listdir(input_folder):
            sample_folder = os.path.join(input_folder, sample)

            chrall_output_folder = os.path.join(output_folder, "chrall")
            if not os.path.exists(chrall_output_folder):
                os.mkdir(chrall_output_folder)

            os.system("cp {}/res_ragtag/ragtag.scaffold.agp {}/{}#ragtag.scaffold.agp".format(sample_folder, chrall_output_folder, sample))
            os.system("cp {}/aSAT.hor.stv_stats.tsv {}/{}#stv_stats.tsv".format(sample_folder, chrall_output_folder, sample))

            for chrom in chrom_order:
                print(sample, chrom)

                chrom_output_folder = os.path.join(output_folder, chrom)
                if not os.path.exists(chrom_output_folder):
                    os.mkdir(chrom_output_folder)


                sample_output_prefix = os.path.join(chrom_output_folder, "{}#".format(sample))

                hor_bed = "{}/aSAT.hor.stv_raw.bed".format(sample_folder)
                monomer_bed = "{}/aSAT.monomer.bed".format(sample_folder)
                sf_monomer_bed = "{}/aSAT.sf.bed".format(sample_folder)
                osat_bed = "{}/oSAT.bed".format(sample_folder)
                rpmask_bed = "{}/rpmask.out".format(sample_folder)


                cmd = "grep {}_RagTag {} | awk -v OFS='\t' ".format(chrom, monomer_bed) + " '{print $1, $2, $3, $4, $6}' " + " > {}monomer.bed".format(sample_output_prefix)
                os.system(cmd)

                cmd = "grep {}_RagTag {} | awk -v OFS='\t' ".format(chrom, hor_bed) + " '{print $1, $2, $3, $4, $6}' " + " > {}stv_hor.bed".format(sample_output_prefix)
                os.system(cmd)

                cmd = "grep {}_RagTag {} | awk -v OFS='\t' ".format(chrom, sf_monomer_bed) + " '{print $1, $2, $3, $4, $6}' " + " > {}sf.bed".format(sample_output_prefix)
                os.system(cmd)

                cmd = "grep {}_RagTag {} | awk -v OFS='\t' ".format(chrom, osat_bed) + " '{print $1, $2, $3, $4, $5}' " + " > {}osat.bed".format(sample_output_prefix)
                os.system(cmd)

                cmd = "grep {}_RagTag {} | awk -v OFS='\t' ".format(chrom, rpmask_bed) + " '{print $5, $6, $7, $11, $9, $10}' " + " > {}rpmask.bed".format(sample_output_prefix)
                os.system(cmd)






