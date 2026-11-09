"""

python step0_ragtag.py -i /home/data/user/songbo/workspace/cpc2_centro/data/part2 -o /home/data/user/songbo/workspace/cpc2_centro/res_cpc2 -r /home/data/user/songbo/workspace/cpc2_centro/data_ref/T2T-CHM13.fa.gz -t 16 -n cn02 cn04 cn05 cn06 -c slurm
python step1_coarse_aSAT.py -i  /home/data/user/songbo/workspace/cpc2_centro/res_cpc2 -s /home/data/user/songbo/workspace/cpc2_centro/pack_aSAT_annot/ -t 16 -n cn02 cn03 cn04 cn05 cn06 -c slurm
python step2_fine_aSAT.py -i  /home/data/user/songbo/workspace/cpc2_centro/res_cpc2 -s /home/data/user/songbo/workspace/cpc2_centro/pack_aSAT_annot/  -t 4 -n cn02 cn03 cn04 cn05 cn06 -c slurm

python step3_hor_aSAT.py -i  /home/data/user/songbo/workspace/cpc2_centro/res_cpc2 -s /home/data/user/songbo/workspace/cpc2_centro/pack_aSAT_annot/ -t 8 -n cn02 cn03 cn04 cn05 cn06 -c slurm
python step4_run_oSAT.py -i  /home/data/user/songbo/workspace/cpc2_centro/res_cpc2 -s /home/data/user/songbo/workspace/cpc2_centro/pack_aSAT_annot/ -t 8 -n cn02 cn03 cn04 cn05 cn06 -c slurm

python step5_proc_oSAT.py -i  /home/data/user/songbo/workspace/cpc2_centro/res_cpc2
python step6_final.py -i  /home/data/user/songbo/workspace/cpc2_centro/res_cpc2 -o /home/data/user/songbo/workspace/cpc2_centro/res_cpc2_final

python step7_to_local.py
python step8_pq_flanking.py

"""