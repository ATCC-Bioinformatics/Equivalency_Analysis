#Collection of scripts involved in the assessment of small variants (SNPs/Indels) between refseq references of ATCC strains relative to the ATCC genome portal strain.

run_referenceAssembly_compile_data.py will run ReferenceAssembly.py in order to generate small variants between ACC read data and related RefSeq reference data and provide annotations such as synonymous/nonsynonymous mutation.

#File paths:#

Input files:

/home/adams/ADAMS/Applicate-Note/equiv_analysis_refAssembly/atcc_product_in_RefSeq.csv 

/home/adams/ADAMS/Applicate-Note/equiv_analysis_refAssembly/June2021/june2021_gcf_list.csv 

Raw reads: /home/shared/sequencing-data ReferenceAssembly path: /mnt/sd1/pipeline/AMGP-reference-based-assembly/ReferenceAssembly.py 

Output files: All GCFs were saved to: /home/adams/ADAMS/Application-Note/equiv_analysis_refAssembly/public_assemblies 

All GBKS were saved to: /home/adams/ADAMS/Application-Note/equiv_analysis_refAssembly/gbk

#Tools:# - Python 3.6.8 - Qualimap v.2.2.1 - BWA 0.7.17 - samtools 1.9 - GATK 4.0.8.1 - Requires Java - VEP 95.1 - tabix 1.10.2 - bcftools 1.9