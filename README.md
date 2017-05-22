# PMI-associated-gene  

Matlab Codes and sample data for our paper: "Systematic analysis of gene expression patterns associated with postmortem interval in human tissues"  

-----------------------------------
Guidance:

1. Make sure your version of MATLAB software is at least 2013.  

2. Edit file paths (e.g., the HOMEPATH) in the scripts based on your settings.  

3. Unzip the data.zip. The folder contains all necessary sample data (Adipose_Subcutaneous) for reproducing the results (part 1 and 4) in our paper.  

  -  genelist.mat --- 18763 Protein Coding Genes
  -  expr_Adipose_Subcutaneous_159_144.mat --- Gene expression matrix
  -  cov_Adipose_Subcutaneous_159_144.mat --- Covariates
  -  peer_Adipose_Subcutaneous_159_144.mat --- PEER algorithm to infer the hidden data structures
  -  smpl_Adipose_Subcutaneous_159_144.mat --- Subject-level variables, including age, gender, BMI, etc., and sample-level variable (PMI)

4. Directly run linear_res/s1_res_regr.m first and then s2_extract_sig_fdr_gene_list.m to identify PMI-associated genes (part 1, results).

5. Directly run PMI_dispersion/s1_levene_test_p_values.m first and then s2_extract_fdr_genes.m to identify PMI-associated DV genes (part 4, results).


Note that we did not attach the genotype data in the PMI-by-genotype analyses, due to the data privacy of the GTEx project. All the data is available upon request from dbGaP (https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000424.v6.p1).

-----------------------------------
Contact: zhuyizhang@bjmu.edu.cn  

