
"""
 This module has been created to analyze the mRNA data from nanostring nSolver 
 following  the guideline from user the manual at https://nanostring.com/wp-content/uploads/Gene_Expression_Data_Analysis_Guidelines.pdf
 and at https://tep.cancer.illinois.edu/files/2020/08/MAN-10030-03_nCounter_Advanced_Analysis_2.0_User_Manual-1.pdf
 It has following functions:

    1. read_rcc       > for reading ".RCC" files
    2. bg_correction  > for back ground correction using negative control probes 
    3. sample_qc      > for quality control of samples
    4. geNorm         > functions  to select house keeping genes for normalization using geNorm algorithm
    6. geCorr         > functions  to select house keeping genes for normalization using correlation based algorithm developed by author Dr. Rudramani Pokhrel
    5. normalize_data > for data normalization using house keeping genes
    6. gMean          > for calculation geometic mean from a list of numbers
    7. DE_analysis    > for differential analysis between two groups
    8. plot_volcano   > for ploting volcano plot

For more detailed of the functions see docstring: function.__doc__
"""


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.stats import linregress as  lr
import glob
import seaborn as sns
from scipy.stats import ttest_ind, ttest_rel
from statsmodels.stats.multitest import fdrcorrection
def read_rcc(path = "./", bc = "mean"):

    """
    Input:  
          path > path to the RCC files
          bc   > back ground correcion methods using negative control probes:
                            "mean"  >> for mean + 2 std 
                            "gmean" >> for geometric mean
                            "max"   >> for maximum

    Output: raw_count > raw count data
            count_bc  > back ground corrected count data
            QC        > quality control checking data
            PC        > positive control Limit of Detection
    """
    
    rcc = glob.glob(path+"*.RCC")

    # list of data frame of each samples
    L_df = []
    
    # List of background corrected data frame 
    L_df_bc = []

    # list QC 
    L_qc = []

    # list positive control concentration and counts to check if we are using accurate POS_F(=0.125) and POS_E(0.5)
    L_pc  = []

    for file in rcc:
        f = open(file, "r")
        lines = f.readlines()
        lines = [line.strip() for line in lines]
        ID = [line.split(",")[1] for line in lines if line.startswith("ID")][0]
        endogenuous = [line.split(",") for line in lines if line.startswith("Endogenous")]
        housekeeping = [line.split(",") for line in lines if line.startswith("Housekeeping")]

        # read endogenous mRNA
        genes = [i[1] for i in endogenuous]
        genes = [i.upper() for i in genes]
        value = [int(i[-1]) for i in endogenuous]
        gene_class = [i[0] for i in endogenuous]
        count = pd.DataFrame({"class": gene_class, ID: value}, index = genes)

        # read house keeping genes
        hs_genes = [i[1] for i in housekeeping]
        hs_genes = [i.upper() for i in hs_genes]
        hs_value = [int(i[-1]) for i in housekeeping]
        hs_class = [i[0] for i in housekeeping]
        count_hk = pd.DataFrame({"class": hs_class, ID: hs_value}, index = hs_genes)


        # concat mRNA and housekeeping genes
        count_df = pd.concat([count,count_hk])
        L_df.append(count_df)

        # Field of view QC
        fov_count = [float(line.split(",")[1]) for line in lines if line.startswith("FovCount")]
        fov = fov_count[1]/fov_count[0] *100

        # Binding Density
        bd = [float(line.split(",")[1]) for line in lines if line.startswith("BindingDensity")]

        # positive control Limit of Detection QC
        # Positive control QC
        pc = [line.split(",") for line in lines if line.startswith("Positive")]
        pc_conc = [float(i[1].split("(")[1].split(")")[0]) for i in pc ]
        pc_name = [i[1].split("(")[0] for i in pc ]
        pc_count = [float(i[3]) for i in pc ]
        pc_df = pd.DataFrame({"pc_conc": pc_conc, "pc_count": pc_count}, index= pc_name)
        pc_df["sample"] = [[ID] * pc_df.shape[0]][0]
        L_pc.append(pc_df)
        pc_df = pc_df.drop("POS_F")
        pos_e = pc_df.loc["POS_E"].pc_count


        # Negative control
        nc = [float(line.split(",")[-1]) for line in lines if line.startswith("Negative")]
        mean_nc = np.mean(nc)
        std_nc = np.std(nc)

        # Limit of detection (LOD) should be greater than 1
        LOD = pos_e/(mean_nc + 2*std_nc)

        max_nc = max(nc)
        max_nc = int(np.round(max_nc))

        nc = [i+0.001 for i in nc]
        gmean_nc = scipy.stats.gmean(nc)
        gmean_nc = int(np.round(gmean_nc))
        
        # let's do background correction
        value_bc = bg_correction(value, mean_nc, max_nc, std_nc, gmean_nc, bc = bc)
        count_bc = pd.DataFrame({"class": gene_class, ID: value_bc}, index = genes)
        
        hs_value_bc = bg_correction(hs_value, mean_nc, max_nc, std_nc, gmean_nc, bc = bc)
        count_hk_bc = pd.DataFrame({"class": hs_class, ID: hs_value_bc}, index = hs_genes)
        
        count_df_bc = pd.concat([count_bc,count_hk_bc])
        L_df_bc.append(count_df_bc)

        # r2 between pc_conc and pc_count should be greater than 0.95
        r = lr(np.log2(pc_df.pc_conc), y=np.log2(pc_df.pc_count), alternative='two-sided')[2]
        r2 = r**2

         # get the qc data
        qc_df = pd.DataFrame({"fov":fov,  "bd":bd, "LOD": LOD, "mean_nc": mean_nc,
                             "max_nc": max_nc, "std_nc": std_nc, "gmean_nc": gmean_nc, "R2": r2})
        qc_df.index = [ID]

        L_qc.append(qc_df)

    # get the raw count data   
    raw_count = pd.concat(L_df, axis = 1)
    raw_count = raw_count.loc[:,~raw_count.columns.duplicated()].copy()
    
    # get the background corrected count
    count_bc = pd.concat(L_df_bc, axis = 1)
    count_bc = count_bc.loc[:,~count_bc.columns.duplicated()].copy()

    # get positive control lines
    PC_df = pd.concat(L_pc, axis =0 )

    # get the qc matrix for each samples
    QC_df = pd.concat(L_qc)
    return raw_count, count_bc,  QC_df, PC_df



def bg_correction(value,  mean_nc, max_nc, std_nc, gmean_nc, bc = "mean"):

    if bc == "mean":
        mean_plus_2sd =  int(np.round(mean_nc + 2*std_nc))
        bc_list = [mean_plus_2sd if i < mean_plus_2sd else i for i in value]
        
    elif bc == "max":
        bc_list = [max_nc if i < max_nc else i for i in value]
        
    elif bc == "gmean":
         bc_list = [gmean_nc if i <gmean_nc else i for i in value]

    return bc_list



def sample_qc(data, fov_min = 75, bd_min = 0.1, bd_max = 2,  R2_min = 0.95, LOD_min = 1):

    """
    Input: data    > QC data from read_rcc ourput
           fov_min > minimum threshold for field of view, the default is 75
           bd_min  > minimum binding density, the default is 0.1
           bd_max  > maximum binding density, the default is 2
           R2_min  > R2 value of positive controls, the default is 0.95
           LOD_min > mimimum threashold of Limit of Detection  QC, default is 1
    Output: two outputs
          QC filtered data
          QC Unfilterd data
    """

    
    sample_before_qc = data.shape[0]
    #samples before qc
    sample_i = data.index
    df = data.copy()
    data = data[(data.fov > fov_min) & (data.bd < bd_max) & (data.bd > bd_min)\
         & (data.R2 > R2_min) & (data.LOD >1) ]
    sample_after_qc = data.shape[0]
    #samples before qc
    sample_f = data.index
    
    print(f"total number of samples befor QC: {sample_before_qc}")
    print(f"total number of samples after QC: {sample_after_qc}")
    if sample_before_qc != sample_after_qc:
        sample_removed = [i for i in sample_i if i not in sample_f]
        print('\ndata failed to pass QC is:')
        print(df.loc[sample_removed, ["fov", "bd", "LOD", "R2"]])
        
    return data, df


def geCorr(df, corr_cutoff = 0.75, number_required = 5):
    # select house keeping genes
    df = df.copy()
    df = df[df["class"] == "Housekeeping"]
    df = df.drop("class",axis = 1)
    corr = df.T.corr(method= "spearman")
    corr = corr[corr> corr_cutoff]
    
    thresh = corr.shape[0]-5
    corr = corr.dropna(thresh = 5, axis = 0)
    corr = corr.dropna(thresh = 5, axis = 1)
    corr = corr.dropna(axis = 0)
    corr = corr.dropna(axis = 1)
    corr_genes = [i for i in corr.index if i in corr.columns]
    mean = df.mean(axis = 1)
    std = df.std(axis = 1)
    df["mean"] = mean
    df["stdev"] = std
    df["cov"] = df.std(axis = 1)/ df.mean(axis = 1) *100
    hs_count_corr = df.loc[corr_genes]
    hs_count_corr = hs_count_corr.sort_values("cov")
    top5 = hs_count_corr.head(number_required)
    hk_genes = top5.index
    print("The table for mean and cov of seleted HK genes is:")
    print(df.loc[hk_genes][["mean","cov"]])
    print("\n The table for correlation coefficient of selected HK genes is:")
    print(corr.loc[hk_genes,hk_genes])
    if len(hk_genes)<5:
        print("seleted HK genes are less than 5. please increase the corr_cutoff.")
    return list(hk_genes)


    
def geNorm(df,  number_required = 5):
    df = df.copy()
    df = df[df["class"] == "Housekeeping"]
    df = df.drop("class",axis = 1)
    index = df.index

    a = 0

    # gene stability measure for each house keeping gene is M
    M = []
    while a < len(index):
        V = []
        for i in range(len(index)):
            if i !=a:

                # calculate log ratio of pair of genes a vs all other genes , and take std 
                v = np.log2(df.iloc[a]/df.iloc[i]).std()
                V.append(v)
        M.append(np.mean(V))
        a = a+1
    df["M_val"] = M
    df = df.sort_values("M_val")
    mean = df.mean(axis = 1)
    std = df.std(axis = 1)
    df["mean"] = mean
    df["stdev"] = std
    df["cov"] = df.std(axis = 1)/ df.mean(axis = 1) *100
    hk_genes = list(df.head(number_required).index)
    print("The table for mean and cov of seleted HK genes is:")
    print(df.loc[hk_genes][["mean","cov", "M_val"]])

    return(hk_genes)


def gMean(list1):
    mult = 1
    for i in list1:
        if i ==0:
            i = 0.001
            mult = mult*i
        else:
            mult = mult*i
    gmean = (mult)**(1/len(list1))
    return gmean
        



def normalize_data(data, method = "geNorm", number_required = 5, corr_cutoff = 0.75 ):

    """
    Two methods for normalizaiont

    "geNorm"  is done using selected house keeping genes
        1. Calculate the geometric mean of the selected housekeeping genes for each sample lane.
        2. Calculate the arithmetic mean of these geometric means for all sample lanes.
        3. Divide this arithmetic mean by the geometric mean of each lane to generate a lane-specific
        normalization factor.
        4. Multiply the counts for every gene by its lane-specific normalization factor.

    "geCorr" is done using 
        1. highly correlated  house keeping genes across samples  
        2. they also have least coefficient of variation across  a sample.

    Input: data            > bc_data or raw_count data
           method          > "geNorm" or "geCorr"
           number_required > desired number of house keeping genes, default is 5
           corr_cutoff     > minimum threshold for correlation if "geCorr" is used as a method default is 0.75

    """

    # select house keeping genes
    df = data[data["class"]=="Housekeeping"].drop("class", axis = 1)

    # select Endogenuous genes
    count_norm = data[data["class"]=="Endogenous"].drop("class", axis = 1)
    
    if method == "geNorm":
        selected_hk_genes = geNorm(df = data, number_required = number_required)
    else:
        selected_hk_genes = geCorr(df = data, number_required = number_required, corr_cutoff = corr_cutoff)
        
    df = df.loc[selected_hk_genes]

    # calculate gmean
    gmean = df.apply(gMean)

    # calculate arithmetic mean
    mean = gmean.mean()
    # normalization factor
    nf = mean/gmean
    
    # normalization QC
    nf = nf[(nf>0.1) & (nf < 10)]
    
    # ploting the normalization factor
    plt.figure(figsize = (10,5))
    sns.barplot( x = nf.index, y = nf.values)
    plt.xticks(rotation = 75, fontsize = 12)
    plt.yticks( fontsize = 15)
    plt.title("Normalization factor", fontsize = 20)
    plt.savefig("normalization_factor.pdf", bbox_inches = "tight")
    plt.show()
    print("Normalization factor should be within 0.1 to 10 scale.")
    print("Otherwise those samples were removed from normalied data.")
    
    count_norm = count_norm[nf.index]
    count_norm = count_norm*nf
    return count_norm


def DE_analysis(df , sample, sample_col, group1 , group0 ,  ttest = "ind"):

    """
    Input: df         > normalized count data with samples in columns and genes in row
           sample     > sample data
           sample_col > column from sample to do DE analysis, eg "treatment"
           group1     > numerator for log fold change example: "treated"
           group0     > denominator for log fold change example: "control"
           ttest      > "ind" >> independent wich is default
                        "pair" >> pairwise

    Output: data with logfold change, p_values and average counts over the groups
    """

    countT = df.T
    genes = countT.columns
    # let's add sample metadata infromation, in this data we have two sample, control (NV) and treatment (NS) mice grups
    # add back these information to count data
    countT[sample_col] = sample[sample_col].values
    countT_mean = countT.groupby(sample_col).mean()
    
    # calculate log2 fold change
    countT_lfc = np.log2(countT_mean.loc[group1].div(countT_mean.loc[group0]))
    
    # let's do independent sample t-test which is performed in log2 transfer data
    countT = df.T
    countT_log = np.log2(countT)
    countT_log[sample_col] = sample
    
    t_value = []
    p_value = []
    
    for col in countT_log.columns[0:-1]:
        df = countT_log[[col, sample_col]]
        first = df[df.treatment == group1]
        second = df[df.treatment == group0]
        if ttest == "ind":
            t, p = ttest_ind(first[col], second[col])
            t_value.append(t)
            p_value.append(p)

        else:
            t, p = ttest_ind(first[col], second[col])
            t_value.append(t)
            p_value.append(p)

        
    stat = pd.DataFrame({"lfc": countT_lfc, "t_value": t_value, "p_value":p_value})
    fdr = fdrcorrection(p_value, alpha=0.05, method='i', is_sorted=False)
    stat["fdr"] = fdr[-1]
    stat["Avg_" + group0] = list(countT_mean.T[group0].values)
    stat["Avg_" + group1] = list(countT_mean.T[group1].values)
    
    return stat




def plot_volcano(data):  

    """
    Input:  data > stat data from DE_analysis

    Output: Figure
    """

    de = []
    sig = []
    data = data.sort_values("lfc", ascending = False)
    for index, row in data.iterrows():
        if row["p_value"] < 0.05:
            if row["lfc"] < -1:
                de.append("lfc < -1")
            elif row["lfc"] > 1:
                de.append("lfc > 1")
            else:
                de.append("p < 0.05")
        else:
            de.append("p > 0.05")

    data["regulation"] = de   
    data["log10_pval"] = -np.log10(data["p_value"])
    data = data.sort_values("regulation")
    color = data.regulation.unique()
    pal = {"lfc < -1": "blue","lfc > 1": "red", "p < 0.05":"green", "p > 0.05":"black" }

    # find significant genes for text lebel
    data2 = data.copy()
    data2 = data2.sort_values("lfc", ascending = False)
    sig1 = list(data2[(data2.p_value < 0.05) &  (data.lfc> 1)].index)
    sig1  = list(sig1[0:3])
    sig2 = list(data2[(data2.p_value < 0.05) &  (data.lfc < -1)].index)
    sig2 = list(sig2[-3:])
    sig = sig1+sig2
    
    # let's plot
    plt.figure(figsize = (6,4))
    sns.scatterplot(data = data, x = "lfc", y = "log10_pval", hue = "regulation", palette = pal)
    plt.ylabel("-log10(pval)", fontsize = 15)
    plt.xlabel("log2(FC)", fontsize = 15)
    plt.xticks( fontsize = 15)
    plt.yticks( fontsize = 15)
    plt.legend(fontsize = 15, loc = (1.01,0.5))
    plt.axhline(y=-np.log10(0.05), color='green', linestyle='--', alpha = 0.5)
    plt.axvline(x=0, color='black', linestyle='--', alpha = 0.5)
    
    # let's put text of significant genes
    for i, gene in enumerate (sig):
        x = data.loc[gene].lfc
        y = data.loc[gene].log10_pval
        plt.text(x+0.05, y+0.05, gene, color = "black", fontsize = 12, fontstyle = "italic", fontweight = "bold")
    plt.savefig("volcano.pdf", bbox_inches = "tight")
    plt.show()



