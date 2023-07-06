# Not originally created for public use and will not be updated further

import numpy as np
from scipy.stats import chi2_contingency
from statsmodels.sandbox.stats.multicomp import multipletests
import pandas as pd


def checkEqual(lst):
    ele = lst[0]
    chk = True

    # Comparing each element with first item
    for item in lst:
        if ele != item:
            chk = False
            break

    if chk:
        return True
    else:
        return False

df = pd.read_csv('UPDATE WITH FILE PATH')

col_names = ["ORF1a","ORF1b", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF10", "Whole Genome"]

count = 1
i = 0
result_columns = []
while i < 13:
    res_col = []

    col_1 = (df.iloc[:, count]).tolist()
    col_2 = (df.iloc[:, count+1]).tolist()

    ab_1 = [col_1[0], col_1[1]]
    ab_2 = [col_2[0], col_2[1]]

    ac_1 = [col_1[0], col_1[2]]
    ac_2 = [col_2[0], col_2[2]]

    bc_1 = [col_1[1], col_1[2]]
    bc_2 = [col_2[1], col_2[2]]

    if checkEqual(col_1) or checkEqual(col_2):
        res_col.append(1)
    else:
        obs = np.array([col_1, col_2])
        res_col.append(chi2_contingency(obs).pvalue)

    if checkEqual(ab_1) or checkEqual(ab_2):
        res_col.append(1)
    else:
        obs = np.array([ab_1, ab_2])
        res_col.append(chi2_contingency(obs).pvalue)

    if checkEqual(ac_1) or checkEqual(ac_2):
        res_col.append(1)
    else:
        obs = np.array([ac_1, ac_2])
        res_col.append(chi2_contingency(obs).pvalue)

    if checkEqual(bc_1) or checkEqual(bc_2):
        res_col.append(1)
    else:
        obs = np.array([bc_1, bc_2])
        res_col.append(chi2_contingency(obs).pvalue)

    result_columns.append(res_col)
    count += 2
    i += 1

pvalue_df = pd.DataFrame()

i = 0
while i < 13:
    pvalue_df[col_names[i]] = result_columns[i]
    i += 1

cor_1 = multipletests(pvalue_df.iloc[0].tolist(), alpha=.05, method='bonferroni')
pvalue_df.loc[len(pvalue_df)] = cor_1[1].tolist()

cor_2 = multipletests(pvalue_df.iloc[1].tolist(), alpha=.05, method='bonferroni')
pvalue_df.loc[len(pvalue_df)] = cor_2[1].tolist()

cor_3 = multipletests(pvalue_df.iloc[2].tolist(), alpha=.05, method='bonferroni')
pvalue_df.loc[len(pvalue_df)] = cor_3[1].tolist()

cor_4 = multipletests(pvalue_df.iloc[3].tolist(), alpha=.05, method='bonferroni')
pvalue_df.loc[len(pvalue_df)] = cor_4[1].tolist()

print(pvalue_df)
pvalue_df.to_csv("Test_Results.csv")