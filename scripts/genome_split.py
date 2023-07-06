# Not originally created for public use and will not be updated further

import pandas as pd

FILE_NAME = 'UPDATE WITH FILE PATH'
RESULT_FILE_NAME = 'UPDATE WITH RESULT FILE'


def create_dataframe(df):
    return pd.DataFrame({'Position': df[0].tolist(), 'A': df[1].tolist(), 'G': df[2].tolist(),
                         'C': df[3].tolist(), 'U': df[4].tolist(), 'GAP': df[5].tolist()})


whole_genome = pd.read_excel(FILE_NAME, skiprows=0, index_col=False)
ORF1ab = pd.read_excel(FILE_NAME, skiprows=266,
                       nrows=21290, header=None, index_col=False)
S = pd.read_excel(FILE_NAME, skiprows=21563, nrows=3822,
                  header=None, index_col=False)
ORF3ab = pd.read_excel(FILE_NAME, skiprows=25393,
                       nrows=828, header=None, index_col=False)
E = pd.read_excel(FILE_NAME, skiprows=26245, nrows=228, header=None, index_col=False)
M = pd.read_excel(FILE_NAME, skiprows=26523, nrows=669, header=None, index_col=False)
ORF6 = pd.read_excel(FILE_NAME, skiprows=27202,
                     nrows=186, header=None, index_col=False)
ORF7ab = pd.read_excel(FILE_NAME, skiprows=27394,
                       nrows=494, header=None, index_col=False)
ORF8 = pd.read_excel(FILE_NAME, skiprows=27894,
                     nrows=366, header=None, index_col=False)
N = pd.read_excel(FILE_NAME, skiprows=28274, nrows=1260,
                  header=None, index_col=False)
ORF10 = pd.read_excel(FILE_NAME, skiprows=29558, nrows=117,
                      header=None, index_col=False)

ORF1ab_df = create_dataframe(ORF1ab)
S_df = create_dataframe(S)
ORF3ab_df = create_dataframe(ORF3ab)
E_df = create_dataframe(E)
M_df = create_dataframe(M)
ORF6_df = create_dataframe(ORF6)
ORF7ab_df = create_dataframe(ORF7ab)
ORF8_df = create_dataframe(ORF8)
N_df = create_dataframe(N)
ORF10_df = create_dataframe(ORF10)

writer = pd.ExcelWriter(RESULT_FILE_NAME)

whole_genome.to_excel(writer, sheet_name='Whole_genome', index=False)
ORF1ab_df.to_excel(writer, sheet_name='ORF1ab', index=False)
S_df.to_excel(writer, sheet_name='S', index=False)
ORF3ab_df.to_excel(writer, sheet_name='ORF3ab', index=False)
E_df.to_excel(writer, sheet_name='E', index=False)
M_df.to_excel(writer, sheet_name='M', index=False)
ORF6_df.to_excel(writer, sheet_name='ORF6', index=False)
ORF7ab_df.to_excel(writer, sheet_name='ORF7ab', index=False)
ORF8_df.to_excel(writer, sheet_name='ORF8', index=False)
N_df.to_excel(writer, sheet_name='N', index=False)
ORF10_df.to_excel(writer, sheet_name='ORF10', index=False)

writer.save()
