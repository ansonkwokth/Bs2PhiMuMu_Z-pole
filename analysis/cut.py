import uproot
import pandas as pd
import numpy as np
import sys
import glob

folder_name = sys.argv[1]

file_name_dt = {"sig.":     folder_name + "ee2Z2Bs2PhiMuMu_reco.root",
                "bb":    folder_name + "ee2Z2b_comb_cutted_reco.root",
                "cc":    folder_name + "ee2Z2c_comb_cutted_reco.root",
                }

print("Loading data")
file_dt = {k: uproot.open(i) for k, i in file_name_dt.items()}

df_dt = {k: pd.DataFrame(np.array(file['t']['features'].array())) for k, file in file_dt.items()}

print("Finished loading data")





def cut_chi2(df):
    n0 = len(df)
    if (n0 == 0):
        return df, 0, 0
    else:
        df = df[df['Chi2'] < 15]
        n1 = len(df)
        return df, n1/n0, n1

def cut_ps(df):
    n0 = len(df)
    if (n0 == 0):
        return df, 0, 0
    else:
        df = df[df['PKp'] > 2]
        df = df[df['PKm'] > 2]
        df = df[df['Pmup'] > 2]
        df = df[df['Pmum'] > 2]
        n1 = len(df)
        return df, n1/n0, n1


def cut_mBs(df):
    n0 = len(df)
    if (n0 == 0):
        return df, 0, 0
    else:
        df = df[(df['mBs'] > 5.3663 - 0.5) & (df['mBs'] < 5.3663 + 0.5)]
        n1 = len(df)
        return df, n1/n0, n1


def cut_mPhi(df):
    n0 = len(df)
    if (n0 == 0):
        return df, 0, 0
    else:
        df = df[(df['mPhi'] > 1) & (df['mPhi'] < 1.04)]
        n1 = len(df)
        return df, n1/n0, n1


def cut_mDimu(df):
    n0 = len(df)
    if (n0 == 0):
        return df, 0, 0
    else:
        #df = df[(df['mDimu'] < 1) | (df['mDimu'] > 1.04)]
        #df = df[(df['mDimu'] < 3.00) | (df['mDimu'] > 3.2)]
        #df = df[(df['mDimu'] < 3.6) | (df['mDimu'] > 3.8)]
        df = df[(df['mDimu'] < 1**0.5) | (df['mDimu'] > 1.08**0.5)]
        df = df[(df['mDimu'] < 9**0.5) | (df['mDimu'] > 10.2**0.5)]
        df = df[(df['mDimu'] < 13**0.5) | (df['mDimu'] > 14.4**0.5)]
        n1 = len(df)
        return df, n1/n0, n1
    


fn_dt = {
        'chi2': cut_chi2,
        'ps': cut_ps,
        'mPhi': cut_mPhi,
        'mDimu': cut_mDimu,
#        'mBs': cut_mBs
        }

cut_eff_lt = []
raw_lt = []
df_cut_dt = {}
for k, df in df_dt.items():
    cut_eff_i = []
    raw_i = []

    for ck, cfn in fn_dt.items():
        df, eff, n1 = fn_dt[ck](df)
        cut_eff_i.append(eff)
        raw_i.append(n1)

    df_cut_dt[k] = df
    cut_eff_lt.append(cut_eff_i)
    raw_lt.append(raw_i)


df_cut = pd.DataFrame(cut_eff_lt, columns=fn_dt.keys(), index=file_name_dt.keys())
df_cut_raw = pd.DataFrame(raw_lt, columns=fn_dt.keys(), index=file_name_dt.keys())

print()
print(df_cut.T)
print()
print(df_cut_raw.T)



#target_folder_name = '../data/preselect/'
#target_folder_name = '../data/preselect/'
#target_folder_name = '../data/preselect_chi2/'
target_folder_name = '../data/preselect_q2/'

for k, df in df_cut_dt.items():
    #df.to_csv(target_file_name_dt[k], index=False)
    outName = file_name_dt[k].replace(folder_name, target_folder_name)
    outName = outName.replace("reco", "preselect")
    print("saving:", outName) 
    outFile = uproot.recreate(outName)
    outFile['t'] = df





