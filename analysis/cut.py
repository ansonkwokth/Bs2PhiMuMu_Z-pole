import uproot
import pandas as pd
import numpy as np
import sys
import glob

folder_name = sys.argv[1]

file_name_dt = {"sig.":     folder_name + "ee2Z2Bs2PhiMuMu_reco.root",
                "comb.":    folder_name + "ee2Z2b_comb_reco.root",
                "Jpsi":     folder_name + "ee2Z2Bs2PhiJpsi_reco.root",
                "psi2S":    folder_name + "ee2Z2Bs2PhiPsi_reco.root",
                "phi1020":  folder_name + "ee2Z2Bs2PhiPhi_reco.root"}

print("Loading data")
file_dt = {k: uproot.open(i) for k, i in file_name_dt.items()}


df_dt = {k: pd.DataFrame(np.array(file['t']['features'].array())) for k, file in file_dt.items()}

print("Finished loading data")





def cut_chi2(df):
    n0 = len(df)
    if (n0 == 0):
        return df, 0
    else:
        df = df[df['Chi2'] < 15]
        n1 = len(df)
        return df, n1/n0


def cut_mBs(df):
    n0 = len(df)
    if (n0 == 0):
        return df, 0
    else:
        df = df[(df['mBs'] > 5.3663 - 0.5) & (df['mBs'] < 5.3663 + 0.5)]
        n1 = len(df)
        return df, n1/n0


def cut_mPhi(df):
    n0 = len(df)
    if (n0 == 0):
        return df, 0
    else:
        df = df[(df['mPhi'] > 1) & (df['mPhi'] < 1.04)]
        n1 = len(df)
        return df, n1/n0


def cut_mDimu(df):
    n0 = len(df)
    if (n0 == 0):
        return df, 0
    else:
        df = df[(df['mDimu'] < 1) | (df['mDimu'] > 1.04)]
        df = df[(df['mDimu'] < 3.00) | (df['mDimu'] > 3.2)]
        n1 = len(df)
        return df, n1/n0
    


fn_dt = {
        'chi2': cut_chi2,
        'mDimu': cut_mDimu,
        'mPhi': cut_mPhi,
        'mBs': cut_mBs
        }

cut_eff_lt = []
df_cut_dt = {}
for k, df in df_dt.items():
    cut_eff_i = []

    for ck, cfn in fn_dt.items():
        df, eff = fn_dt[ck](df)
        cut_eff_i.append(eff)

    df_cut_dt[k] = df
    cut_eff_lt.append(cut_eff_i)


df_cut = pd.DataFrame(cut_eff_lt, columns=fn_dt.keys(), index=file_name_dt.keys())

print()
print(df_cut)
print()



target_folder_name = '../data/preselect/'
"""
target_file_name_dt = {"sig.":      target_folder_name + "/Bs_phimumu_reco.csv",
                        "comb.":    target_folder_name + "/Bs_comb_reco.csv",
                        "Jpsi":     target_folder_name + "/Bs_phiJpsi_reco.csv",
                        "psi2S":    target_folder_name + "/Bs_phipsi2S_reco.csv",
                        "phi1020":  target_folder_name + "/Bs_phiphi_reco.csv"}
"""

for k, df in df_cut_dt.items():
    #df.to_csv(target_file_name_dt[k], index=False)
    outName = file_name_dt[k].replace(folder_name, target_folder_name)
    outName = outName.replace("reco", "preselect")
    print("saving:", outName) 
    outFile = uproot.recreate(outName)
    outFile['t'] = df





