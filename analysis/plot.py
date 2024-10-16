import uproot
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#import seaborn as sns
import sys
import glob

FOLDER_NAME = sys.argv[1]

OUT_PATH = "../outputs/"

CUTTED = 'preselect' in FOLDER_NAME
print("\n"*10)
print(CUTTED)


if not CUTTED: 
    """ Before Cutting """
    FILE_NAME_DT = {"sig.":     FOLDER_NAME + "ee2Z2Bs2PhiMuMu_reco.root",
                    "bb":       FOLDER_NAME + "ee2Z2b_comb_cutted_reco.root",
                    "cc":       FOLDER_NAME + "ee2Z2c_comb_cutted_reco.root",
                    }

    FILE_DT = {k: uproot.open(i) for k, i in FILE_NAME_DT.items()}
    DF_DT = {k: pd.DataFrame(np.array(file['t']['features'].array())) for k, file in FILE_DT.items()}

    YIELD_DT = {"sig.":     2.76e5,
                "bb":       4.34e9,
                "cc":       2.82e8,
                }

'''
else: 
    """ After Cutting """
    FILE_NAME_DT = {"sig.":     FOLDER_NAME + "ee2Z2Bs2PhiMuMu_preselect.root",
                    "bb":       FOLDER_NAME + "ee2Z2b_comb_cutted_preselect.root",
                    "cc":       FOLDER_NAME + "ee2Z2c_comb_cutted_preselect.root",
                    }

    FILE_DT = {k: uproot.open(i) for k, i in FILE_NAME_DT.items()}
    DF_DT = {k: pd.DataFrame(np.array(FILE_DT[k]['t'].arrays())) for k, file in FILE_DT.items()}

    YIELD_DT = {"sig.":     2.72e5,
                "bb":       2.15e8,
                "cc":       8.95e6
                }
    #YIELD_DT = {"sig.":     1.46e5,
    #            "bb":       5.34e5,
    #            }
'''


for k, df in DF_DT.items():
    print(k)
    print(k, len(df), f'{YIELD_DT[k]/len(df):.3e}')
print()

LEGEND_DT = {"sig.":     r'$B_s\to \phi \mu \mu$',
             "bb":       r'$Z\to b\bar{b}$',
             "cc":       r'$Z\to c\bar{c}$',
             "Jpsi":     r'$B_s\to \phi J/\psi$',
             "psi2S":    r'$B_s\to \phi \psi $',
             "phi1020":  r'$B_s\to \phi \phi$'
             }

UNIT_DT = { 'mPhi':     'GeV',
            'mBs':      'GeV',
            'EBs':      'GeV',
            'mBs_reso': 'GeV',
            'mDimu':    'GeV',
            'q2Dimu':   r'GeV$^2$',
            'PKp':      'GeV',
            'PKm':      'GeV',
            'PK':       'GeV',
            'Pmup':     'GeV',
            'Pmum':     'GeV',
            'Pmu':      'GeV',
            'D0Kp':     'mm',
            'D0Km':     'mm',
            'D0mup':    'mm',
            'D0mum':    'mm',
            'DZKp':     'mm',
            'DZKm':     'mm',
            'DZmup':    'mm',
            'DZmum':    'mm',
            'D0min':    'mm',
            'DZmin':    'mm',
            'Chi2':     '',
            'Chi2_KK':  '',
            'cosTheta_dimu': '',
            'cosTheta_dimu_neg': '',
            'cosTheta_Bs': '',
            'DV':       'mm',
            'DV_truth': 'mm',
            'Delta_t':  's',
            'dx':       'mm',
            'dy':       'mm',
            'dz':       'mm',
            'dl':       'mm'
            }

TICK_DT = { 'mPhi':     r'$m_{KK}$',
            'mBs':      r'$m_{KK\mu\mu}$',
            'EBs':      r'$E_{KK\mu\mu}$',
            'mBs_reso': r'$m_{KK\mu\mu}$',
            'mDimu':    r'$m_{\mu\mu}$',
            'q2Dimu':   r'$q^2_{\mu\mu}$',
            'PKp':      r'$|\vec{p}_{K^+}|$',
            'PKm':      r'$|\vec{p}_{K^-}|$',
            'PK':       r'$|\vec{p}_{K^\pm}|$',
            'Pmup':     r'$|\vec{p}_{\mu^+}|$',
            'Pmum':     r'$|\vec{p}_{\mu^-}|$',
            'Pmu':      r'$|\vec{p}_{\mu^\pm}|$',
            'D0Kp':     r'$d_0(K^+)$',
            'D0Km':     r'$d_0(K^-)$',
            'D0mup':    r'$d_0(\mu^+)$',
            'D0mum':    r'$d_0(\mu^-)$',
            'DZKp':     r'$z_0(K^+)$',
            'DZKm':     r'$z_0(K^-)$',
            'DZmup':    r'$z_0(\mu^+)$',
            'DZmum':    r'$z_0(\mu^-)$',
            'D0min':    r'$d_0^{\rm min}$',
            'DZmin':    r'$z_0^{\rm min}$',
            'Chi2':     r'$\chi_{KK\mu\mu}^2$',
            'Chi2_KK':  r'$\chi_{KK}^2$',
            'cosTheta_dimu': r'$\cos{\theta}$',
            'cosTheta_dimu_neg': r'$\cos{\theta}$',
            'cosTheta_Bs': r'$\cos{\theta}$',
            'DV':       r'$|\vec{x}^{\rm reco}_{\rm DV}|$',
            'DV_truth': r'$|\vec{x}^{\rm truth}_{\rm DV}-\vec{x}^{\rm truth}_{\rm PV}|$',
            'Delta_t':  r'$\tau$',
            'dx':       r'$x^{\rm reco}_{\rm DV} - x^{\rm truth}_{\rm DV}$',
            'dy':       r'$y^{\rm reco}_{\rm DV} - y^{\rm truth}_{\rm DV}$',
            'dz':       r'$z^{\rm reco}_{\rm DV} - z^{\rm truth}_{\rm DV}$',
            'dl':       r'$|\vec{x}^{\rm reco}_{\rm DV} - \vec{x}^{\rm truth}_{\rm DV}|$'
            }

RANGE_DT = {'mPhi':         [0.98, 1.15],
            'mBs':          [4, 6],
            'EBs':          [0, 50],
            'mBs_reso':     [5.366-0.025, 5.366+0.025],
            'mDimu':        [0, 5],
            'q2Dimu':       [0, 25],
            'PKp':          [0, 40],
            'PKm':          [0, 40],
            'PK':           [0, 40],
            'Pmup':         [0, 40],
            'Pmum':         [0, 40],
            'Pmu':          [0, 40],
            'D0Kp':         [-7, 7],
            'D0Km':         [-7, 7],
            'D0mup':        [-16, 16],
            'D0mum':        [-16, 16],
            'DZKp':         [-20, 20],
            'DZKm':         [-20, 20],            
            'DZmup':        [-40, 40],
            'DZmum':        [-40, 40],
            'D0min':        [0, 16],
            'DZmin':        [0, 40],
            'Chi2':         [0, 35],
            'Chi2_KK':      [0, 35],
            'cosTheta_dimu':     [-1, 1],
            'cosTheta_dimu_neg':     [-1, 0],
            'cosTheta_Bs':  [-1, 1],
            'DV':           [0, 30],
            'DV_truth':     [0, 30],
            'Delta_t':      [0, 1e-11],
            'dx':           [-0.1, 0.1],
            'dy':           [-0.1, 0.1],
            'dz':           [-0.1, 0.1],
            'dl':           [0, 0.1]
            }






def cal_minD(df):
    df['D0min'] = df[['D0Kp', 'D0Km', 'D0mup', 'D0mum']].abs().min(axis=1)
    df['DZmin'] = df[['DZKp', 'DZKm', 'DZmup', 'DZmum']].abs().min(axis=1)
    return df

def cal_q2(df):
    df['q2Dimu'] = df['mDimu']**2 
    return df

def cal_error(df):
    df['dx'] = df['DV_X'] - df['DV_X_truth']
    df['dy'] = df['DV_Y'] - df['DV_Y_truth']
    df['dz'] = df['DV_Z'] - df['DV_Z_truth']
    df['dl'] = (df['dx']**2 + df['dy']**2 + df['dz']**2)**0.5 
    return df

for k, df in DF_DT.items():
    df = cal_q2(df)
    df = cal_minD(df)
    df['mBs_reso'] = df['mBs'].copy()
    df['cosTheta_dimu_neg'] = df['cosTheta_dimu'].copy()

    DF_DT[k] = df

DF_DT['sig.'] = cal_error(DF_DT['sig.'])
















def find_plot_range(DF_DT, label):
    min_ = 99999
    max_ = -99999

    for k, df in DF_DT.items():
        if (min_ > df[label].min()):
            min_ = df[label].min()
        if (max_ < df[label].max()):
            max_ = df[label].max()
    if (min_ > 0):
        min_range = min_ * 0.9
    else:
        min_range = min_ * 1.1
    if (max_ > 0):
        max_range = max_ * 1.1
    else:
        max_range = max_ * 0.9
    return [min_range, max_range]




def save_fig(DF_DT, label, range_dt, AU=False):
    # print out the status 
    if range_dt is not None: 
        print("Plotting", label, "\trange_dt:", range_dt[label])
    else: 
        print("Plotting", label)

    # book the figure and plot
    f = plt.figure(figsize=(8,6))
    ax1 = plt.subplot(111)

    # range of the plot
    range_ = find_plot_range(DF_DT, label) if range_dt is None else RANGE_DT[label]
    bins_ = 50
 
    ############################
    #        reslutionn        #
    ############################
    if label=='dx' or label=='dy' or label=='dz':
        df = DF_DT['sig.']
        sigm = df[(df[label] > RANGE_DT[label][0]) & (df[label] < RANGE_DT[label][1])][label].std()
        text_ = r'$\sigma=${:.4f}'.format(sigm)
        bins_ = 50
    elif label=='DV' or label=='dl':
        df = DF_DT['sig.']
        mea = df[(df[label] > RANGE_DT[label][0]) & (df[label] < RANGE_DT[label][1])][label].mean()
        text_ = r'sig. mean$=${:.4f}'.format(mea)
        bins_ = 50
    elif label=='mBs_reso':
        df = DF_DT['sig.']
        sigm = df[(df[label] > RANGE_DT[label][0]) & (df[label] < RANGE_DT[label][1])][label].std()
        text_ = r'$\sigma=${:.4f}'.format(sigm)
        bins_ = 50
    else: 
        text_ = None

    ############################
    #     truth level info     #
    ############################
    if (label =='DV_truth' or label=='dx' or label=='dy' or label=='dz' or label=='dl' or label=='mBs_reso'):
        DF_DT = {'sig.': DF_DT['sig.']}

    # change the y-axis to log scale for these variables
    #log_label_lt = ['mBs', 'EBs', 'mPhi', 'mDimu', 'q2Dimu', 'DV', 'Chi2', 'Chi2_KK', 'DV_truth', 'DV', 'dx', 'dy', 'dz', 'dl', 'cosTheta_dimu', 'cosTheta_dimu_neg', 
    #                 'PKp', 'PKm', 'Pmup', 'Pmum', 'D0Kp', 'cosTheta_Bs']
    for k, df in DF_DT.items():
        # log_sc = True if label in log_label_lt else False
        log_sc = True
        scale_f = 1 if AU else YIELD_DT[k]
        _, b, _ = ax1.hist(df[label], weights=[scale_f/len(df)]*len(df),
                  bins=bins_, range=range_, lw=3, histtype='step', alpha=0.5, label=LEGEND_DT[k],
                  log=log_sc)

    bins_size = (range_[1] - range_[0]) / bins_
    print(bins_size)

    ax1.text(0.8, 0.2, text_, fontsize=20, horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes)
    if AU: 
        ax1.set_ylabel('A.U.', fontsize=22)
        #ax1.set_yticks([])
    else: 
        ax1.set_ylabel(f'Events / {bins_size} {UNIT_DT[label]}', fontsize=22)

    if label == 'Chi2':
        plt.axvspan(15, range_[1], alpha=0.15, color='grey')
    elif label == 'PMup' or label == 'PMum' or label == 'PKp' or label == 'PKm':
        plt.axvspan(range_[0], 2, alpha=0.15, color='grey')
    elif label == 'q2Dimu':
        plt.axvspan(1, 1.0816, alpha=0.15, color='grey')
        plt.axvspan(9, 10.24, alpha=0.15, color='grey')
        plt.axvspan(12.96, 14.44, alpha=0.15, color='grey')
    elif label == 'mPhi':
        plt.axvspan(range_[0], 1, alpha=0.15, color='grey')
        plt.axvspan(1.04, range_[1], alpha=0.15, color='grey')
#    elif label == 'mBs':
#        plt.axvspan(range_[0], 4.8663, alpha=0.15, color='grey')
#        plt.axvspan(5.8663, range_[1], alpha=0.15, color='grey')




    x_label = f'{TICK_DT[label]}' if UNIT_DT[label] == '' else f'{TICK_DT[label]} [{UNIT_DT[label]}]'
    ax1.set_xlabel(x_label, fontsize=22)
    ax1.tick_params(axis='both', which='major', labelsize=17.5)
    plt.margins(x=0)
    plt.legend(fontsize=22)

    plt.tight_layout()
    out_name = OUT_PATH + label 
    if (CUTTED):
        out_name = out_name + "_cut"
    if (range_dt is None): 
        out_name = out_name + "_full"


    f.savefig(out_name + ".jpg")




def plot_osc_er(DF_DT, er=1, bins_=10):
    def get_er_ar(n_s, n_b, bins_):
        # er = np.sqrt(n_s)
        er_n = np.sqrt(n_s + n_b)
        er = np.divide(er_n, n_s, out=np.zeros_like(er_n), where=n_s!=0) * n_s
        half_bin_width = 0.5*(bins_[1:] + bins_[:-1])
        return er, half_bin_width


    
    #label = "DV"
    def get_Dt(df):
        df['gamma'] = df['EBs'] / (df['mBs'] * 3e8 * 3e8) 
        df['beta'] = df['PBs'] * 3e8 / df['EBs']
        df['Delta_t'] = df['DV'] * 0.001 / (df['beta'] * df['gamma'] * 3e8) / 3e8
        return df




    label = "Delta_t"

    f, axs = plt.subplots(2, 1, figsize=(8,6), sharex=True, gridspec_kw={'height_ratios': [3, 1]})
    f.subplots_adjust(hspace=0)
    ax1 = axs[0]
    #bins_ = 30
    df_sig = DF_DT['sig.']
    df_sig = get_Dt(df_sig)
    df_sig['w'] = np.array([1/len(df_sig)]* len(df_sig)) * YIELD_DT['sig.']
    #df_ii = df_sig[(df_sig.SS_truth == 1) & (df_sig.BBbar_truth == 1)]
    df_ii = df_sig[(df_sig.BBbar_truth == 1)]
    ns1, bins, _ = ax1.hist(df_ii[label], weights=df_ii['w'],
                            bins=bins_, range=RANGE_DT[label], lw=2, histtype='step', alpha=0.5,
                            label=r'$B\to K^{*} \mu \mu $', color='C0')
    #df_ii = df_sig[(df_sig.SS_truth == -1) & (df_sig.BBbar_truth == 1)]
    df_ii = df_sig[(df_sig.BBbar_truth == -1)]
    ns2, bins, _ = ax1.hist(df_ii[label], weights=df_ii['w'], 
                            bins=bins_, range=RANGE_DT[label], lw=2, histtype='step', alpha=0.5,
                            label=r'$B\to \bar{B}\to \bar{K}^{*} \mu \mu $', color='C1')
    """
    df_ii = df_sig[(df_sig.SS_truth == 1) & (df_sig.BBbar_truth == -1)]
    ns3, bins, _ = ax1.hist(df_ii[label], weights=df_ii['w'], 
                            bins=bins_, range=RANGE_DT[label], lw=2, histtype='step', alpha=0.5,
                            label=r'$\bar{B}\to \bar{K}^{*} \mu \mu $', color='C2')
    df_ii = df_sig[(df_sig.SS_truth == -1) & (df_sig.BBbar_truth == -1)]
    ns4, bins, _ = ax1.hist(df_ii[label], weights=df_ii['w'],
                            bins=bins_, range=RANGE_DT[label], lw=2, histtype='step', alpha=0.5,
                            label=r'$\bar{B}\to B \to K^{*} \mu \mu $', color='C3')
    """

    if er == 0:
        if CUTTED: 
            rescale = {"bb":     0,
                       "Jpsi":      0, 
                       "psi2S":     0,
                       "phi1020":   0}
            """
            rescale = {"bb":     1,
                       "Jpsi":      1, 
                       "psi2S":     10,
                       "phi1020":   10000}
            """
        else:
            rescale = {"bb":    0,
                       "cc":    0}
            """
            rescale = {"bb":     1e-5,
                       "Jpsi":      1e-2, 
                       "psi2S":     1e-1,
                       "phi1020":   2e2}
            """

        color_dt = {"bb":   'C4',
                    "cc":   'C5'}
        hist_lb = {"bb":    r'$Z\to b \bar{b}'+r'$\times$'+f"{rescale['bb']:.1E}",
                   "cc":    r'$Z\to c \bar{c}'+r'$\times$'+f"{rescale['cc']:.1E}"}
        for k, df in DF_DT.items():
            if (k == 'sig.'): 
                continue
            df_bkg = DF_DT[k]
            df_bkg = get_Dt(df_bkg)
            df_bkg['w'] = np.array([1/len(df_bkg)]* len(df_bkg)) * YIELD_DT[k] * rescale[k]
            _, _, _ = ax1.hist(df_bkg[label], weights=df_bkg['w'],
                                    bins=bins_, range=RANGE_DT[label], lw=2, histtype='step', alpha=0.5,
                                    label=hist_lb[k], color=color_dt[k])



    ax2 = axs[1]
    a = ns1-ns2
    b = ns1+ns2
    d = np.divide(a, b, out=np.zeros_like(a), where=b!=0)
    bins_mid = (bins[:-1] + bins[1:]) / 2
    ax2.plot(bins_mid, d, color='blue')
    if er == 1:
        lt_nKstar = []
        lt_nKstarbar = []
        for k, df in DF_DT.items():
            if (k == 'sig.'): 
                continue
            df_bkg = DF_DT[k]
            df_bkg = get_Dt(df_bkg)
            df_bkg['w'] = np.array([1/len(df_bkg)]* len(df_bkg)) * YIELD_DT[k]
            df_Kstar = df_bkg[df_bkg.Kstarbar == 1]
            df_Kstarbar = df_bkg[df_bkg.Kstarbar == -1]
            nKstar, bins = np.histogram(df_Kstar[label], weights=df_Kstar['w'], bins=bins_, range=RANGE_DT[label])
            lt_nKstar.append(nKstar)
            nKstarbar, bins = np.histogram(df_Kstarbar[label], weights=df_Kstarbar['w'], bins=bins_, range=RANGE_DT[label])
            lt_nKstarbar.append(nKstarbar)
        nKstar = np.sum(lt_nKstar)
        nKstarbar = np.sum(lt_nKstarbar)


        er1, hbins1 = get_er_ar(ns1, nKstar, bins)
        er2, hbins2 = get_er_ar(ns2, nKstarbar, bins)
        er3, hbins3 = get_er_ar(ns3, nKstarbar, bins)
        er4, hbins4 = get_er_ar(ns4, nKstar, bins)
        #er1, hbins1 = get_er_ar(ns1, ns4 + nKstar, bins)
        #er2, hbins2 = get_er_ar(ns2, ns2 + nKstarbar, bins)
        #er3, hbins3 = get_er_ar(ns3, ns3 + nKstarbar, bins)
        #er4, hbins4 = get_er_ar(ns4, ns1 + nKstar, bins)
        ax1.errorbar(hbins1, ns1, yerr=er1, fmt="none", color='C0', capsize=4)
        ax1.errorbar(hbins2, ns2, yerr=er2, fmt="none", color='C1', capsize=4)
        ax1.errorbar(hbins3, ns3, yerr=er3, fmt="none", color='C2', capsize=4)
        ax1.errorbar(hbins4, ns4, yerr=er4, fmt="none", color='C3', capsize=4)
        a_er = (er1**2 + er2**2 + er3**2 + er4**2)**0.5
        b_er = a_er
        a_d = np.divide(a_er, a, out=np.zeros_like(a_er), where=a!=0)
        b_d = np.divide(b_er, b, out=np.zeros_like(b_er), where=b!=0)
        d_er = np.abs(d) * (a_d**2 + b_d**2)**0.5
        ax2.errorbar(bins_mid, d, yerr=d_er, fmt="none", color='blue', capsize=4)

    """

    ax1.set_ylim(0, ax1.get_ylim()[1])
    ax1.set_ylabel('# Events', fontsize=22)
    # ax1.set_ylabel('A.U.', fontsize=22)
    # ax1.set_yticks([])
    ax1.margins(x=0)
    ax1.tick_params(axis='y', which='major', labelsize=17.5)
    ax1.legend(fontsize=17)
    ax1.yaxis.get_major_formatter().set_useOffset(True)
    """

    ax2.axhline(0, color='black')
    ax2.set_ylabel('Asymmetry', fontsize=22)
    ax2.set_xlabel(TICK_DT[label], fontsize=22)
    ax2.tick_params(axis='both', which='major', labelsize=17.5)
    ax2.margins(x=0)
    #ax2.set_ylim(-1.2, 1.2)
    ax2.set_ylim(-0.2, 0.2)
    plt.tight_layout()
    if er==1:
        f.savefig(OUT_PATH + label + "_Mix_err.jpg")
    elif er==0:
        f.savefig(OUT_PATH + label + "_Mix.jpg")

    if er == 1:
        ns1_f = np.array(np.sum(ns1))
        ns2_f = np.array(np.sum(ns2))
        ns3_f = np.array(np.sum(ns3))
        ns4_f = np.array(np.sum(ns4))
        nKstar_f = np.array(np.sum(nKstar))
        nKstarbar_f = np.array(np.sum(nKstarbar))
        a_f = ns1_f+ns3_f-ns2_f-ns4_f
        b_f = ns1_f+ns2_f+ns3_f+ns4_f
        d_f = np.divide(a_f, b_f, out=np.zeros_like(a_f), where=b_f!=0)
        er1_f, _ = get_er_ar(ns1_f, nKstar_f, bins[:1])
        er2_f, _ = get_er_ar(ns2_f, nKstarbar_f, bins[:1])
        er3_f, _ = get_er_ar(ns3_f, nKstarbar_f, bins[:1])
        er4_f, _ = get_er_ar(ns4_f, nKstar_f, bins[:1])
        a_f_er = (er1_f**2 + er2_f**2 + er3_f**2 + er4_f**2)**0.5
        b_f_er = a_f_er
        a_f_d = np.divide(a_f_er, a_f, out=np.zeros_like(a_f_er), where=a_f!=0)
        b_f_d = np.divide(b_f_er, b_f, out=np.zeros_like(b_f_er), where=b_f!=0)
        d_f_er = np.abs(d_f) * (a_f_d**2 + b_f_d**2)**0.5
        print(d_f, d_f*d_f_er)
        print("rel. precision of A_CP:", d_f_er)







print()
print("weights")
#for k, v in YIELD_DT.items():
#    print(k, v/len(DF_DT[k]))





# Money plots

save_fig(DF_DT, 'Chi2', RANGE_DT)
save_fig(DF_DT, 'PKp', RANGE_DT)
save_fig(DF_DT, 'PKm', RANGE_DT)
save_fig(DF_DT, 'Pmup', RANGE_DT)
save_fig(DF_DT, 'Pmum', RANGE_DT)


"""
save_fig(DF_DT, 'mPhi', None)
save_fig(DF_DT, 'mPhi', RANGE_DT)
save_fig(DF_DT, 'Chi2', None)
save_fig(DF_DT, 'Chi2_KK', None)
save_fig(DF_DT, 'Chi2_KK', RANGE_DT)
save_fig(DF_DT, 'cosTheta_dimu', RANGE_DT)
save_fig(DF_DT, 'cosTheta_dimu_neg', RANGE_DT)
save_fig(DF_DT, 'cosTheta_Bs', None)
save_fig(DF_DT, 'cosTheta_Bs', RANGE_DT)
save_fig(DF_DT, 'DV', None)
save_fig(DF_DT, 'DV', RANGE_DT)
save_fig(DF_DT, 'EBs', None)
save_fig(DF_DT, 'EBs', RANGE_DT)
save_fig(DF_DT, 'mBs', None)
save_fig(DF_DT, 'mBs', RANGE_DT)
save_fig(DF_DT, 'mDimu', None)
save_fig(DF_DT, 'mDimu', RANGE_DT)
save_fig(DF_DT, 'q2Dimu', None)
save_fig(DF_DT, 'q2Dimu', RANGE_DT)
save_fig(DF_DT, 'DV_truth', None)
save_fig(DF_DT, 'DV_truth', RANGE_DT)


save_fig(DF_DT, 'dx', RANGE_DT)
save_fig(DF_DT, 'dy', RANGE_DT)
save_fig(DF_DT, 'dz', RANGE_DT)
save_fig(DF_DT, 'dl', RANGE_DT)
save_fig(DF_DT, 'mBs_reso', RANGE_DT)
save_fig(DF_DT, 'D0Kp', None)
save_fig(DF_DT, 'D0Kp', RANGE_DT)
save_fig(DF_DT, 'D0Km', None)
save_fig(DF_DT, 'D0Km', RANGE_DT)
save_fig(DF_DT, 'D0mup', None)
save_fig(DF_DT, 'D0mup', RANGE_DT)
save_fig(DF_DT, 'D0mum', None)
save_fig(DF_DT, 'D0mum', RANGE_DT)

save_fig(DF_DT, 'DZKp', None)
save_fig(DF_DT, 'DZKp', RANGE_DT)
save_fig(DF_DT, 'DZKm', None)
save_fig(DF_DT, 'DZKm', RANGE_DT)
save_fig(DF_DT, 'DZmup', None)
save_fig(DF_DT, 'DZmup', RANGE_DT)
save_fig(DF_DT, 'DZmum', None)
save_fig(DF_DT, 'DZmum', RANGE_DT)
save_fig(DF_DT, 'D0min', None)
save_fig(DF_DT, 'D0min', RANGE_DT)
save_fig(DF_DT, 'DZmin', None)
save_fig(DF_DT, 'DZmin', RANGE_DT)

plot_osc_er(DF_DT, 0, 40)
#plot_osc_er(DF_DT)
"""



