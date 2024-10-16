# plot the money plots

import uproot
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pyplot
import sys
import glob

# ../data/reco
FOLDER_NAME = sys.argv[1]

OUT_PATH = "../outputs/money_plots/"

CUTTED = 'preselect' in FOLDER_NAME

print("\n"*10)
print("plotting the cutted dataset?", CUTTED)

if not CUTTED: 
    """ Before Cutting """
    # just after the reconstruction 
    FILE_NAME_DT = {"sig.":     FOLDER_NAME + "ee2Z2Bs2PhiMuMu_reco.root",
                    "bb":       FOLDER_NAME + "ee2Z2b_comb_cutted_reco.root",
                    "cc":       FOLDER_NAME + "ee2Z2c_comb_cutted_reco.root",
                    }

    FILE_DT = {k: uproot.open(i) for k, i in FILE_NAME_DT.items()}
    DF_DT = {k: pd.DataFrame(np.array(file['t']['features'].array())) for k, file in FILE_DT.items()}

    # Yields (N_FS)
    # TODO: might need to update the number, after the theory prediciton
    YIELD_DT = {"sig.":     1.16e5,
                "bb":       4.34e9,
                "cc":       2.82e8,
                }

for k, df in DF_DT.items():
    print(k, len(df), f'{YIELD_DT[k]/len(df):.3e}')
print()


# {{{ config
# legnend
LEGEND_DT = {"sig.":     r'$B_s^0\to \phi \mu^+ \mu^-$',
             "bb":       r'$Z\to b\bar{b}$',
             "cc":       r'$Z\to c\bar{c}$',
             }


# unit of the y-axis
UNIT_DT = { 'mPhi':     'GeV',
            'q2Dimu':   r'GeV$^2$',
            'PKp':      'GeV',
            'PKm':      'GeV',
            'PK':       'GeV',
            'Pmup':     'GeV',
            'Pmum':     'GeV',
            'Pmu':      'GeV',
            'Chi2':     '',
            }


# the x-tick
TICK_DT = { 'mPhi':     r'$m_{\phi}$',
            'q2Dimu':   r'$q^2$',
            'PKp':      r'$|\vec{p}_{K^+}|$',
            'PKm':      r'$|\vec{p}_{K^-}|$',
            'PK':       r'$|\vec{p}_{K^\pm}|$',
            'Pmup':     r'$|\vec{p}_{\mu^+}|$',
            'Pmum':     r'$|\vec{p}_{\mu^-}|$',
            'Pmu':      r'$|\vec{p}_{\mu^\pm}|$',
            'Chi2':     r'$\chi^2$',
            }


# the x-range
RANGE_DT = {
            'Chi2':     [0, 75],
            'PKp':      [0, 0.67*50],
            'PKm':      [0, 0.67*50],
            'PK':       [0, 0.67*50],
            'Pmup':     [0, 0.67*50],
            'Pmum':     [0, 0.67*50],
            'q2Dimu':   [0, 25],
            'mPhi':     [1-0.005*3, 1.04+0.005*39],
            }


# color of the hist
COLORS_DT = {
            'sig.': 'C0',
            'bb':   'C1',
            'cc':   'C2'
            }
# }}}



# {{{ helper functions, and do some calculation in the df
def cal_q2(df):
    df['q2Dimu'] = df['mDimu']**2 
    return df

def cal_error(df):
    df['dx'] = df['DV_X'] - df['DV_X_truth']
    df['dy'] = df['DV_Y'] - df['DV_Y_truth']
    df['dz'] = df['DV_Z'] - df['DV_Z_truth']
    df['dl'] = (df['dx']**2 + df['dy']**2 + df['dz']**2)**0.5 
    return df

# calculate all q2 for signal & bkg
for k, df in DF_DT.items():
    df = cal_q2(df)
    DF_DT[k] = df

# calculate the vertex error for only signal
DF_DT['sig.'] = cal_error(DF_DT['sig.'])

# }}}



# plot the figures
def save_fig(DF_DT, label, range_dt, AU=False):
    # {{{
    out_name = OUT_PATH + label 
    print("Plotting", label)

    isP = 0 # if this is plotting momentum: If so, then we need to consider the +/- type
    if label == 'PK':
        label = 'PKp'
        isP = 1
    if label == 'Pmu':
        label = 'Pmup'
        isP = 1
        
    # book the figure and plot
    f = plt.figure(figsize=(8,6))
    ax1 = plt.subplot(111)

    # range of the plot
    range_ = find_plot_range(DF_DT, label) if range_dt is None else RANGE_DT[label]
    bins_ = 50
    log_sc = True   

    for k, df in DF_DT.items():
        scale_f = 1 if AU else YIELD_DT[k]  # normalization scale
        _, b, _ = ax1.hist(df[label], weights=[scale_f/len(df)]*len(df),
                  bins=bins_, range=range_, lw=2, histtype='step', alpha=0.5, label=LEGEND_DT[k],
                  log=log_sc)
        if isP: # for plotting momenta (the -ve is in "--" line style)
            _, b, _ = ax1.hist(df[label.replace('p', 'm')], weights=[scale_f/len(df)]*len(df),
                      bins=bins_, range=range_, lw=2, histtype='step', alpha=0.5,
                      log=log_sc, linestyle='--', color=COLORS_DT[k])

    ax1.yaxis.set_major_locator(matplotlib.ticker.LogLocator(numticks=999))
    ax1.yaxis.set_minor_locator(matplotlib.ticker.LogLocator(numticks=999, subs="auto"))

    if AU: 
        ax1.set_ylabel('A.U.', fontsize=22)
    else: 
        bins_size = (range_[1] - range_[0]) / bins_
        ax1.set_ylabel(f'Events / {bins_size:.2g} {UNIT_DT[label]}', fontsize=22)

    # add the grey region indicating the cut
    if label == 'Chi2':
        plt.axvspan(15, range_[1], alpha=0.15, color='grey')
    elif label == 'Pmup' or label == 'Pmum' or label == 'PKp' or label == 'PKm':
        plt.axvspan(range_[0], 2, alpha=0.15, color='grey')
    elif label == 'q2Dimu':
        plt.axvspan(1, 1.08, alpha=0.15, color='grey')
        plt.axvspan(9, 10.2, alpha=0.15, color='grey')
        plt.axvspan(13, 14.4, alpha=0.15, color='grey')
    elif label == 'mPhi':
        plt.axvspan(range_[0], 1, alpha=0.15, color='grey')
        plt.axvspan(1.04, range_[1], alpha=0.15, color='grey')

    # manually adjust the ylim in case the legend is overlapped
    if label == 'Pmup' or label == 'Pmum':
        plt.ylim(1e2, 1e10)

    x_label = f'{TICK_DT[label]}' if UNIT_DT[label] == '' else f'{TICK_DT[label]} [{UNIT_DT[label]}]'
    if isP: x_label = f"{TICK_DT[label.replace('p', '')]}" if UNIT_DT[label] == '' else f"{TICK_DT[label.replace('p', '')]} [{UNIT_DT[label]}]"
    ax1.set_xlabel(x_label, fontsize=22)
    ax1.tick_params(axis='both', which='major', labelsize=17.5)
    ax1.text(1, 1.03, r'FCC-$ee$ Simulation (IDEA)', horizontalalignment='right',
            verticalalignment='center', transform=ax1.transAxes, fontsize=17.5);
    plt.margins(x=0)
    plt.legend(fontsize=17.5)
    plt.tight_layout()

    # save
    if (CUTTED):
        out_name = out_name + "_cut"
    if (range_dt is None): 
        out_name = out_name + "_full"
    f.savefig(out_name + ".jpg")
    # f.savefig(out_name + ".pdf")

    # }}}



# plot the vertex resolution
def save_vertex_reso_fig(DF_DT, range_=[-0.1*1000, 0.1*1000], AU=False, chi2_cut=None):
    # {{{
    print("range:", range_)
    k = 'sig.'

    # book the figure and plot
    f = plt.figure(figsize=(8,6))
    ax1 = plt.subplot(111)

    ############################
    #        reslutionn        #
    ############################
    df = DF_DT[k]
    n0 = len(df)
    if chi2_cut is not None:
        df = df[df.Chi2 < chi2_cut]

    bins_ = 50
    # text_ = None
    log_sc = True
    scale_f = 1 if AU else YIELD_DT[k]
    label_lt = ['dx', 'dy', 'dz']
    colors_lt = ['C3', 'C8', 'C9']
    for i, label in enumerate(label_lt):
        sigm = df[(df[label] > range_[0]/1000) & (df[label] < range_[1]/1000)][label].std() * 1000
        text_ = rf'$\sigma_{label[-1]}=${sigm:.2g} $\mu$m'

        values = df.loc[:, label]*1000
        _, b, _ = ax1.hist(values, weights=[scale_f/len(df)]*len(df),
                  bins=bins_, range=range_, lw=2, histtype='step', alpha=0.5,
                  log=log_sc, label=text_, color=colors_lt[i])

    ax1.yaxis.set_major_locator(matplotlib.ticker.LogLocator(numticks=999))
    ax1.yaxis.set_minor_locator(matplotlib.ticker.LogLocator(numticks=999, subs="auto"))

    bins_size = (range_[1] - range_[0]) / bins_
    ax1.set_ylabel(f'Events / {bins_size:.2g} ' + r'$\mu$m', fontsize=22)
    x_label = r'$\vec{s}_i^{\rm reco} - \vec{s}_i^{\rm truth}$ [$\mu$m]'
    ax1.set_xlabel(x_label, fontsize=22)
    ax1.tick_params(axis='both', which='major', labelsize=17.5)
    ax1.text(1, 1.03, r'FCC-$ee$ Simulation (IDEA)', horizontalalignment='right',
            verticalalignment='center', transform=ax1.transAxes, fontsize=17.5);
    plt.margins(x=0)
    plt.legend(fontsize=17.5)
    plt.tight_layout()
    out_name = OUT_PATH + 'vertexReso'

    # save
    if chi2_cut is not None:
        out_name = out_name + "_chi2_" + str(chi2_cut)
    f.savefig(out_name + ".jpg")
    # f.savefig(out_name + ".pdf")

    # }}}




# {{{ Plotting 
save_vertex_reso_fig(DF_DT)
save_vertex_reso_fig(DF_DT, chi2_cut=15)

save_fig(DF_DT, 'Chi2', RANGE_DT)
save_fig(DF_DT, 'PKp', RANGE_DT)
save_fig(DF_DT, 'PKm', RANGE_DT)
save_fig(DF_DT, 'Pmup', RANGE_DT)
save_fig(DF_DT, 'Pmum', RANGE_DT)
save_fig(DF_DT, 'q2Dimu', RANGE_DT)
save_fig(DF_DT, 'mPhi', RANGE_DT)
save_fig(DF_DT, 'PK', RANGE_DT)
save_fig(DF_DT, 'Pmu', RANGE_DT)
# }}} 
