import os
import matplotlib as mpl
mpl.use('pdf')
from matplotlib import pyplot as plt
from matplotlib import rc
#.Allow for using TeX mode in matplotlib Figures
rc('font',**{'family':'sans-serif','sans-serif':['Computer Modern Roman']})
rc('text', usetex=True)
plt.rcParams['text.latex.preamble']=[r"\usepackage{lmodern}"]

ratio=5.0/7.0
fig_width_pt = 3*246.0  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inch
golden_mean = ratio if ratio != 0.0 else (np.sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_width*golden_mean      # height in inches
fig_size =  [fig_width,fig_height]

params = {'text.usetex' : True,
        'axes.labelsize': 24,
        'font.size': 24,
        'legend.fontsize': 10,
        'xtick.labelsize': 20,
        'ytick.labelsize': 20,
        'font.family' : 'lmodern',
        'text.latex.unicode': True,
        'axes.grid' : True,
        'text.usetex': True,
        'figure.figsize': fig_size}
plt.rcParams.update(params)

import uproot
import glob
import pandas as pd
from tqdm import tqdm
import numpy as np
from sklearn.model_selection import train_test_split, StratifiedKFold, KFold
from sklearn.metrics import classification_report, roc_curve, auc, roc_auc_score, precision_score
import xgboost as xgb
from sklearn import metrics
import numpy as np
from skopt import gp_minimize
from skopt.space import Real, Integer
from skopt.utils import use_named_args
from skopt.plots import plot_convergence, plot_evaluations, plot_objective
from scipy import interp
import h5py
import time
#import ROOT
import sys
from matplotlib.backends.backend_pdf import PdfPages



rseed=6
np.random.seed(rseed)

#vardir = {}
#vardir["tau_pt"] = {'bin':np.linspace(0, 15, 20),'xtitle':'tau pT (GeV)'}
#vardir["tau_vprob"] = {'bin':np.linspace(0, 1, 20), 'xtitle':'Tau vertex prob.'}
#vardir["tau_fls3d"] = {'bin':np.linspace(0, 10, 20), 'xtitle':'Tau flight length sig.'}
#vardir["tau_sumofdnn"] = {'bin':np.linspace(0, 3, 20), 'xtitle':'Tau sum of dnn'}
#vardir["tau_sumofdnn_others"] = {'bin':np.linspace(0, 3, 20), 'xtitle':'Tau sum of dnn'}
#
#vardir["b_mindoca"] = {'bin':np.linspace(0, 0.01, 20), 'xtitle':'B mindoca'}
#
#vardir["b_pt"] = {'bin':np.linspace(0, 45, 20), 'xtitle':'B pT (GeV)'}
#vardir["b_eta"] = {'tree':'tree',  'nbin':20, 'xmin':-3., 'xmax':3., 'xtitle':'B eta'}
#vardir["b_iso"] = {'bin':np.linspace(0, 200, 20), 'xtitle':'B iso'}
#vardir["b_nch"] = {'bin':np.linspace(0, 120, 120), 'xtitle':'Num of charged hadrons'}
#vardir["b_iso_ntracks"] = {'bin':np.linspace(0, 10, 10), 'xtitle':'B iso ntracks'}
#vardir["b_vprob"] = {'bin':np.linspace(0, 1, 20), 'xtitle':'B vertex prob.'}
#vardir["b_alpha"] = {'bin':np.linspace(0.95, 1, 20), 'xtitle':'B cos(alpha)'}
#vardir["b_pvips"] = {'bin':np.linspace(0, 30, 20), 'xtitle':'B pvip'}
#vardir["b_fls3d"] = {'bin':np.linspace(0, 10, 20), 'xtitle':'B flight length sig.'}
#vardir["b_lips"] = {'bin':np.linspace(-20, 20, 20), 'xtitle':'B lip'}
#vardir["dz_b_pv"] = {'bin':np.linspace(-1, 1, 20), 'xtitle':'Deltaz (B, bbPV) (cm)'}
#vardir["dz_jpsi_tau"] = {'bin':np.linspace(-1, 1, 20), 'xtitle':'Deltaz (J/psi, tau) (cm)'}
#
#vardir["ncand"] = {'bin':np.linspace(0, 10, 10), 'xtitle':'Nr. of tau cand.'}
#vardir["met"] = {'bin':np.linspace(0, 50, 20), 'xtitle':'puppi MET'}



def ensureDir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)


def get_df(root_file_name, branches):
    f = uproot.open(root_file_name)
    if len(f.allkeys()) == 0:
        return pd.DataFrame()
    df = pd.DataFrame(uproot.open(root_file_name)["tree"].arrays(branches, namedecode="utf-8"))
    return df

def get_label(name):
    if name == 0:
        return "background"
    else:
        return "signal"

def plot_roc_curve(df, score_column, tpr_threshold=0.0, ax=None, color=None, linestyle='-', label=None):
    print('Plotting ROC...')
    if ax is None:
        ax = plt.gca()
    if label is None:
        label = score_column
    fpr, tpr, thresholds = roc_curve(df["isSignal"], df[score_column], drop_intermediate=True)
    roc_auc = roc_auc_score(df["isSignal"], df[score_column])
    roc_pauc = roc_auc_score(df["isSignal"], df[score_column], max_fpr=1.0e-2)
    print("auc: {}, pauc: {}".format(roc_auc, roc_pauc))
    mask = tpr > tpr_threshold
    fpr, tpr = fpr[mask], tpr[mask]
    #ax.semilogx(fpr, tpr, label=label, color=color, linestyle=linestyle)
    ax.plot(fpr, tpr, label=label, color=color, linestyle=linestyle)

def get_signal_efficiency(df, score_column, working_point):
    df_sig = df.query("matchedToGenEle == 1")
    k = len(df_sig[df_sig[score_column] >= working_point])
    n = len(df_sig)
    return 1.*k/n if n != 0 else np.nan

def get_background_efficiency(df, score_column, working_point):
    df_bkg = df.query("matchedToGenEle == 0")
    k = len(df_bkg[df_bkg[score_column] >= working_point])
    n = len(df_bkg)
    return 1.*k/n if n != 0 else np.nan

#def get_signal_efficiency_unc(df, score_column, working_point, bUpper):
#    df_sig = df.query("matchedToGenEle == 1")
#    k = len(df_sig[df_sig[score_column] >= working_point])
#    n = len(df_sig)
#    teff = ROOT.TEfficiency()
#    return teff.Bayesian(n, k, 0.683, 1.0, 1.0, bUpper, True) if n != 0 else np.nan
#
#def get_background_efficiency_unc(df, score_column, working_point, bUpper):
#    df_bkg = df.query("matchedToGenEle == 0")
#    k = len(df_bkg[df_bkg[score_column] >= working_point])
#    n = len(df_bkg)
#    teff = ROOT.TEfficiency()
#    return teff.Bayesian(n, k, 0.683, 1.0, 1.0, bUpper, True) if n != 0 else np.nan

def fpreproc(dtrain, dtest, param):
    label = dtrain.get_label()
    ratio = float(np.sum(label == 0)) / np.sum(label == 1)
    param['scale_pos_weight'] = ratio
    return (dtrain, dtest, param)

def pauc(predt, dtrain):
    y = dtrain.get_label()
    return 'pauc', roc_auc_score(y, predt, max_fpr=1.0e-2)

#def prec(predt, dtrain):
#    y = dtrain.get_label()
#    return 'prec', precision_score(y, predt)

space  = [Integer(5, 8, name='max_depth'),
         Real(0.01, 0.2, name='eta'),
         Real(0.0, 10.0, name='gamma'),
         Integer(1.0, 10.0, name='min_child_weight'),
         Real(0.5, 1.0, name='subsample'),
         Real(0.1, 1.0, name='colsample_bytree'),
         Real(0.0, 10.0, name='alpha'),
         Real(0.0, 10.0, name='lambda'),
         ]

@use_named_args(space)
def objective(**X):
    global best_auc, best_auc_std, best_params
    print("New configuration: {}".format(X))
    params = X.copy()
    params['objective'] = 'binary:logitraw'
    #params['eval_metric'] = 'auc'
    params['nthread'] = 6
    params['silent'] = 1
    cv_result = xgb.cv(params, dmatrix_train, num_boost_round=n_boost_rounds, nfold=5, shuffle=True, stratified=True, maximize=True, early_stopping_rounds=75, fpreproc=fpreproc, feval=pauc)
#    cv_result = xgb.cv(params, dmatrix_train, num_boost_round=n_boost_rounds, nfold=5, shuffle=True, stratified=True, maximize=True, early_stopping_rounds=75, fpreproc=fpreproc, feval=prec)
    ave_auc = cv_result['test-pauc-mean'].iloc[-1]
    ave_auc_std = cv_result['test-pauc-std'].iloc[-1]
    print("Average pauc: {}+-{}".format(ave_auc, ave_auc_std))
    if ave_auc > best_auc:
      best_auc = ave_auc
      best_auc_std = ave_auc_std
      best_params = X.copy()
    print("Best pauc: {}+-{}, Best configuration: {}".format(best_auc, best_auc_std, best_params))
    return -ave_auc

def train(xgtrain, xgtest, hyper_params=None):
    watchlist = [(xgtrain, 'train'), (xgtest, 'eval')]
    params = hyper_params.copy()
    label = xgtrain.get_label()
    ratio = float(np.sum(label == 0)) / np.sum(label == 1)
    params['scale_pos_weight'] = ratio
    params['objective'] = 'binary:logitraw'
    #params['eval_metric'] = 'auc'
    params['nthread'] = 10
    params['silent'] = 1
    params['seed'] = rseed
    results = {}
    model = xgb.train(params, xgtrain, num_boost_round=n_boost_rounds, evals=watchlist, evals_result=results, maximize=True, early_stopping_rounds=75, verbose_eval=False, feval=pauc)
#    model = xgb.train(params, xgtrain, num_boost_round=n_boost_rounds, evals=watchlist, evals_result=results, maximize=True, early_stopping_rounds=75, verbose_eval=False, feval=prec)
    best_iteration = model.best_iteration + 1
    if best_iteration < n_boost_rounds:
        print("early stopping after {0} boosting rounds".format(best_iteration))
    return model, results

def train_cv(X_train_val, Y_train_val, X_test, Y_test, w_train_val, hyper_params=None):
    xgtrain = xgb.DMatrix(X_train_val, label=Y_train_val, weight=w_train_val)
    xgtest  = xgb.DMatrix(X_test , label=Y_test )
    model, results = train(xgtrain, xgtest, hyper_params=hyper_params)
    Y_predict = model.predict(xgtest, ntree_limit=model.best_ntree_limit)
    fpr, tpr, thresholds = roc_curve(Y_test, Y_predict, drop_intermediate=True)
    roc_auc = roc_auc_score(Y_test, Y_predict, max_fpr=1.0e-2)
    print("Best pauc: {}".format(roc_auc))
    return model, fpr, tpr, thresholds, roc_auc, results


def get_label(name):
    if name == 0:
      return "Bkg"
        #return "EB" #"B+ to K+ ee"
    if name == 1:
      return "Signal" #"B0 to K* ee"


#def plot_hist(df, column, bins=None, logscale=False, ax=None, title=None):
#    if ax is None:
#        ax = plt.gca()
#    #ax.hist(np.clip(df[column], bins[0], bins[-1]), bins=bins, histtype='step', normed=True)
#    for name, group in df.groupby("isSignal"):
#        print(name)
#        #ax.hist(np.clip(group[column], bins[0], bins[-1]), bins=bins, histtype='step', label=get_label(name), normed=True)
#        ax.hist(group[column], bins=bins, histtype='step', label=get_label(name), normed=True)
#        #ax.hist(group[column], bins=bins, histtype='step', label='{0}, Mean: {1:.2f}'.format(get_label(name), np.mean(group[column])), normed=True)
#    ax.set_ylabel("density")
#    ax.set_xlabel(vardir[column]['xtitle'])
#    ax.legend()
#    ax.set_title(title)
#    if logscale:
#        ax.set_yscale("log", nonposy='clip')



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="A simple ttree plotter")
    parser.add_argument("-s", "--signal", dest="signal", default="/work/ytakahas/work/BsTauTau/CMSSW_10_2_10/src/job/BcJpsiTau_large_BcJpsiTauNu_truth_2020-04-06-142954/Myroot_new.root", help="Signal file")
    parser.add_argument("-b", "--background", dest="background", default="/work/ytakahas/work/BsTauTau/CMSSW_10_2_10/src/job/BJpsiX_BcJpsiTauNu_mc_2020-04-06-112822/Myroot_new.root")
    parser.add_argument("-f", "--suffix", dest="suffix", default=None, help="Suffix of the output name")
    parser.add_argument("-o", "--optimization", dest="optimization", action='store_true', help="Perform Bayesian optimization")
    parser.add_argument("-d", "--dir", dest="dir", default="normal", help="output directory", required=True)
    args = parser.parse_args()

    pdir='Plots_' + args.dir
    mdir='model_' + args.dir

    ensureDir(pdir)
    ensureDir(mdir)

#    features = ['b_pt', 'b_eta', 'b_alpha', 'b_vprob', 'b_iso', 'b_iso_ntracks', 'b_lips', 'b_pvips', 'b_mindoca', 'dr_b_pv', 'dr_jpsi_tau', 'tau_fls3d', 'tau_vprob', 'tau_sumofdnn', 'ncand', 'estar']
#    features = ['b_pt', 'b_eta', 'b_alpha', 'b_vprob', 'b_iso', 'b_iso_ntracks', 'b_lips', 'b_pvips', 'b_mindoca', 'dr_b_pv', 'dr_jpsi_tau', 'tau_fls3d', 'tau_vprob', 'tau_sumofdnn', 'ncand', 'estar', 'nch_after_dnn', 'ptmiss']


# model_wmass
    features = None

#    if args.dir=='normal':

# First trial ...
#    features = ['b_pt', 'b_eta', 'b_alpha', 'b_vprob', 'b_iso_0p7', 'b_iso_ntracks_0p7', 'b_lips', 'b_pvips', 'b_mindoca', 'dr_b_pv', 'dr_jpsi_tau', 'tau_fls3d', 'tau_vprob', 'tau_sumofdnn_old', 'ncand', 'estar']

#    features = ['b_pt', 'b_eta', 'b_alpha', 'b_vprob', 'tau_iso_0p7', 'tau_iso_ntracks_0p7', 'b_lips', 'b_pvips', 'b_mindoca', 'dr_b_pv', 'dr_jpsi_tau', 'tau_fls3d', 'tau_vprob', 'tau_sumofdnn', 'tau_sumofdnn_1prong', 'tau_sumofdnn_otherB', 'tau_sumofdnn_pu', 'ncand', 'estar']

#    features = ['b_pt', 'b_eta', 'b_alpha', 'b_vprob', 'tau_iso_0p7', 'tau_iso_ntracks_0p7', 'b_lips', 'b_pvips', 'b_mindoca', 'dr_b_pv', 'tau_fls3d', 'tau_vprob', 'tau_fls3d_wjpsi', 'ncand', 'estar']
#    features = ['b_pt', 'b_eta', 'b_alpha', 'b_vprob', 'tau_iso_0p7',  'b_lips', 'b_pvips', 'b_mindoca', 'dr_b_pv', 'tau_fls3d', 'tau_vprob', 'tau_fls3d_wjpsi',  'tau_sumofdnn','tau_sumofdnn_1prong', 'tau_sumofdnn_otherB', 'tau_sumofdnn_pu', 'ncand', 'estar']
    features = ['b_pt', 'b_eta', 'b_alpha', 'b_vprob', 'tau_iso_0p7',  'b_lips', 'b_pvips', 'b_mindoca', 'dr_b_pv', 'tau_fls3d', 'tau_vprob', 'tau_fls3d_wjpsi',  'tau_sumofdnn', 'ncand', 'estar']
#    features = ['b_pt', 'b_eta', 'b_alpha', 'b_vprob', 'tau_iso_0p7',  'b_lips', 'b_pvips', 'b_mindoca', 'dr_b_pv', 'tau_fls3d', 'tau_vprob', 'tau_fls3d_wjpsi',  'tau_sumofdnn', 'ncand', 'estar']






    print('features = ', features)
    features = sorted(features)

    branches_bg = features + ['weight']
    branches_sig = features + ['weight', 'tau_isRight_3prong']
#    branches_bg = features
#    branches_sig = features + ['tau_isRight_3prong']
#    branches = features

    ddf = {}
    ddf['sig'] = get_df(args.signal, branches_sig)
    ddf['bkg'] = get_df(args.background, branches_bg)

    ddf['sig'].replace([np.inf, -np.inf], 0.0, inplace=True)
    ddf['bkg'].replace([np.inf, -np.inf], 0.0, inplace=True)

    # only allow 3prong
    print(ddf['sig'])
    ddf['sig'].query("tau_isRight_3prong==1", inplace=True)
    print(ddf['sig'])
    
    nSig = ddf['sig'].shape[0]
#    nBkg = 30000
    nBkg = 100000
#    nBkg = 300000
#    nBkg = 1000000
    #nSig = 10000
    #nBkg = 10000

    print('nSig=', nSig, 'nBkg=', nBkg)
    ddf['sig'] = ddf['sig'].sample(frac=1)[:nSig]
    ddf['bkg'] = ddf['bkg'].sample(frac=1)[:nBkg]

    # add isSignal variable
    ddf['sig']['isSignal'] = 1
    ddf['bkg']['isSignal'] = 0


    df = pd.concat([ddf['sig'],ddf['bkg']]).sort_index(axis=1).sample(frac=1).reset_index(drop=True)
#    df['weights'] = np.where(df['isSignal'], 1.0/df['BToKEE_fit_massErr'].replace(np.nan, 1.0), 1.0)
#    df['weights'] = np.where(df['isSignal']==0, df['weight'], 1.0)
#    df['weights'] = np.where(df['isSignal']==0, 1.0, 1.0)
    df['weights'] = df['weight']

    X = df[features]
    y = df['isSignal']
    W = df['weights']

    suffix = args.suffix
    n_boost_rounds = 800
    n_calls = 80
    n_random_starts = 40
    do_bo = args.optimization
    do_cv = True
    best_params = {'colsample_bytree': 0.8380017432637168, 'subsample': 0.7771020436861611, 'eta': 0.043554653675279234, 'alpha': 0.13978587730419964, 'max_depth': 3, 'gamma': 0.5966218064835417, 'lambda': 1.380893119219306}

    # split X and y up in train and test samples
    X_train, X_test, y_train, y_test = train_test_split(X, y, stratify=y, test_size=0.20, random_state=rseed)
#    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.20, shuffle=None, random_state=42)
#    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.20, random_state=42)

    # Get the number of positive and nevative training examples in this category
    n_pos = np.sum(y_train == 1)
    n_neg = np.sum(y_train == 0)

    print("training on {0} signal and {1} background".format(n_pos, n_neg))
   
    # get the indices correspondng to the testing and training samples
    idx_train = X_train.index
    idx_test = X_test.index

    w_train = W.loc[idx_train]

    dmatrix_train = xgb.DMatrix(X_train.copy(), label=np.copy(y_train), feature_names=[f.replace('_','-') for f in features], weight=np.copy(w_train))
    dmatrix_test  = xgb.DMatrix(X_test.copy(), label=np.copy(y_test), feature_names=[f.replace('_','-') for f in features])

    # Bayesian optimization
    if do_bo:
        begt = time.time()
        print("Begin Bayesian optimization")
        best_auc = 0.0
        best_auc_std = 0.0
        best_params = {}
        res_gp = gp_minimize(objective, space, n_calls=n_calls, n_random_starts=n_random_starts, verbose=True, random_state=36)
        print("Finish optimization in {}s".format(time.time()-begt))
        plt.figure()
        plot_convergence(res_gp)
        plt.savefig(pdir + '/training_resultis_bo_convergencePlot_xgb_{}.pdf'.format(suffix))
        #plt.figure()
        #plot_evaluations(res_gp)
        #plt.savefig('training_resultis_bo_evaluationsPlot_xgb_{}.pdf'.format(suffix))
        #plt.figure()
        #plot_objective(res_gp)
        #plt.savefig('training_resultis_bo_objectivePlot_xgb_{}.pdf'.format(suffix))

    # Get the cv plots with the best hyper-parameters
    if do_bo or do_cv:
        print("Get the cv plots with the best hyper-parameters")
        tprs = []
        aucs = []
        figs, axs = plt.subplots()
        #cv = KFold(n_splits=5, shuffle=True)
        cv = StratifiedKFold(n_splits=5, shuffle=True)
        mean_fpr = np.logspace(-6, 0, 100)

        iFold = 0
        for train_idx, test_idx in cv.split(X_train, y_train):
            X_train_cv = X_train.iloc[train_idx]
            X_test_cv = X_train.iloc[test_idx]
            Y_train_cv = y_train.iloc[train_idx]
            Y_test_cv = y_train.iloc[test_idx]
            w_train_cv = w_train.loc[X_train_cv.index]

            model, fpr, tpr, thresholds, roc_auc, results = train_cv(X_train_cv, Y_train_cv, X_test_cv, Y_test_cv, w_train_cv, hyper_params=best_params)
            epochs = len(results['train']['pauc'])
            x_axis = range(0, epochs)
            fig, ax = plt.subplots()
            ax.plot(x_axis, results['train']['pauc'], label='Train')
            ax.plot(x_axis, results['eval']['pauc'], label='Test')
            ax.legend()
            plt.ylabel('PAUC')
            plt.xlabel('Epoch')
            plt.title('Fold: {}'.format(iFold))
            fig.savefig(pdir + '/training_results_learning_curve_cv_fold_{}_{}.pdf'.format(suffix,iFold), bbox_inches='tight')

            tprs.append(interp(mean_fpr, fpr, tpr))
            tprs[-1][0] = 0.0
            aucs.append(roc_auc)
            axs.plot(fpr, tpr, lw=1, alpha=0.3, label='ROC fold %d (PAUC = %0.2f)' % (iFold, roc_auc))
            iFold += 1

        mean_tpr = np.mean(tprs, axis=0)
        mean_tpr[-1] = 1.0
        mean_auc = auc(mean_fpr, mean_tpr)
        std_auc = np.std(aucs)
        axs.plot(mean_fpr, mean_tpr, color='b', label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc), lw=2, alpha=.8)
        axs.plot(np.logspace(-5, 0, 1000), np.logspace(-5, 0, 1000), linestyle='--', lw=2, color='k', label='Random chance')

        std_tpr = np.std(tprs, axis=0)
        tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
        tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
        axs.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2, label=r'$\pm$ 1 std. dev.')
        axs.set_xscale('log')
        axs.set_xlim([1.0e-5, 1.0])
        axs.set_ylim([0.0, 1.0])
        axs.set_xlabel('False Alarm Rate')
        axs.set_ylabel('Signal Efficiency')
        axs.set_title('Cross-validation Receiver Operating Curve')
        axs.legend(loc="lower right") 
        figs.savefig(pdir + '/training_results_roc_cv_{}.pdf'.format(suffix), bbox_inches='tight')


    xgboost_params = best_params.copy()

    # Re-train the whole dataset with the best hyper-parameters (without doing any cross validation)
    print('Training full model...')
    model, results = train(dmatrix_train, dmatrix_test, hyper_params=xgboost_params)
  
    # We want to know if and when the training was early stopped.
    # `best_iteration` counts the first iteration as zero, so we increase by one.
    best_iteration = model.best_iteration + 1
    if best_iteration < n_boost_rounds:
        print("Final model: early stopping after {0} boosting rounds".format(best_iteration))
    print("")
    
    model.save_model(mdir + "/xgb_fulldata_{}.model".format(suffix))

    df.loc[idx_train, "score"] = model.predict(dmatrix_train, ntree_limit=model.best_ntree_limit)
    df.loc[idx_test, "score"] = model.predict(dmatrix_test, ntree_limit=model.best_ntree_limit)
   
    df.loc[idx_train, "test"] = False
    df.loc[idx_test, "test"] = True

    print("")
    print("Final model: Best hyper-parameters: {}, ntree_limit: {}".format(best_params, model.best_ntree_limit))
    print("")

    #df_train = df.query("not test")
    #df_test = df.query("test")
    df_train = df[np.logical_not(df['test'])]
    df_test = df[df['test']]
    #df_test.to_csv('training_results_testdf_{}.csv'.format(suffix))

    fpr, tpr, thresholds = roc_curve(df_test["isSignal"], df_test["score"], drop_intermediate=True)
    roc_dict = {'fpr': fpr, 'tpr': tpr, 'thresholds': thresholds}
    roc_df = pd.DataFrame(data=roc_dict)
    roc_df.to_csv(mdir + '/training_results_roc_csv_{}.csv'.format(suffix))

    epochs = len(results['train']['pauc'])
    x_axis = range(0, epochs)
    fig, ax = plt.subplots()
    ax.plot(x_axis, results['train']['pauc'], label='Train')
    ax.plot(x_axis, results['eval']['pauc'], label='Test')
    ax.legend()
    plt.ylabel('PAUC')
    plt.xlabel('Epoch')
    fig.savefig(pdir + '/training_results_learning_curve_{}.pdf'.format(suffix), bbox_inches='tight')


    fig, ax = plt.subplots()
    plot_roc_curve(df_test, "score", ax=ax, label="XGB")
    ax.plot(np.logspace(-5, 0, 1000), np.logspace(-5, 0, 1000), linestyle='--', color='k')
    ax.set_xlim([1.0e-5, 1.0])
    ax.set_ylim([0.0, 1.0])
    ax.set_xscale('log')
    ax.set_xlabel("False Alarm Rate")
    ax.set_ylabel("Signal Efficiency")
    ax.set_title('Receiver Operating Curve')
    ax.legend(loc='lower right')
    fig.savefig(pdir + '/training_results_roc_curve_{}.pdf'.format(suffix), bbox_inches='tight')

    plt.figure()
    xgb.plot_importance(model)
    plt.savefig(pdir + '/training_results_feature_importance_{}.pdf'.format(suffix), bbox_inches='tight')

    ###################### plotting ###################
#    import pdb; pdb.set_trace()
#####    with PdfPages('Plots/dist.pdf') as pdf:
#####      for var in features:
#####        #if var != "ele_pt": continue
#####        print("plotting {}...".format(var))
#####        #fig, axes = plt.subplots()
#####        fig, axes = plt.subplots()
######        plot_hist(df.sample(n=1000000), var, bins=bins, ax=axes[0], title='All')
######        plot_hist(df.sample(n=1000000), var, bins=np.linspace(0, 5,100), ax=axes[0], title='All')
######        axes[0].set_xlabel('alpha')
######        axes[0].set_ylabel('a.u.')
#####        plot_hist(df, var, bins=vardir[var]['bin'], ax=axes, title='')
######        plot_hist(df[pf_selection], var, bins=bins, ax=axes[1], title='PF-PF')
######        plot_hist(df[low_selection].sample(n=1000000), var, bins=bins, ax=axes[2], title='LowPt-LowPt')
#####        #fig, axes = plt.subplots(1, 2, figsize=(15, 5))
#####        #plot_hist(df.query("abs(BToKEE_fit_eta) <= 1.57"), var, bins=bins, ax=axes[0], title="Barrel")
#####        #plot_hist(df, var, bins=bins, ax=axes[1], title="All")
#####        pdf.savefig(fig, bbox_inches='tight')
#####        plt.close()
#####        print("finished plotting {}...".format(var))

#    training_features = ['BToKEE_fit_l1_normpt', 'BToKEE_fit_l1_eta', 'BToKEE_l1_dxy_sig',
#              'BToKEE_fit_l2_normpt', 'BToKEE_fit_l2_eta', 'BToKEE_l2_dxy_sig',
#              'BToKEE_fit_k_normpt', 'BToKEE_fit_k_eta', 'BToKEE_k_DCASig',
#              'BToKEE_fit_normpt', 'BToKEE_fit_eta', 'BToKEE_svprob', 'BToKEE_fit_cos2D', 'BToKEE_l_xy_sig', 'BToKEE_dz',
#              ]
#    training_features += ['BToKEE_minDR', 'BToKEE_maxDR']
#    training_features += ['BToKEE_l1_iso04_rel', 'BToKEE_l2_iso04_rel', 'BToKEE_k_iso04_rel', 'BToKEE_b_iso04_rel']
#    training_features += ['BToKEE_l1_pfmvaId_lowPt', 'BToKEE_l2_pfmvaId_lowPt', 'BToKEE_l1_pfmvaId_highPt', 'BToKEE_l2_pfmvaId_highPt']
#    training_features += ['BToKEE_ptImbalance']
    #training_features += ['BToKEE_l1_mvaId', 'BToKEE_l2_mvaId']
#    plot_corr(df.query('Category == 0')[features], args.outputfile.replace('.pdf','') + '_corr_sig.pdf')
#    plot_corr(df.query('Category == 1')[features], args.outputfile.replace('.pdf','') + '_corr_bkg.pdf')
    
#    plot_pairgrid(df[features+['Category']].sample(n=1000), 'Plots/pairgrid.pdf')
