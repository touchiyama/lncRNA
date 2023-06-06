import warnings
warnings.filterwarnings('ignore')

import os
import sys

system = sys.version_info
if system[0] < 3 or system[1] < 9:
    sys.stderr.write("ERROR: Your python version is: "+ str(system[0]) + "." + str(system[1]) +"\n" + "LncDC requires python 3.9 or newer!\n" )
    sys.exit(1)


import re
import sys
import gzip
import pathlib
import pickle
import subprocess
import numpy as np
import pandas as pd
import argparse
import xgboost as xgb
import matplotlib.pyplot as plt
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from imblearn.over_sampling import SMOTE
from imblearn.under_sampling import RandomUnderSampler
from sklearn.model_selection import KFold
from sklearn.model_selection import RandomizedSearchCV
import train_SIF_PF_extraction

seed = 666

def command(cmd):
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f'Failed: {e}')
        print('Exit..')
        sys.exit(1)

def pred_eval(y, y_hat, title, outfile, plotflg=False):
    """RMSEと相関係数を計算

    Args:
        y: 真値
        y_hat: 予測値
        title (_type_): _description_
        plotflg: プロットのOn/ Off(デフォルト: False)

    Returns:
        rmse: RMSE
        r: 相関係数
    """
    rmse = np.linalg.norm(y - y_hat)/np.sqrt(len(y))
    r = np.corrcoef(y, y_hat)[0,1]

    if plotflg:
        fig, ax = plt.subplots()
        plt.title(f'{title}')
        plt.xlabel('Prediction(y_test_model)')
        plt.ylabel('Reference(y_test)')
        plt.scatter(y, y_hat, color='black')

        xmin, xmax = ax.get_xlim()
        ymin, ymax = ax.get_ylim()

        plt.plot([xmin,xmax], [ymin,ymax], color='gray', linestyle='dashed')

        r_text = f'r={r:.2f}'
        rmse_text = f'rmse={rmse :.2f}'

        posx = (xmax - xmin)*0.01 + xmin
        posy_1 = (ymax - ymin)*0.95 + ymin
        posy_2 = (ymax - ymin)*0.90 + ymin
        ax.text(posx, posy_1, r_text)
        ax.text(posx, posy_2, rmse_text)

        plt.savefig(outfile, bbox_inches='tight', dpi=300)

    return rmse, r

def calc_score(pred, target):
    """分類性能の評価

    Args:
        pred (_type_): _description_
        target (_type_): _description_

    Returns:
        _type_: _description_
    """
    # 真陽性、真陰性、偽陽性、偽陰性
    TP, TN, FP, FN = 0, 0, 0, 0

    for i in range(0, len(pred)):
        if pred[i] == 0 and target[i] == 0:
            TN +=1
        if pred[i] == 1 and target[i] == 1:
            TP +=1
        if pred[i] == 0 and target[i] == 1:
            FN +=1
        if pred[i] == 1 and target[i] == 0:
            FP +=1

    if (TP + FN != 0) and (FP + TN != 0):
        sensitivity = TP / (TP + FN)
        spcificity = TN / (FP + TN)
        # PPV = TP / (TP + FP)
    else:
        if (TP + FN  == 0) and (FP + TN == 0):
            sensitivity = 0
            spcificity = 0
        elif (TP + FN  == 0) and (FP + TN != 0):
            sensitivity = 0
            spcificity = TN / (FP + TN)
        else:
            sensitivity = TP / (TP + FN)
            spcificity  = 0

    g_mean = np.sqrt(sensitivity*spcificity)

    return sensitivity, spcificity, g_mean

def under_over_process(x, y, njobs):
    ratio = y.value_counts()[0] / y.value_counts()[1]
    # undersample first
    if(ratio < 2/3 or ratio > 1.5):
        rus = RandomUnderSampler(sampling_strategy=2/3,random_state=seed)
        x_underSampled, y_underSampled = rus.fit_resample(x, y)
    else:
        x_underSampled = x
        y_underSampled = y

    # then oversample
    smote = SMOTE(random_state=seed, n_jobs = njobs)
    x_resampled, y_resampled = smote.fit_resample(x_underSampled, y_underSampled)
    return x_resampled, y_resampled

def XGB_performance(df, njobs, n_split=10, n_cv=10):
    """
    # 不均衡データの場合は、学習用データと検証用データに入っている少数クラスの組み合わせによって、
    性能が大きく変化する。そこでデータセットをランダムに10回組み替えて、モデル学習と性能評価を繰り返す。
    最終的な性能は、10回の計算での感度、特異度、G-meanそれぞれの平均値と標準偏差で評価。


    Args:
        df (_type_): _description_
        njobs (_type_): _description_
        n_split (int, optional): _description_. Defaults to 10.
        n_cv (int, optional): _description_. Defaults to 10.

    Returns:
        _type_: _description_
    """
    # datasetの準備 ---
    # 学習用データと検証用データを10回ランダムに組み替える
    kf = KFold(n_splits=n_split, shuffle=True, random_state=seed)

    # 正解ラベルの抽出 ---
    Y = df.iloc[:, -1]
    X = df.drop(df.columns[-1], axis=1)

    # データの補完 ---
    print("Imputing missing values ...")
    # imputation of missing values
    imputer_mean = SimpleImputer(missing_values=np.nan, strategy='mean')
    imputer_mean.fit(X[X.columns])
    X[X.columns] = imputer_mean.transform(X[X.columns])
    # write imputer to a pickle file
    # pkl_filename = output_prefix + '_imputer_SIF_PF.pkl'
    # print("Writing the imputer to file: " + str(pkl_filename))

    # 学習と検証の繰り返し ---
    for train_index, test_index in kf.split(X):
        X_train = X.iloc[train_index, :].values
        Y_train = Y.iloc[train_index]
        X_test = X.iloc[test_index, :].values
        Y_test = Y.iloc[test_index]

    # SMOTEサンプリング ---
    X_train, Y_train = under_over_process(X_train, Y_train, njobs)

    # データの標準化 ---
    # サンプリング後の学習用データセットにおける基本統計量を用いて検証用データセットを標準化することで、
    # モデル学習後の性能評価で学習用と検証用のデータ分布の乖離を防ぐ
    mean = np.mean(X_train, axis = 0)
    std = np.std(X_train, axis = 0, ddof = 1)
    X_train = (X_train - mean) / std
    X_test = (X_test - mean) / std

    # モデルの学習、ランダムリサーチ
    xgb_model = xgb.XGBClassifier(
        objective='binary:logistic',
        random_state=seed,
        eval_metric='logloss'
    )

    # https://qiita.com/c60evaporator/items/a9a049c3469f6b4872c6
    rand_param = {
        "colsample_bytree": uniform(0.7, 0.3),
        "gamma": uniform(0, 0.5),
        "learning_rate": uniform(0.03, 0.3), # default 0.1
        "max_depth": randint(2, 6), # default 3
        "n_estimators": randint(100, 150), # default 100
        "subsample": uniform(0.6, 0.4)
    }
    rand_search = RandomizedSearchCV(
        xgb_model, # model selection
        rand_param,
        scoring='neg_root_mean_squared_error',
        cv=n_cv,
        n_jobs=njobs
    )
    rand_search.fit(X_train, Y_train)

    # ベストモデルのハイパーパラメータ ---
    XGB_best = rand_search.best_estimator_
    XGB_best.fit(X_train, Y_train)

    # 検証用データのラベル予測 ---
    Y_pred = XGB_best.predict(X_test)

    # 性能評価 ---
    sense, spec, g_mean = calc_score(Y_pred, Y_test)

    # 平均値と標準偏差を求める
    # ハイパーパラメータも平均値か？ <- rand_search.best_params_の平均値

    # return sense, spec, g_mean



def RS_XGBoost(X_train, y_train, X_test, y_test, outfile, grid_param, n_cv=10):
    """ランダムリサーチの実行

    Args:
        X_train (_type_): _description_
        y_train (_type_): _description_
        X_test (_type_): _description_
        y_test (_type_): _description_
        outfile (_type_): _description_
        grid_param (_type_): _description_
        n_cv (int, optional): _description_. Defaults to 10.

    Returns:
        _type_: _description_
    """

    xgb_model = xgb.XGBClassifier(
        objective='binary:logistic',
        random_state=seed,
        eval_metric='logloss'
    )

    params = {
        "colsample_bytree": uniform(0.7, 0.3),
        "gamma": uniform(0, 0.5),
        "learning_rate": uniform(0.03, 0.3), # default 0.1
        "max_depth": randint(2, 6), # default 3
        "n_estimators": randint(100, 150), # default 100
        "subsample": uniform(0.6, 0.4)
    }

    rand_search = RandomizedSearchCV(
        xgb_model, # model selection
        grid_param,
        scoring='neg_root_mean_squared_error',
        cv=n_cv,
        n_jobs=-1
    )
    rand_search.fit(X_train, y_train)

    """
    name = os.path.splitext(os.path.basename(outfile))[0]
    with open(outfile, 'w') as wf:
        logger.info(f'@{name}')
        wf.write(f'@{name}\n')
        logger.info('Accuracy Score (train) : ', (-1) * rand_search.score(X_train, y_train))
        wf.write(f'Accuracy Score (train) : {(-1) * rand_search.score(X_train, y_train)}\n')
        logger.info('Accuracy Score (test) : ', (-1) * rand_search.score(X_test, y_test))
        wf.write(f'Accuracy Score (test) : {(-1) * rand_search.score(X_test, y_test)}\n')
        logger.info('Best Parameters : ', rand_search.best_params_)
        wf.write(f'Best Parameters : {rand_search.best_params_}\n')
        logger.info('Best CV Score : ', (-1) * rand_search.best_score_) # trainデータに対するCVの平均CV精度
        wf.write(f'Best CV Score : {(-1) * rand_search.best_score_}\n')
        logger.info('Best Estimator : ', rand_search.best_estimator_)
        wf.write(f'Best Estimator : {rand_search.best_estimator_}\n')
    """

    pls_best = rand_search.best_estimator_

    # 回帰モデルの評価 ---
    pls = pls_best.fit(X_train, y_train)
    pls_coef = np.hstack(pls.coef_)
    y_test_pls = X_test @ pls_coef
    outpng = os.path.splitext(outfile)[0] + '_eval.png'

    # gmeansで評価 ---
    #rmse, r = pred_eval(y_test, y_test_pls, name, outpng, plotflg=True)

    return pls_coef, rmse, r

def file_exist(filename,parser):
    if filename == None:
        parser.print_help()
        sys.exit(1)
    else:
        if not os.path.isfile(filename):
            sys.stderr.write("ERROR: No such file: " + filename + "\n")
            sys.exit(1)
        else:
            print("File "+ filename + " exist.")

def fasta_format_check(filename):
    if filename.endswith((".gz", ".Z", ".z")):
        fd = gzip.open(filename, 'rt')
    else:
        fd = open(filename, 'r')

    filetype = 1
    for line in fd:
        if line.startswith('#'):
            continue
        elif line.strip() == None:
            continue
        elif line.startswith('>'):
            filetype = 0
            break
        else:
            filetype = 1
            break

    fd.close()
    if not filetype == 0:
        sys.stderr.write("ERROR: Your " + filename + "is not in fasta format \n")
        sys.exit(1)
    else:
        print(filename + " format checking: PASS")

def get_splice_site_score(five_ss_score_file, three_ss_score_file):
    """_summary_

    Args:
        five_ss_score_file (_type_): _description_
        three_ss_score_file (_type_): _description_
        # file = '/Users/tomoyauchiyama/lncRNA/maxentpy/gencode.v43.lncRNA/gencode_v43_lncRNA_score5.txt'
        # file = '/Users/tomoyauchiyama/lncRNA/maxentpy/gencode.v43.lncRNA/gencode_v43_lncRNA_score3.txt'
    Returns:
        _type_: _description_
    """
    df5 = pd.read_csv(five_ss_score_file, sep='\t', header=None)
    df_split = df5[0].str.split('::', expand=True)
    df5 = pd.concat([df_split, df5[[1, 2]]], axis=1)
    df5.columns = ['transcript_id', '5ss_pos', '5ss_seq', '5ss_score']

    df3 = pd.read_csv(three_ss_score_file, sep='\t', header=None)
    df_split = df3[0].str.split('::', expand=True)
    df3 = pd.concat([df_split, df3[[1, 2]]], axis=1)
    df3.columns = ['transcript_id', '3ss_pos', '3ss_seq', '3ss_score']

    df_ss_score = pd.merge(df5, df3, how='left', on='transcript_id')

    return df_ss_score

def get_fasta(gtf_file, ref_genome, outdir, gffread_path, CDS):
    fname = os.path.splitext(os.path.basename(gtf_file))[0]
    if CDS:
        fasta = os.path.join(outdir, fname + '_cds.fa')
        cmd = [gffread_path, gtf_file, '-g', ref_genome, '-x', fasta]
        command(cmd)
    else:
        fasta = os.path.join(outdir, fname + '.fa')
        cmd = [gffread_path, gtf_file, '-g', ref_genome, '-w', fasta]
        command(cmd)

    return fasta

def load_fasta(filename):
    '''
    load_fasta
    Inputs:
        filename - name of FASTA file to load
    Returns:
        a list of sequences
    '''

    seq = {}
    #sequences = []
    #seq = ''
    if filename.endswith((".gz", ".Z", ".z")):
        fd = gzip.open(filename, 'rt')
    else:
        fd = open(filename, 'r')

    for line in fd:
        if line.startswith('>'):
            transcript_id = re.search('>(.*)', line.strip().split()[0]).group(1)
            #if seq != '':
                #sequences.append(seq.replace('a','A').replace('t','T').replace('g','G').replace('c','C'))
            #seq = ''
        else:
            #seq += line.strip()
            if seq.get(transcript_id):
                seq[transcript_id] += line.strip()
            else:
                seq[transcript_id] = line.strip()
    #if seq != '':
        #sequences.append(seq.replace('a','A').replace('t','T').replace('g','G').replace('c','C'))
    fd.close()

    df_seq = pd.DataFrame([seq.keys(), seq.values()], index=['transcript_id', 'Sequence']).T

    return df_seq

def main():
    parser = argparse.ArgumentParser(
        description='LncDC: a machine learning based tool for long non-coding RNA detection from RNA-Seq data'
    )
    parser.add_argument(
        '-v','--version', action = 'version', version = '%(prog)s version:1.3.5'
    )
    parser.add_argument(
        '-m','--mrna',
        help = 'The file with mRNA sequences in fasta format. The fasta file could be regular text file or gzip compressed file (*.gz).',
        type = str, required = True, default = None
    )
    parser.add_argument(
        '-c','--cds',
        help = 'The CDS sequences of the mRNAs in fasta format. The fasta file could be regular text file or gzip compressed file (*.gz). The order and number of the CDS sequences should be the same as the mRNA sequences.',
        type = str, required = True, default = None
    )
    parser.add_argument(
        '-l','--lncrna',
        help = 'The file with lncRNA sequences in fasta format. The fasta file could be regular text file or gzip compressed file (*.gz).',
        required = True, default = None
    )
    parser.add_argument(
        '-s','--splicescore',
        help = 'The file with splice score in both mRNA and lncRNA.',
        required = True, default = None
    )
    parser.add_argument(
        '-o','--output',
        help = 'The prefix of the output files, including a hexamer table, a prediction model, an imputer and a scaler. If the -r parameter turned on, the output files will also include secondary structure kmer tables.',
        type = str, default = 'self_train'
    )
    parser.add_argument(
        '-r','--secondary',
        help = '(Optional) Turn on to train a model with secondary structure features. This will generate secondary structure kmer tables. Default: turned off.',
        action = 'store_true'
    )
    parser.add_argument(
        '-t','--thread',
        help = '(Optional) The number of threads assigned to use. Set -1 to use all cpus. Default value: -1.',
        type = int, default = -1
    )
    parser.add_argument(
        '-f','--fraction',
        help = '(Optional) The fraction of random sampling. Set 1 to use all dataset. Default value: 1.',
        type = float, default = 1
    )

    args = parser.parse_args()
    mrna = args.mrna
    cds = args.cds
    lncrna = args.lncrna
    splicescore = args.splicescore
    ss_feature = args.secondary
    output_prefix = args.output
    frac = args.fraction
    thread = args.thread

    if thread == -1:
        thread = os.cpu_count()
    else:
        thread = thread

    if ss_feature:
        try:
            import RNA
        except:
            sys.stderr.write("ViennaRNA is not installed! \n")
            sys.stderr.write("ViennaRNA is required for secondary structure feature extraction. \nYou can install it by conda: conda install -c bioconda viennarna \n")
            sys.exit(1)

    print("Process Start.")
    print("Checking if the training files exist ...")
    file_exist(mrna, parser)
    file_exist(cds, parser)
    file_exist(lncrna, parser)

    print()
    print("Checking if the input files are in gtf format ...")
    #format_check(mrna)
    #format_check(cds)
    #format_check(lncrna)

    # create any parent directories if needed
    filepath = pathlib.Path().resolve()
    output_prefix = os.path.join(filepath, output_prefix)

    if not os.path.isfile(output_prefix):
        os.makedirs(os.path.dirname(output_prefix), exist_ok = True)

    # ['transcript_id', '5ss_seq', '5ss_score', '3ss_seq', '3ss_score']
    df_ss_score = pd.read_csv(splicescore, sep='\t')
    df_ss_score = df_ss_score.drop(columns=['5ss_seq', '3ss_seq'])

    """
    # conduct gffread manually each other ---
    mrna = get_fasta(gtf_file, ref_genome, outdir, gffread_path, False)
    cds = get_fasta(gtf_file, ref_genome, outdir, gffread_path, False)
    lncrna = get_fasta(gtf_file, ref_genome, outdir, gffread_path, False)
    """

    print()
    print("Checking if the training files are in fasta format ...")
    fasta_format_check(mrna)
    fasta_format_check(cds)
    fasta_format_check(lncrna)

    # load mrna data
    mrna_data = load_fasta(mrna)
    # load cds data
    cds_data = load_fasta(cds)

    print()
    print("Initializing dataframe ...")
    # initialize a mrna dataframe
    mrna_dataset = pd.merge(mrna_data, cds_data, how='left', on='transcript_id')
    # rename column name
    mrna_dataset = mrna_dataset.rename(columns={'Sequence_x': 'Sequence', 'Sequence_y': 'CDS_seq'})
    mrna_dataset.loc[:, 'type'] = 'mrna'

    # load lncrna data
    lncrna_dataset = load_fasta(lncrna)
    lncrna_dataset.loc[:, 'CDS_seq'] = 0
    lncrna_dataset.loc[:, 'type'] = 'lncrna'

    """"
    mrna_dataset = pd.DataFrame(index=range(len(mrna_data)), columns=['Sequence', 'type', 'CDS_seq'])
    for i in range(mrna_dataset.index.size):
        mrna_dataset.loc[i, 'Sequence'] = mrna_data[i]
        mrna_dataset.loc[i,'type'] = 'mrna'
        mrna_dataset.loc[i, 'CDS_seq'] = cds_data[i]

    # load lncrna data
    lncrna_data = load_fasta(lncrna)

    # initialize a lncrna dataframe
    lncrna_dataset = pd.DataFrame(index=range(len(lncrna_data)), columns=['Sequence', 'type', 'CDS_seq'])
    for i in range(lncrna_dataset.index.size):
        lncrna_dataset.loc[i, 'Sequence'] = lncrna_data[i]
        lncrna_dataset.loc[i,'type'] = 'lncrna'
        lncrna_dataset.loc[i, 'CDS_seq'] = 0
    """

    # combine to a single dataframe
    dataset_seq = pd.concat([mrna_dataset, lncrna_dataset])
    # reset the index of the dataframe
    dataset_seq.reset_index(drop=True, inplace=True)
    # join splicing score -----
    dataset = pd.merge(dataset_seq, df_ss_score, how='left', on=['transcript_id', 'type'])

    print("Total Number of transcripts loaded: " + str(dataset.index.size))

    # remove the sequences that have non [ATGC] inside
    for i in range(dataset.index.size):
        if len(re.findall(r'[^ATGC]',dataset.loc[i,'Sequence'])) > 0:
            dataset.loc[i,'Sequence'] = float('NaN')
    dataset.dropna(how = 'any', inplace = True)
    # reset the index of the dataframe
    dataset.reset_index(drop = True, inplace = True)

    print("Calculating transcript lengths ...")
    print()
    # Calculate the length of the transcripts
    for i in range(dataset.index.size):
        dataset.loc[i,'Transcript_length'] = len(dataset.loc[i, 'Sequence'])


    # IF only use SIF and PF features
    if not ss_feature:
        # Filter out sequence length less than 200nt or more than 20000nt
        print()
        dataset = dataset[(dataset['Transcript_length'] >= 200) & (dataset['Transcript_length'] <= 20000)]
        print("Removing Non-valid transcripts (sequence that have non-ATGCatgc letters & sequence length less than 200 nt) ...")
        dataset = dataset.reset_index(drop=True)
        print("Number of valid transcripts for training before random sampling: " + str(dataset.index.size))
        dataset = dataset.sample(frac=frac, random_state=seed).reset_index(drop=True)
        print("Number of valid transcripts for training after random sampling: " + str(dataset.index.size))

        if dataset.index.size == 0:
            sys.stderr.write("No valid transcripts detected! \n")
            sys.exit(1)

        print()
        print("Extracting SIF and PF features ...")

        # extract features
        # drop transcript id, add 5' ss score and 3' ss score to dataset
        dataset = train_SIF_PF_extraction.sif_pf_extract(dataset, thread, output_prefix)
        columns = [
            'type','Transcript_length', '5ss_score', '3ss_score', 'GC_content', 'Fickett_score',
            'ORF_T0_length', 'ORF_T1_length', 'ORF_T2_length','ORF_T0_coverage', 'ORF_T1_coverage', 'ORF_T3_coverage',
            'Hexamer_score_ORF_T0', 'Hexamer_score_ORF_T1', 'Hexamer_score_ORF_T2', 'Hexamer_score_ORF_T3',
            'RCB_T0', 'RCB_T1', 'ORF_T0_PI', 'ORF_T0_MW', 'ORF_T0_aromaticity', 'ORF_T0_instability',
            'ORF_T1_MW', 'ORF_T1_instability', 'ORF_T2_MW', 'ORF_T3_MW'
        ]
        dataset = dataset[columns]
        x_train = dataset.drop(['type','Transcript_length'], axis = 1)
        y_train = dataset[['type']]

        print("Imputing missing values ...")
        # imputation of missing values
        imputer_mean = SimpleImputer(missing_values=np.nan, strategy='mean')
        imputer_mean.fit(x_train[x_train.columns])
        x_train[x_train.columns] = imputer_mean.transform(x_train[x_train.columns])

        # write imputer to a pickle file
        pkl_filename = output_prefix + '_imputer_SIF_PF.pkl'
        print("Writing the imputer to file: " + str(pkl_filename))

        with open(pkl_filename, 'wb') as file:
            pickle.dump(imputer_mean, file)

        print("Standardizing ...")
        # standardization
        scaler = StandardScaler()
        scaler.fit(x_train[x_train.columns])
        x_train[x_train.columns] = scaler.transform(x_train[x_train.columns])
        pkl_filename = output_prefix + '_scaler_SIF_PF.pkl'
        print("Writing the scaler to file: " + str(pkl_filename))
        with open(pkl_filename, 'wb') as file:
            pickle.dump(scaler, file)

        # build standardized training dataset
        train_std = pd.concat([x_train, y_train], axis=1)

        print("Balancing training data ...")
        # balancing training data
        x_resampled, y_resampled = under_over_process(
            train_std.drop(['type'], axis=1),
            train_std['type'], thread)

        y_resampled = y_resampled.map({'mrna':1,'lncrna':0})

        # RandomReseach ----
        #


        print("Model fitting ...")
        # fit a model
        xgb_model = xgb.XGBClassifier(alpha = 10, colsample_bytree = 0.8, gamma = 0.01, learning_rate = 0.1,
                                     max_depth = 9, min_child_weight = 4, n_estimators = 1000,
                                     objective = 'binary:logistic', subsample = 0.8, verbosity = 0,
                                     random_state = seed, n_jobs = thread)
        xgb_model.fit(x_resampled, y_resampled)
        pkl_filename = output_prefix + '_xgb_model_SIF_PF.pkl'
        print("Writing the model to file: " +str(pkl_filename))
        with open(pkl_filename, 'wb') as file:
            pickle.dump(xgb_model, file)

        print("Done!")

    else:
        # Use SIF + PF + SSF
        import train_SSF_extraction
        print()
        # Filter out sequence length less than 200nt or more than 20000nt
        dataset = dataset[(dataset['Transcript_length'] >= 200) & (dataset['Transcript_length'] <= 20000)]
        print("Removing Non-valid transcripts (sequence that have non-ATGCatgc" + " letters & sequence length less than 200 nt) ...")
        dataset = dataset.reset_index(drop=True)
        print("Number of valid transcripts for training before random sampling: " + str(dataset.index.size))
        dataset = dataset.sample(frac=frac, random_state=seed).reset_index(drop=True)
        print("Number of valid transcripts for training after random sampling: " + str(dataset.index.size))

        if dataset.index.size == 0:
            sys.stderr.write("No valid transcripts detected! \n")
            sys.exit(1)

        print()
        print("Extracting SIF and PF features ...")

        # extract features
        # drop transcript id, add 5' ss score and 3' ss score to dataset
        dataset = train_SIF_PF_extraction.sif_pf_extract(dataset, thread, output_prefix)
        columns = [
            'type', 'Sequence','Transcript_length', '5ss_score', '3ss_score', 'GC_content', 'Fickett_score',
            'ORF_T0_length', 'ORF_T1_length', 'ORF_T2_length', 'ORF_T0_coverage', 'ORF_T1_coverage', 'ORF_T3_coverage',
            'Hexamer_score_ORF_T0', 'Hexamer_score_ORF_T1', 'Hexamer_score_ORF_T2', 'Hexamer_score_ORF_T3',
            'RCB_T0', 'RCB_T1', 'ORF_T0_PI', 'ORF_T0_MW', 'ORF_T0_aromaticity', 'ORF_T0_instability',
            'ORF_T1_MW', 'ORF_T1_instability', 'ORF_T2_MW', 'ORF_T3_MW'
        ]
        dataset = dataset[columns]

        # extract Secondary Structure Features
        # drop transcript id, add 5' ss score and 3' ss score to dataset
        dataset = train_SSF_extraction.ssf_extract(dataset, thread, output_prefix)
        full_columns = [
            'type', 'Transcript_length', '5ss_score', '3ss_score', 'GC_content', 'Fickett_score',
            'ORF_T0_length', 'ORF_T1_length','ORF_T2_length', 'ORF_T0_coverage', 'ORF_T1_coverage', 'ORF_T3_coverage',
            'Hexamer_score_ORF_T0','Hexamer_score_ORF_T1', 'Hexamer_score_ORF_T2', 'Hexamer_score_ORF_T3',
            'RCB_T0','RCB_T1', 'SS_score_k1','SS_score_k2','SS_score_k3','SS_score_k4',
            'SS_score_k5','GC_content_paired_ss',
            'ORF_T0_PI', 'ORF_T0_MW', 'ORF_T0_aromaticity', 'ORF_T0_instability',
            'ORF_T1_MW', 'ORF_T1_instability', 'ORF_T2_MW', 'ORF_T3_MW'
        ]
        dataset = dataset[full_columns]
        x_train = dataset.drop(['type', 'Transcript_length'], axis = 1)
        y_train = dataset[['type']]

        print("Imputing missing values ...")
        # imputation of missing values
        imputer_mean = SimpleImputer(missing_values=np.nan, strategy='mean')
        imputer_mean.fit(x_train[x_train.columns])
        x_train[x_train.columns] = imputer_mean.transform(x_train[x_train.columns])

        # write imputer to a pickle file
        pkl_filename = output_prefix + '_imputer_SIF_PF_SSF.pkl'
        print("Writing the imputer to file: " + str(pkl_filename))

        with open(pkl_filename, 'wb') as file:
            pickle.dump(imputer_mean, file)

        print("Standardizing ...")
        # standardization
        scaler = StandardScaler()
        scaler.fit(x_train[x_train.columns])
        x_train[x_train.columns] = scaler.transform(x_train[x_train.columns])
        pkl_filename = output_prefix + '_scaler_SIF_PF_SSF.pkl'
        print("Writing the scaler to file: " + str(pkl_filename))
        with open(pkl_filename, 'wb') as file:
            pickle.dump(scaler, file)

        # build standardized training dataset
        train_std = pd.concat([x_train, y_train], axis=1)

        print("Balancing training data ...")
        # balancing training data
        x_resampled, y_resampled = under_over_process(
            train_std.drop(['type'], axis=1),
            train_std['type'], thread)

        y_resampled = y_resampled.map({'mrna':1,'lncrna':0})

        print("Model fitting ...")
        # fit a model
        xgb_model = xgb.XGBClassifier(alpha=10, colsample_bytree=0.8, gamma=0.01, learning_rate=0.1,
                                      max_depth=9, min_child_weight=4, n_estimators=1000,
                                      objective='binary:logistic', subsample=0.8, verbosity = 0,
                                      random_state=seed, n_jobs=thread)
        xgb_model.fit(x_resampled, y_resampled)
        pkl_filename = output_prefix + '_xgb_model_SIF_PF_SSF.pkl'
        print("Writing the model to file: " +str(pkl_filename))
        with open(pkl_filename, 'wb') as file:
            pickle.dump(xgb_model, file)

        print("Done!")

if __name__ == '__main__':
    main()
