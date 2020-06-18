# -*- coding:utf-8 -*-
import pandas as pd
import pyfaidx
import vcf
import re


class SlidingWindow(object):
    def __init__(self, start, stop, win_len, slide_len):
        self.start = start
        self.stop = stop
        self.win_len = win_len
        self.slide_len = slide_len
        self.i_start = self.start
        self.i_stop = self.start + self.win_len - 1

    def __iter__(self):
        return self

    def __next__(self):
        i_start = self.i_start
        i_stop = self.i_stop
        self.i_start = self.i_start + self.slide_len
        self.i_stop = self.i_start + self.win_len - 1
        if i_start >= self.stop or i_stop > self.stop:
            raise StopIteration
        else:
            if self.i_start >= self.stop or self.i_stop >= self.stop:
                return i_start, self.stop
            else:
                return i_start, i_stop


def get_vcf_dp_gt_ratio_flt(vcf_file):
    chrom, pos, ref, alt, depth, genotype, ratio, flt, mutype = [], [], [], [], [], [], [], [], []
    vcf_reader = vcf.Reader(filename=vcf_file)
    for record in vcf_reader:
        try:
            if record.samples[0]['AD']:
                if record.samples[0]['GT'] != '1/2':
                    ratio.append(float(record.samples[0]['AD'][1])/(float(record.samples[0]['AD'][0])+float(record.samples[0]['AD'][1])))
                else:
                    ratio.append(float(record.samples[0]['AD'][2])/(float(record.samples[0]['AD'][1])+float(record.samples[0]['AD'][2])))
                chrom.append(record.CHROM)
                pos.append(record.POS)
                ref.append(record.REF)
                if len(record.ALT) == 1:
                    alt.append(str(record.ALT[0]))
                else:
                    alt.append(','.join([str(record.ALT[0]), str(record.ALT[1])]))
                depth.append(record.INFO['DP'])
                genotype.append(record.samples[0]['GT'])
                if record.FILTER:
                    flt.append('Filter')
                else:
                    flt.append('Pass')
                if record.is_snp:
                    mutype.append("snp")
                elif record.is_indel:
                    mutype.append("indel")
                else:
                    mutype.append("")
        except:
            pass
    dpgt_dict = {'#CHROM': chrom, 'POS': pos, 'REF': ref, 'ALT': alt, 'DP': depth, 'GT': genotype, 'Ratio': ratio, 'FILTER': flt, 'Type': mutype}
    dpgt_df = pd.DataFrame(dpgt_dict, columns=['#CHROM', 'POS', 'REF', 'ALT', 'DP', 'GT', 'Ratio', 'FILTER', 'Type'])
    return dpgt_df


def split_df(df, split_num):
    df.reset_index(drop=True, inplace=True)
    df_list = list()
    step = round(df.shape[0]/split_num)
    for i in range(split_num):
        if i == 0:
            df_list.append(df.loc[0: step-1])
        elif i == split_num-1:
            df_list.append(df.loc[step*i:])
        else:
            df_list.append(df.loc[step*i:step*(i+1)-1])
    return df_list


def bgi_anno_2_vcf_format(in_df, reference):
    df = in_df.copy()
    df['#Chr'] = df['#Chr'].astype('str')
    if len(df[df['#Chr'].str.startswith('chr')]):
        df['#CHROM'] = df['#Chr']
    else:
        df['#CHROM'] = 'chr' + df['#Chr']
    df.loc[df['#CHROM'] == 'chrMT', '#CHROM'] = 'chrM_NC_012920.1'
    df['ID'] = '.'
    df['QUAL'] = '.'
    df['FILTER'] = '.'
    df['INFO'] = '.'
    df['MuType'] = 'delins'
    df.loc[df['Ref'] == '.', 'MuType'] = 'ins'
    df.loc[df['Call'] == '.', 'MuType'] = 'del'
    df.loc[(df['Ref'].map(len) == 1) & (df['Call'].map(len) == 1) & (df['Ref'] != '.') & (df['Call'] != '.'), 'MuType'] = 'snp'
    df['POS'] = df['Stop']
    df.loc[df['MuType'] == 'del', 'POS'] = df.loc[df['MuType'] == 'del', 'Start']
    df.loc[df['MuType'] == 'delins', 'POS'] = df.loc[df['MuType'] == 'delins', 'Start']
    df['REF'] = df['Ref']
    df['ALT'] = df['Call']
    fa = pyfaidx.Fasta(reference)
    for i in range(df.shape[0]):
        if df.loc[i, 'MuType'] == 'ins':
            base = str(fa.get_seq(df.loc[i, '#CHROM'], df.loc[i, 'POS'], df.loc[i, 'POS'])).upper()
            df.loc[i, 'REF'] = base
            df.loc[i, 'ALT'] = base + df.loc[i, 'ALT']
        elif df.loc[i, 'MuType'] == 'del':
            base = str(fa.get_seq(df.loc[i, '#CHROM'], df.loc[i, 'POS'], df.loc[i, 'POS'])).upper()
            df.loc[i, 'ALT'] = base
            df.loc[i, 'REF'] = base + df.loc[i, 'REF']
        elif df.loc[i, 'MuType'] == 'delins':
            base = str(fa.get_seq(df.loc[i, '#CHROM'], df.loc[i, 'POS'], df.loc[i, 'POS'])).upper()
            df.loc[i, 'REF'] = base + df.loc[i, 'REF']
            df.loc[i, 'ALT'] = base + df.loc[i, 'ALT']
        else:
            pass
    a = df[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']].copy()
    a.sort_values(by=['#CHROM', 'POS'], ascending=True, inplace=True)
    b = df[['#Chr', 'Start', 'Stop', 'Ref', 'Call', '#CHROM', 'POS', 'REF', 'ALT']].copy()
    df.drop(columns=['ID', 'QUAL', 'FILTER', 'INFO', 'MuType'], inplace=True)
    return a, b, df


def vcf_format_2_bgi_anno(in_df):
    df = in_df.copy()
    df['MuType'] = 'delins'
    df.loc[(df['REF'].str.len() == 1) & (df['ALT'].str.len() == 1), 'MuType'] = 'snv'
    df.loc[(df['REF'].str.len() == 1) & (df['ALT'].str.len() > 1), 'MuType'] = 'ins'
    df.loc[(df['REF'].str.len() > 1) & (df['ALT'].str.len() == 1), 'MuType'] = 'del'
    df['#CHROM'] = df['#CHROM'].astype('str')
    if len(df[df['#CHROM'].str.startswith('chr')]):
        df['#Chr'] = df['#CHROM']
    else:
        df['#Chr'] = 'chr' + df['#CHROM']
    df.loc[df['#CHROM'] == 'chrM_NC_012920.1', '#Chr'] = 'chrMT'
    df['Start'] = df['POS'].astype(int) - 1
    df['Stop'] = df['POS'].astype(int)
    df['Ref'] = df['REF']
    df['Call'] = df['ALT']
    df.loc[df['MuType'] == 'ins', 'Ref'] = '.'
    df.loc[df['MuType'] == 'ins', 'Call'] = df.loc[df['MuType'] == 'ins', 'ALT'].str[1:]
    df.loc[df['MuType'] == 'del', 'Call'] = '.'
    df.loc[df['MuType'] == 'del', 'Ref'] = df.loc[df['MuType'] == 'del', 'REF'].str[1:]
    df.loc[df['MuType'] == 'ins', 'Start'] = df.loc[df['MuType'] == 'ins', 'Stop']
    df.loc[df['MuType'] == 'del', 'Start'] = df.loc[df['MuType'] == 'del', 'Stop']
    df.loc[df['MuType'] == 'del', 'Stop'] = df.loc[df['MuType'] == 'del', 'Stop'] + df.loc[df['MuType'] == 'del', 'Ref'].str.len()
    df.loc[df['MuType'] == 'delins', 'Stop'] = df.loc[df['MuType'] == 'delins', 'Start'] + df.loc[df['MuType'] == 'delins', 'Ref'].str.len()
    df.loc[df['MuType'] == 'delins', 'Start'] = df.loc[df['MuType'] == 'delins', 'POS']
    df.loc[df['MuType'] == 'delins', 'Ref'] = df.loc[df['MuType'] == 'delins', 'REF'].str[1:]
    df.loc[df['MuType'] == 'delins', 'Call'] = df.loc[df['MuType'] == 'delins', 'ALT'].str[1:]
    a = df[['#CHROM', 'POS', 'REF', 'ALT', '#Chr', 'Start', 'Stop', 'Ref', 'Call']].copy()
    df.drop(columns=['MuType'], inplace=True)
    return a, df


def interpretation(spliceai_list, chrom, pos):
    threshold = 0.11
    ds = spliceai_list[0:4]
    dp = spliceai_list[4:8]
    ds = [float(i) for i in ds]
    dp = [int(i) for i in dp]
    ds_dic = {0: 'acceptor gain', 1: 'acceptor loss', 2: 'donor gain', 3: 'donor loss'}
    ds_index = [i for i, j in enumerate(ds) if j >= threshold]
    if len(ds_index) > 0:
        pred = []
        for k in ds_index:
            if dp[k] < 0:
                comment = str(chrom) + ':' + str(pos + dp[k]) + ' (=' + str(pos) + str(dp[k]) + ') ' + ds_dic[k] + ' ' + str(ds[k])
            else:
                comment = str(chrom) + ':' + str(pos + dp[k]) + ' (=' + str(pos) + '+' + str(dp[k]) + ') ' + ds_dic[k] + ' ' + str(ds[k])
            pred.append(comment)
        pred_content = ';'.join(pred)
        return pred_content
    else:
        return '.'


def judge_stream(difference):
    if difference >= 0:
        return 'upstream'
    else:
        return 'downstream'


def biology(spliceai_list):
    comment = '.'
    threshold = 0.2
    ds = spliceai_list[0:4]
    dp = spliceai_list[4:8]
    ds = [float(i) for i in ds]
    dp = [int(i) for i in dp]
    ds_tuple = [(i, x) for i, x in enumerate(ds)]
    ds_tuple.sort(key=lambda x: x[1], reverse=True)
    ds_index = [ds_tuple[0][0], ds_tuple[1][0]]
    if ds_tuple[0][1] >= threshold and ds_tuple[1][1] >= threshold:
        if {0, 1}.issubset(set(ds_index)):
            comment = "Alternative 3' ss usage. Use of a cryptic site {} nt {} from 3' ss".format(dp[0]-dp[1], judge_stream(dp[0]-dp[1]))
        elif {2, 3}.issubset(set(ds_index)):
            comment = "Alternative 5' ss usage. Use of a cryptic site {} nt {} from 5' ss".format(dp[0] - dp[1], judge_stream(dp[0]-dp[1]))
        elif {1, 3}.issubset(set(ds_index)):
            comment = 'Exon skipped'
        elif {0, 2}.issubset(set(ds_index)):
            comment = 'Intron retention'
    elif ds_tuple[0][1] < threshold and ds_tuple[1][1] < threshold:
        comment = 'No change'
    return comment


def spliceai_interpre(spliceai_res):
    spliceai_res_na = spliceai_res[(spliceai_res['SpliceAI Pred'] == '.') | (spliceai_res['SpliceAI Pred'].isna())].copy()
    spliceai_res_is = spliceai_res[(spliceai_res['SpliceAI Pred'] != '.') & (~spliceai_res['SpliceAI Pred'].isna())].copy()
    spliceai_res_na['SpliceAI Interpretation'] = '.'
    if not spliceai_res_is.empty:
        spliceai_res_is['SpliceAI Interpretation'] = spliceai_res_is['SpliceAI'].str.extract('SpliceAI=(.*?)$')
        spliceai_res_is['SpliceAI Interpretation'] = spliceai_res_is['SpliceAI Interpretation'].str.split(',').str[0]
        spliceai_res_is['SpliceAI Interpretation'] = spliceai_res_is['SpliceAI Interpretation'].str.split('|').str[2:10]
        spliceai_res_is['SpliceAI Interpretation'] = spliceai_res_is.apply(lambda x: interpretation(x['SpliceAI Interpretation'],
                                                                                                    x['#CHROM'], x['POS']), axis=1)
    spliceai_res_final = spliceai_res_is.append(spliceai_res_na, sort=False)
    return spliceai_res_final


def spliceai_biology(spliceai_res):
    spliceai_res_na = spliceai_res[(spliceai_res['SpliceAI Pred'] == '.') | (spliceai_res['SpliceAI Pred'].isna())].copy()
    spliceai_res_is = spliceai_res[(spliceai_res['SpliceAI Pred'] != '.') & (~spliceai_res['SpliceAI Pred'].isna())].copy()
    spliceai_res_na['SpliceAI Biology'] = '.'
    if not spliceai_res_is.empty:
        spliceai_res_is['SpliceAI Biology'] = spliceai_res_is['SpliceAI'].str.extract('SpliceAI=(.*?)$')
        spliceai_res_is['SpliceAI Biology'] = spliceai_res_is['SpliceAI Biology'].str.split(',').str[0]
        spliceai_res_is['SpliceAI Biology'] = spliceai_res_is['SpliceAI Biology'].str.split('|').str[2:10]
        spliceai_res_is['SpliceAI Biology'] = spliceai_res_is.apply(lambda x: biology(x['SpliceAI Biology']), axis=1)
    spliceai_res_final = spliceai_res_is.append(spliceai_res_na, sort=False)
    return spliceai_res_final


def spliceai_max_score(spliceai_res):
    spliceai_res_na = spliceai_res[(spliceai_res['SpliceAI'] == '.') | (spliceai_res['SpliceAI'].isna())].copy()
    spliceai_res_is = spliceai_res[(spliceai_res['SpliceAI'] != '.') & (~spliceai_res['SpliceAI'].isna())].copy()
    spliceai_res_na['SpliceAI Pred'] = '.'
    if not spliceai_res_is.empty:
        spliceai_res_is['SpliceAI Pred'] = spliceai_res_is['SpliceAI'].str.extract('SpliceAI=(.*?)$')
        spliceai_res_is['SpliceAI Pred'] = spliceai_res_is['SpliceAI Pred'].str.split(',').str[0]
        spliceai_res_is['SpliceAI Pred'] = spliceai_res_is['SpliceAI Pred'].str.split('|').str[2:6]
        spliceai_res_is['SpliceAI Pred'] = spliceai_res_is.apply(lambda x: max(x['SpliceAI Pred']), axis=1)
    spliceai_res_final = spliceai_res_is.append(spliceai_res_na, sort=False)
    return spliceai_res_final


def pub_af(df, paf):
    cols = ['GnomAD EAS AF', 'GnomAD AF', '1000G AF', 'ESP6500 AF', 'ExAC EAS AF', 'ExAC AF']
    flt = df[df.apply(lambda x: (x[cols[0]] == '.' or float(x[cols[0]]) > paf) or
                                (x[cols[1]] == '.' or float(x[cols[1]]) > paf) or
                                (x[cols[2]] == '.' or float(x[cols[2]]) > paf) or
                                (x[cols[3]] == '.' or float(x[cols[3]]) > paf) or
                                (x[cols[4]] == '.' or float(x[cols[4]]) > paf) or
                                (x[cols[5]] == '.' or float(x[cols[5]]) > paf), axis=1)].copy()
    return flt


def splice_filter(strings, range_len):
    splice_range = int(min(re.findall(r'[\+-](\d+)[_a-zA-Z]', strings)))
    if splice_range <= range_len:
        return 'in splice ' + str(range_len)
    else:
        return 'out splice ' + str(range_len)


a = pd.read_csv('clinvar20190219plp-hgmd201901dm-intron-afless0.05.txt', sep='\t', low_memory=False)
a['splice range'] = a.apply(lambda x: splice_filter(x['cHGVS'], 20), axis=1)
a.to_csv('clinvar20190219plp-hgmd201901dm-intron-afless0.05.splice.txt', sep='\t', index=False)

# a_df = pd.read_csv('E:\\fangzhonghai\\work\\org\\BGI-WORK\\fangzhonghai\\BGI-Work\\Anno\\clinvar.2019-02-19.vcf', sep='\t', skiprows=27)
# a_df1, a_df2 = vcf_format_2_bgi_anno(a_df)
# a_df1.to_csv('test.trans.txt', sep='\t', index=False)

spliceai_df = pd.read_csv('E:\\fangzhonghai\\work\\fzh\\NA12878-2.chr1.spliceai.vcf', sep='\t', skiprows=range(28))
cols = ['#CHROM', 'POS', 'REF', 'ALT', 'SpliceAI']
spliceai_df.rename(columns={'INFO': 'SpliceAI'}, inplace=True)
spliceai_result = spliceai_df[cols].copy()
spliceai_result1 = spliceai_max_score(spliceai_result)
# spliceai_result1 = pd.read_csv('look.txt', sep='\t')
spliceai_result2 = spliceai_interpre(spliceai_result1)
spliceai_result3 = spliceai_biology(spliceai_result2)
spliceai_result3.to_excel('spliceai.test.xlsx', index=False)
