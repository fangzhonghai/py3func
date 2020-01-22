# -*- coding:utf-8 -*-
import pandas as pd
import pyfaidx


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
    df.loc[df['MuType'] == 'delins', 'POS'] = df.loc[df['MuType'] == 'delins', 'Start'] + 1
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
    df.loc[(df['MuType'] == 'ins') & (df['REF'].str[0] != df['ALT'].str[0]), 'MuType'] = 'delins'
    df.loc[(df['MuType'] == 'del') & (df['REF'].str[0] != df['ALT'].str[0]), 'MuType'] = 'delins'
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
    a = df[['#CHROM', 'POS', 'REF', 'ALT', '#Chr', 'Start', 'Stop', 'Ref', 'Call']].copy()
    df.drop(columns=['MuType'], inplace=True)
    return a, df


a_df = pd.read_csv('E:\\fangzhonghai\\work\\org\\BGI-WORK\\fangzhonghai\\BGI-Work\\Anno\\clinvar.2019-02-19.vcf', sep='\t', skiprows=27)
a_df1, a_df2 = vcf_format_2_bgi_anno(a_df)
a_df1.to_csv('test.trans.txt', sep='\t', index=False)
