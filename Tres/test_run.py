# python geneset_score.py -E /sibcb2/bioinformatics/ImmuneDNB/Tres_data/sc_cohorts/Breast.GSE156728.10x.pickle.gz -G Tres.kegg -S signature_name_file.txt -O Breast.GSE156728.10x

import argparse
import time
import resource
import os, sys, pandas, numpy, pathlib
import CytoSig

geneset_file = 'Tres.kegg'
expression = read_expression('/sibcb2/bioinformatics/ImmuneDNB/Tres_data/sc_cohorts/Breast.GSE156728.10x.pickle.gz')
result_response = pandas.read_csv('Breast.GSE156728.10x_Prolifertion.tsv', sep='\t', index_col=0)
result_signaling = pandas.read_csv('Breast.GSE156728.10x_Signaling.tsv', sep='\t', index_col=0)

response_key = 'Proliferation'
signaling_key = 'TGFB1'

flag_group = ['.'.join(v.split('.')[:2]) for v in expression.columns]
expression_group = expression.groupby(flag_group, axis=1)
merge = []

for title, expression_sub in expression_group:
	print(title)

y = (result_response.loc[response_key]).loc[expression_sub.columns]
y = y.to_frame(name=None)
N = expression_sub.shape[1]

X = pandas.DataFrame(numpy.zeros((N,1)), columns = ['pivot'], index = expression_sub.columns)
X.loc[:,'pivot'] = result_signaling.loc[signaling_key, expression_sub.columns]

result = CytoSig.ridge_significance_test(X, y, alpha=10000, alternative="two-sided", nrand=1000, verbose=1)
