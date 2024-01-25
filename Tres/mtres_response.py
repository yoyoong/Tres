import argparse
import time
import resource
import os, sys, pandas, numpy, pathlib
import CytoSig

parser = argparse.ArgumentParser()
parser.add_argument('-E', "--expression_file", type=str, default=None, required=True, help="Gene expression file.")
parser.add_argument('-G', "--genesets_GMT_file", type=str, default=None, required=True, help="Background gene sets in GMT format.")
parser.add_argument('-S', "--signature_name_file", type=str, default=None, required=True, help="Names of the signatures, one name in one line.")
parser.add_argument('-O', "--output_tag", type=str, default=None, required=True, help="Prefix for output files.")
args = parser.parse_args()

def profile_geneset_signature(expression, geneset_file, name_file):
    # all gene sets
    signature = []
    fin = open(geneset_file)
    for l in fin:
        fields = l.strip().split('\t')
        s = fields[2:]
        signature.append(pandas.Series(numpy.ones(len(s)), index=s, name=fields[0]))
    fin.close()

    # gene set names
    names = []
    fin = open(name_file)
    for l in fin:
        fields = l.strip().split('\t')
        names.append(fields[0])
    fin.close()

    # pandas data frame
    signature = pandas.concat(signature, axis=1, join='outer', sort=False)
    signature.fillna(0, inplace=True)
    
    common = expression.index.intersection(signature.index)
    signature, expression = signature.loc[common], expression.loc[common]
    
    background = signature.mean(axis=1)
    background.name = 'study bias'
    
    X = signature.loc[:, names].mean(axis=1)
    X.name = 'Proliferation'
    
    X = pandas.concat([background, X], axis=1, join='inner')
    
    # regression
    result = CytoSig.ridge_significance_test(X, expression, alpha=0, alternative="two-sided", nrand=0, verbose=1)
    
    return result[2]

def read_expression(input_file):
    # read input
    try:
        f = os.path.basename(input_file)
        if f.find('.pickle') >= 0:
            print(f)
            expression = pandas.read_pickle(input_file)
        else:
            expression = pandas.read_csv(input_file, sep='\t', index_col=0)
    except:
        sys.stderr.write('Fail to open input file %s\n' % input_file)
        sys.exit(1)
    
    # gene and sample names must be unique
    assert expression.index.value_counts().max() == 1
    assert expression.columns.value_counts().max() == 1
    
    print('input matrix dimension', expression.shape)
    return expression

expression = read_expression(args.expression_file)
result = profile_geneset_signature(expression, args.genesets_GMT_file, args.signature_name_file)

result.to_csv(args.output_tag + '_Prolifertion.tsv', sep='\t')
