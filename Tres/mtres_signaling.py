import argparse
import time
import resource
import os, sys, pandas, numpy, pathlib
import CytoSig

parser = argparse.ArgumentParser()
parser.add_argument('-E', "--expression_file", type=str, default=None, required=True, help="Gene expression file.")
parser.add_argument('-M', "--model_matrix_file", type=str, default=None, required=True, help="Quantitative signatures for cytokines.")
parser.add_argument('-O', "--output_tag", type=str, default=None, required=True, help="Prefix for output files.")
args = parser.parse_args()

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

def compute_signaling(expression, model_matrix_file):
	# read model matrix file
    if not os.path.exists(model_matrix_file):
        sys.stderr.write('Cannot find model matrix file %s.\n' % model_matrix_file)
        sys.exit(1)
    else:
        signature = pandas.read_csv(model_matrix_file, sep='\t', index_col=0)

    # compute signaling activity
    try:
        result_signaling = CytoSig.ridge_significance_test(signature, expression, alpha=1E4, nrand=1000, verbose=1)
    except ArithmeticError:
        sys.stderr.write('CytoSig regression failed.\n')
        sys.exit(1)
    
    # get the z-scores
    result_signaling = result_signaling[2]
    return result_signaling

expression = read_expression(args.expression_file)
result_signaling = compute_signaling(expression, args.model_matrix_file)
result_signaling.to_csv(args.output_tag + '_Signaling.tsv', sep='\t')
