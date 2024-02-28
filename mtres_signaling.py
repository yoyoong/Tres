import argparse
import os
import pandas
import sys
import warnings

import CytoSig

from Util import read_expression

warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser()
parser.add_argument('-E', "--expression_file", type=str, required=False, help="Gene expression file.",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/0.Tres_data/sc_cohorts/Nasopharyngeal.GSE162025.10x.pickle.gz')
parser.add_argument('-M', "--model_matrix_file", type=str, required=False, help="Quantitative signatures for cytokines.",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/1.model_data/signature.centroid.expand')
parser.add_argument('-D', "--output_file_directory", type=str, required=False, help="Directory for output files.",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.qc_data/Signaling')
parser.add_argument('-O', "--output_tag", type=str, required=False, help="Prefix for output files.", default='Nasopharyngeal.GSE162025.10x')
args = parser.parse_args()

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
signaling_filename = os.path.join(args.output_file_directory, f'{args.output_tag}.Signaling.csv')
result_signaling.to_csv(signaling_filename, sep='\t')
print("Process end!")
