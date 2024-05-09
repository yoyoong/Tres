import argparse
import os
import pandas as pd
import sys
import warnings
from tqdm.autonotebook import tqdm
import CytoSig

from Util import read_expression

warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser()
parser.add_argument('-E', "--expression_path", type=str, required=False, help="Gene expression file.",
                    default='/sibcb2/bioinformatics2/yangyue/project/240125_Tres/senecence/SPEED2_signaling.tsv')
parser.add_argument('-M', "--model_matrix_file", type=str, required=False, help="Quantitative signatures for cytokines.",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/0.model_file/signature.centroid.expand')
parser.add_argument('-D', "--output_file_directory", type=str, required=False, help="Directory for output files.",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/1.neutrophil_data/Gao2024')
parser.add_argument('-O', "--output_tag", type=str, required=False, help="Prefix for output files.", default='Gao2024.signaling')
args = parser.parse_args()

expression_path = args.expression_path
model_matrix_file = args.model_matrix_file
output_file_directory = args.output_file_directory
output_tag = args.output_tag

print("Process start!")
def compute_signaling(expression, model_matrix_file):
	# read model matrix file
    if not os.path.exists(model_matrix_file):
        sys.stderr.write('Cannot find model matrix file %s.\n' % model_matrix_file)
        sys.exit(1)
    else:
        signature = pd.read_csv(model_matrix_file, sep='\t', index_col=0)

    # compute signaling activity
    try:
        result_signaling = CytoSig.ridge_significance_test(signature, expression, alpha=1E4, nrand=1000, verbose=1)
    except ArithmeticError:
        sys.stderr.write('CytoSig regression failed.\n')
        sys.exit(1)
    
    # get the z-scores
    result_signaling = result_signaling[2]
    return result_signaling

result = pd.DataFrame()
if not os.path.isdir(expression_path):
    expression = read_expression(expression_path)
    result = compute_signaling(expression, model_matrix_file)
else:
    expression_list = sorted(os.listdir(expression_path))
    for expression_file in tqdm(expression_list, desc="Calculate"):
        expression = read_expression(os.path.join(expression_path, expression_file))
        sub_result = compute_signaling(expression, model_matrix_file)
        result = pd.concat([result, sub_result], axis=1)
        print(f"{expression_file} calculate end.")

signaling_filename = os.path.join(output_file_directory, f'{output_tag}.csv')
result.to_csv(signaling_filename, sep='\t')
print("Process end!")
