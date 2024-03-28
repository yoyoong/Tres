import argparse
import time
import resource
import os, sys, pandas, numpy, pathlib

parser = argparse.ArgumentParser()
parser.add_argument('-E', "--expression_file", type=str, default=None, required=False, help="Gene expression file.")
parser.add_argument('-R', "--reaction_file", type=str, default=None, required=True, help="A tsv file with metabolite and gene association.")
parser.add_argument('-S', "--scoring_function", type=str, default='mean', required=False, help="Scoring function, one of [mean, median, max].")
parser.add_argument('-O', "--output_tag", type=str, default=None, required=False, help="Prefix for output files.")
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

def read_reaction_table(reaction_file):
	metmap = {}
	fin = open(reaction_file)
	for l in fin:
		fields = l.strip().split('\t')
		if len(fields) < 5:
			continue
		group = fields[4]
		if group not in ['product', 'substrate']:
			continue
		genes = set()
		chunks = fields[3].split(';')
		for block in chunks:
			genes.add(block.split('[')[0].strip())
		hmdb = fields[1]
		if hmdb not in metmap:
			metmap[hmdb]= {group: genes}
		elif group not in metmap[hmdb]:
			metmap[hmdb][group] = genes
		else:
			metmap[hmdb][group].update(genes)
	fin.close()
	return metmap

def estimate_metabolite(metmap, expression, method = 'mean'):
	for hmdb in metmap:
		if 'product' not in metmap[hmdb]:
			continue
		pGenes = metmap[hmdb]['product']
		pIndex = expression.index.isin(pGenes)
		if not any(pIndex):
			continue
		if method == 'mean':
			score = expression.loc[pIndex, :].mean()
		elif method == 'median':
			score = expression.loc[pIndex, :].median()
		elif method == 'max':
			score = expression.loc[pIndex, :].max()
		else:
			raise KeyError('Methods shoud be one of mean, median and max.')
		if 'substrate' in metmap[hmdb]:
			nGenes = metmap[hmdb]['substrate']
			nIndex = expression.index.isin(nGenes)
			if any(nIndex):
				if method == 'mean':
					score -= expression.loc[nIndex, :].mean()
				elif method == 'median':
					score -= expression.loc[nIndex, :].median()
				elif method == 'max':
					score -= expression.loc[nIndex, :].max()
				else:
					raise KeyError('Methods shoud be one of mean, median and max.')
		# 
