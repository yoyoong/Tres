import os, sys
import pandas
import pickle
def read_expression(input_file):
    f = os.path.basename(input_file)
    if f.find('.pickle') >= 0:
        expression = pandas.read_pickle(input_file)
    else:
        expression = pandas.read_csv(input_file, sep=',', index_col=0)
    # read input
    # try:
    #     f = os.path.basename(input_file)
    #     if f.find('.pickle') >= 0:
    #         expression = pandas.read_pickle(input_file)
    #     else:
    #         expression = pandas.read_csv(input_file, sep=',', index_col=0)
    # except:
    #     sys.stderr.write('Fail to open input file %s\n' % input_file)
    #     sys.exit(1)

    # gene and sample names must be unique
    assert expression.index.value_counts().max() == 1
    assert expression.columns.value_counts().max() == 1

    # print('expression matrix dimension', expression.shape)
    return expression