import os, sys, pathlib, pandas


def strip_cancer_type_list(lst):
    return [v.split('_')[0].split('.')[0] if v.find('Liver') < 0 else v.split('_')[0].split('.', 1)[1] for v in lst]

def median_merge(output_path = '/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/test'):
    qthres = 0.05
    frac_thres = 1e-3

    output = os.path.join(output_path, 'merge')

    lst = []

    pivots = ['TGFB1', 'TRAIL', 'PGE2']

    for pivot in pivots:
        merge = []

        # extract CD8 T cells
        dataset_Tumor_CD8 = ['Post_T_CD8']
        for cell_pivot in dataset_Tumor_CD8:
            result = pandas.read_csv('/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/test/output.run', sep='\t', index_col=0)

            # focus on current pivot
            flag = [(v.split('.')[-1] == pivot) for v in result.columns]
            assert sum(flag) > 0

            result = result.loc[:, flag]
            result = result.loc[result.isnull().mean(axis=1) < 1]

            # extract t-values
            flag = [(v.find('t.' + cell_pivot) == 0) for v in result.columns]
            assert sum(flag) > 0

            result_t = result.loc[:, flag]

            # strip t and pivot
            result_t.columns = ['.'.join(v.split('.')[1:-1]) for v in result_t.columns]
            assert result_t.columns.value_counts().max() == 1  # if sample ID extraction is correct, must have no redundancy

            # extract q-values
            flag = [(v.find('q.' + cell_pivot) == 0) for v in result.columns]
            assert sum(flag) > 0

            result_q = result.loc[:, flag]

            # strip t and pivot
            result_q.columns = ['.'.join(v.split('.')[1:-1]) for v in result_q.columns]
            assert result_q.columns.value_counts().max() == 1  # if sample ID extraction is correct, must have no redundancy

            flag = (result_q < qthres).mean() > frac_thres

            if flag.sum() == 0:
                # nothing to include
                continue

            result = result_t.loc[:, flag]

            merge.append(result)

        merge = pandas.concat(merge, axis=1, join='outer', sort=False)
        assert merge.columns.value_counts().max() == 1

        merge.to_csv(output + '.' + pivot, sep='\t', index_label=False)
        lst.append(merge)

    # create an average score as Tres score
    common_col = None

    for result in lst:
        print(result.shape)
        # get single cells that are profiled by all signals
        if common_col is None:
            common_col = result.columns
        else:
            common_col = common_col.intersection(result.columns)

    print(common_col.shape)
    for i, result in enumerate(lst): lst[i] = result.loc[:, common_col]

    # create median signature across three immuno-suppressive signals
    merge = pandas.concat(lst, axis=1, join='inner')

    cnt_map = merge.columns.value_counts()
    assert sum(cnt_map != len(pivots)) == 0

    merge = merge.groupby(merge.columns, axis=1).median()
    print(merge.shape)

    merge.to_csv(output + '.Median', sep='\t', index_label=False)

    # create cancer-level merge
    # merge = pandas.read_csv(output + '.Median', sep='\t', index_col=0)
    dddd = merge.isnull().mean(axis=1)
    merge = merge.loc[merge.isnull().mean(axis=1) < 0.5]
    print(merge.shape)

    # don't merge liver cancer
    flag = strip_cancer_type_list(merge.columns)

    merge = merge.groupby(flag, axis=1).median()
    merge.to_csv(output, sep='\t', index_label=False)

    # create a median signature across all dataset for prediction purpose later
    merge = merge.dropna().median(axis=1)
    merge.name = 'Tres'
    merge.to_csv(output + '.signature', sep='\t', index_label=False)

if __name__ == '__main__': median_merge()