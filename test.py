import scanpy as sc
import numpy as np
import h5py
from scipy.sparse import coo_matrix

filename = '/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.qc_data/AEL_GSE142213_expression.h5'
h5file = h5py.File(filename, "r")
h5file.visit(print)

# adata = sc.read_10x_h5(filename)
adata = sc.read_hdf(filename)

matrix = h5file['matrix']
data = matrix['data']
shape = np.array(matrix['shape'])
indices = np.array(matrix['indices'])
indptr = np.array(matrix['indptr'])

row = []
for index in range(1, len(indptr)):
    row = row + [index - 1] * indptr[index]

data = coo_matrix((data[:], (indices, row)), shape=shape, dtype=float).toarray()

data_exp = np.exp(data)
col_sum = np.sum(data_exp, axis=0)
print(matrix)