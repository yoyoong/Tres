from scipy.sparse import csr_matrix
from scipy.stats import norm
import numpy as np
import pandas as pd
import sys, os
from tqdm.autonotebook import tqdm
import torch
from torch import Tensor

def upperlower(data, boxsize):
    (n1, n2) = data.shape  # n1 gene; n2 sample
    upper = np.zeros((n1, n2), dtype=np.float32)
    lower = np.zeros((n1, n2), dtype=np.float32)

    for i in tqdm(range(0, n1), desc="Get upperlower"):
        s1 = sorted(data[i, :])
        s2 = data[i, :].argsort()
        sum = int(np.sum(np.sign(s1)))
        n3 = n2 - sum
        h = int(boxsize / 2 * sum + 0.5) # traditional rounding e.g. 2.4->2 2.5->3
        k = 0
        while k <= (n2 - 1):
            s = 0
            while (k + s + 1 <= (n2 - 1)) and (s1[k + s + 1] == s1[k]):
                s = s + 1
            if s >= h:
                upper[i, s2[k:k + s + 1]] = data[i, s2[k]]
                lower[i, s2[k:k + s + 1]] = data[i, s2[k]]
            else:
                upper[i, s2[k:k + s + 1]] = data[i, s2[min(n2 - 1, k + s + h)]]
                lower[i, s2[k:k + s + 1]] = data[i, s2[max(n3 * (n3 > h), k - h)]]

            k = k + s + 1
    return (upper, lower)

def getCSNMatrix(gem: Tensor, upper: Tensor, lower: Tensor,
                 index: int, # cell index
                 is_weight: float, # whether the network's edge is weighted
                 has_zero: bool, # whether the network's edge include the zero expression genes
                 alpha: float, # CSN alpha value
                 device): # calculate device, cpu or cuda
    (n1, n2) = gem.shape
    eps = torch.finfo(torch.float32).eps
    B = torch.zeros((n1, n2), dtype=torch.float32).to(device)
    for j in range(0, n2):
        upper_cutoff = (gem[:, j] < upper[:, index]) | (torch.isclose(gem[:, j], upper[:, index], rtol=1e-4))
        lower_cutoff = (gem[:, j] > lower[:, index]) | (torch.isclose(gem[:, j], lower[:, index], rtol=1e-4))
        if has_zero:
            B[:, j] = upper_cutoff * lower_cutoff
        else:
            B[:, j] = upper_cutoff * lower_cutoff * (gem[:, index] > 0)
    a = B.sum(axis=1)
    a = torch.reshape(a, (n1, 1))
    temp = (torch.mm(B, B.T) * n2 - torch.mm(a, a.T)) / torch.sqrt(torch.mm(a, a.T) * torch.mm((n2 - a), (n2 - a).T) / (n2 - 1) + eps)
    # temp[temp < 0] = 0
    temp.fill_diagonal_(0)
    matrix = csr_matrix(temp.cpu().detach().numpy()).tocoo()
    matrix = matrix.multiply(matrix > norm.ppf(1 - alpha))
    if not is_weight:
        matrix = (matrix > norm.ppf(1 - alpha))
    return matrix