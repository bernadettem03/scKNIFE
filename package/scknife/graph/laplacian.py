import numpy as np
import networkx as nx
import scipy.sparse as sp
from scipy.sparse import diags


def build_graph_laplacian(G):
    A = nx.adjacency_matrix(G)
    D = diags(np.ravel(A.sum(axis=1)))
    L = D - A
    return A, D, L
