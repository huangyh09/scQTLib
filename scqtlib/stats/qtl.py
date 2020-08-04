import limix
import numpy as np
from statsmodels.sandbox.stats.multicomp import multipletests

from .glm import glm_LRT

def qtl_scan(y, GT, depth, M=None, interact=None, method='NB', **kwargs):
    n_gene = y.shape[1]
    n_sample = y.shape[0]
    
    pval_GT = np.array([None] * y.shape[1], float)
    fdr_GT  = np.array([None] * y.shape[1], float)
    pval_inter = np.array([None] * y.shape[1], float)
    fdr_inter  = np.array([None] * y.shape[1], float)
    
    for i in range(y.shape[1]): # SNP-gene pair
        if np.max(np.unique(GT[:, i], return_counts=True)[1]) > (0.95 * n_sample):            
            continue
            
        if method == 'NB':
            pval_GT[i], pval_inter[i] = qtl_NB(
                y[:, i], GT[:, i], depth, M, interact, **kwargs)[:2]
        else:
            pval_GT[i], pval_inter[i] = qtl_limix(
                y[:, i], GT[:, i], depth, M, interact, **kwargs)[:2]

    _idx = pval_GT >= 0
    print("%d out %d tests valid for FDR correction" %(sum(_idx), len(_idx)))
    if sum(_idx) > 1:
        fdr_GT[_idx] = multipletests(pval_GT[_idx], method='fdr_bh')[1]
        fdr_inter[_idx] = multipletests(pval_inter[_idx], method='fdr_bh')[1]

    return pval_GT, fdr_GT, pval_inter, fdr_inter


def qtl_NB(y, GT, depth, M=None, interact=None, **kwargs):
    """
    """
    if M is not None:
        M1 = np.append(np.log(depth + 1).reshape(-1, 1), M, axis=1)
    else:
        M1 = np.log(depth + 1).reshape(-1, 1)
    
    # if M is not None:
    #     M1 = np.append(depth.reshape(-1, 1), M, axis=1)
    # else:
    #     M1 = depth.reshape(-1, 1)

    try:
        nb_res1 = glm_LRT(y, GT, M=M1, **kwargs)
        pval_GT = nb_res1.p_LRT
    except:
        pval_GT, nb_res1 = None, None
        
    # Test for interaction
    pval_inter, nb_res2 = None, None
    if interact is not None:
        M2 = np.append(GT.reshape(-1, 1), M1, axis=1)
        try:
            nb_res2 = glm_LRT(y, GT * interact, M=M2, **kwargs)
            pval_inter = nb_res2.p_LRT
        except:
            pass
    
    return pval_GT, pval_inter, nb_res1, nb_res2


def qtl_limix(y, GT, depth=None, M=None, interact=None, add_intercept=True,
              family='normal'):
    """
    family: normal, poisson, binomial, etc.
    """
    Ken = None if depth is None else np.diag(1.0 / depth)
    if add_intercept:
        if M is not None:
            M1 = np.append(M, np.ones((y.shape[0], 1)), axis=1)
        else:
            M1 = np.ones((y.shape[0], 1))
    else:
        M1 = M.copy()
    
    qtl1 = limix.qtl.scan(GT.reshape(-1, 1), y.reshape(-1, 1), 
                          family, M=M1, K=Ken, verbose=False)
    pval_GT = np.array(qtl1.stats)[0, 4]
        
    pval_inter, qtl2 = None, None
    if interact is not None:
        M2 = np.append(GT.reshape(-1, 1), M1, axis=1)
        qtl2 = limix.qtl.scan((GT * interact).reshape(-1,1), y.reshape(-1, 1), 
                              family, M=M2, K=Ken, verbose=False)
        pval_inter = np.array(qtl2.stats)[0, 4]
    
    return pval_GT, pval_inter, qtl1, qtl2


    
# def qtl_scan_limix(y, GT, depth, M=None, interact=None):
#     """
#     """
#     n_gene = y.shape[1]
#     n_sample = y.shape[0]
    
#     pval_GT = np.array([None] * y.shape[1], float)
#     fdr_GT  = np.array([None] * y.shape[1], float)
#     pval_inter = np.array([None] * y.shape[1], float)
#     fdr_inter  = np.array([None] * y.shape[1], float)
    
#     Ken = np.diag(1.0 / depth)
    
#     for i in range(y.shape[1]): # SNP-gene pair
#         if np.max(np.unique(GT[:, i], return_counts=True)[1]) > (0.95 * n_sample):            
#             continue
            
        
            
#         M1 = M.copy()
#         M2 = np.append(GT[:, i:i+1], M1, axis=1)

#         qtl = limix.qtl.scan(GT[:, i:i+1], y[:, i+1], 
#                              'normal', M=M1, K=Ken, verbose=False)
#         pval_GT[i] = np.array(qtl.stats)[:,4]

#         qtl = limix.qtl.scan((GG[:, i] * disease).reshape(-1,1), y[:, i+1], 
#                              'normal', M=M2, K=Ken, verbose=False)
#         pval_inter[i] = np.array(qtl.stats)[:,4]

#     _idx = pval_GT >= 0
#     print(sum(_idx), len(_idx), pval_GT.shape)
#     if sum(_idx) > 1:
#         fdr_GT[_idx] = multipletests(pval_GT[_idx], method='fdr_bh')[1]
#         fdr_inter[_idx] = multipletests(pval_inter[_idx], method='fdr_bh')[1]

#     return pval_GT, fdr_GT, pval_inter, fdr_inter
    

# def qtl_scan_NB(y, GT, depth, M=None, interact=None):
#     """
#     """
#     n_gene = y.shape[1]
#     n_sample = y.shape[0]
    
#     pval_GT = np.array([None] * y.shape[1], float)
#     fdr_GT  = np.array([None] * y.shape[1], float)
#     pval_inter = np.array([None] * y.shape[1], float)
#     fdr_inter  = np.array([None] * y.shape[1], float)
    
#     for i in range(y.shape[1]): # SNP-gene pair
#         if np.max(np.unique(GT[:, i], return_counts=True)[1]) > (0.95 * n_sample):            
#             continue
            
#         if M is not None:
#             M1 = np.append(np.log(depth + 1).reshape(-1, 1), M, axis=1)
#         else:
#             M1 = np.log(depth + 1).reshape(-1, 1)
#         M2 = np.append(GT[:, i:i+1], M1, axis=1)

#         try:
#             nb_res1 = glm_LRT(y[:, i], GT[:, i], M=M1)
#             nb_res2 = glm_LRT(y[:, i], GT[:, i] * interact, M=M2)
#         except:
#             continue

#         pval_GT[i] = nb_res1.p_LRT
#         pval_inter[i] = nb_res2.p_LRT

#     _idx = pval_GT >= 0
#     print(sum(_idx), len(_idx), pval_GT.shape)
#     if sum(_idx) > 1:
#         fdr_GT[_idx] = multipletests(pval_GT[_idx], method='fdr_bh')[1]
#         fdr_inter[_idx] = multipletests(pval_inter[_idx], method='fdr_bh')[1]

#     return pval_GT, fdr_GT, pval_inter, fdr_inter

