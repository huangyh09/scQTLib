import numpy as np
import matplotlib.pyplot as plt

vega_20_scanpy = [
    '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', 
    '#e377c2',  # '#7f7f7f' removed grey
    '#bcbd22', '#17becf', '#aec7e8', '#ffbb78', '#98df8a', '#ff9896',
    '#c5b0d5', '#c49c94', '#f7b6d2',  # '#c7c7c7' removed grey
    '#dbdb8d', '#9edae5', '#ad494a', '#8c6d31']

wzhu_colors = ["#4796d7", "#f79e54", "#79a702", "#df5858", "#556cab", 
               "#de7a1f", "#ffda5c", "#4b595c", "#6ab186", "#bddbcf", 
               "#daad58", "#488a99", "#f79b78", "#ffba00"]

wzhu_gradient = ["#2e496e", "#375e97", "#556cab", "#2987ba", "#488a99",
                 "#20948c", "#5f958f", "#6ab186", "#bddbcf"]

def scatter_adata(adata, mode='tsne', key=None, size=0.5, 
                  id_list=None, show_num=True, color_list=None, 
                  figsize=None, label_on=False, label_size=10):
    """
    Plot adata
    ==========
    
    Example
    -------
    fig = plt.figure(figsize=[5, 4], dpi=300)
    scatter_adata(adata[adata.obs['donor_ids0']=="A_singlet"], 
                  key='vdj_type', size=0.2, 
                  id_list=["B cell", "T cell", "BCR_TCR", "Others"],
                  color_list=['#1f77b4', '#d62728', '#333333', '#c7c7c7'])
    plt.legend(markerscale=10, prop={'size': 8})
    plt.grid(False)
    plt.xlabel("t-SNE 1")
    plt.ylabel("t-SNE 2")
    plt.xticks([])
    plt.yticks([])
    plt.show()
    """
    X = adata.obsm['X_%s' %(mode.lower())]
    ids = adata.obs[key]
    if id_list is None:
        id_list = np.unique(ids)
    for ii in range(len(id_list)):
        idx = np.where(ids == id_list[ii])[0]
        if color_list is not None and len(color_list) > ii:
            color_use = color_list[ii]
        else:
            color_use = wzhu_colors[ii % len(wzhu_colors)]
        if show_num:
            label_use = id_list[ii] +": %d" %len(idx)
        else:
            label_use = id_list[ii]
        plt.scatter(X[idx, 0], X[idx, 1], s=size, label=label_use, 
                    color=color_use)
        if label_on:
            plt.text(np.median(X[idx, 0]), np.median(X[idx, 1]), 
                     label_use, horizontalalignment='center', 
                     verticalalignment='center', fontsize=label_size)


            
def Gboxplot(G, x, gname=None, SNPname=None):
    """Boxplot for quantitive traits vs genotype
    TODO: add documentation!
    """
    import seaborn as sns
    
    idx = (G >= 0) * (x == x)
    sns.boxplot(G[idx], x[idx], showfliers=False)
    sns.swarmplot(G[idx], x[idx], color=".25")
    
    if gname is not None:
        plt.ylabel("log10(CPM + 1), %s" %gname)
    else:
        plt.ylabel("log10(CPM + 1)")
        
    if SNPname is not None:
        SNP_val = SNPname.split("_")
        if len(SNP_val) == 4:
            _GT = [SNP_val[2] + SNP_val[2] + ": %d" %(sum(G == 0)),
                   SNP_val[2] + SNP_val[3] + ": %d" %(sum(G == 1)),
                   SNP_val[3] + SNP_val[3] + ": %d" %(sum(G == 2))]
            plt.xticks([0, 1, 2], _GT)
        plt.xlabel(SNPname)
    plt.xlim(-0.5, 2.5)
    
