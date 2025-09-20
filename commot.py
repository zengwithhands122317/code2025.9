# pip install commot
import os
import gc
import ot
import pickle
import anndata
import scanpy as sc
import pandas as pd
import numpy as np
from scipy import sparse
from scipy.stats import spearmanr, pearsonr
from scipy.spatial import distance_matrix
import matplotlib.pyplot as plt
import commot as ct

## 读取数据
#ALT+SHIFT+E
adata= sc.read_visium(path="E:/article4.5/smrsp2/HRCC-1",count_file ='filtered_feature_bc_matrix.h5')
adata.var_names_make_unique()

sc.pp.normalize_total(adata, inplace=True)

sc.pp.log1p(adata)

## 读取通讯数据库
df_cellchat = ct.pp.ligand_receptor_database(species='human', signaling_type='Secreted Signaling', database='CellChat')
#df_cellchat = ct.pp.ligand_receptor_database(species='mouse', signaling_type='Secreted Signaling', database='CellChat')

df_cellchat_filtered = ct.pp.filter_lr_database(df_cellchat, adata, min_cell_pct=0.05)


## 可视化------------


## 每次做一个！！
LR=np.array([['SPP1', 'ITGAV_ITGB1', 'SPP1']],dtype=str)
df_ligrec = pd.DataFrame(data=LR)

ct.tl.spatial_communication(adata,
database_name='cellchat', df_ligrec=df_ligrec, dis_thr=500, heteromeric=True, pathway_sum=True)


adata.obsm['commot-CellChat-sum-sender']
#Plot the amount of sent and received signal, for example, of the LR pair Fgf1-Fgfr1.
pts = adata.obsm['spatial']
s = adata.obsm['commot-cellchat-sum-sender']['s-ITGAV_ITGB1-TNFRSF1A-SPP1']
r = adata.obsm['commot-cellchat-sum-receiver']['r-ITGAV_ITGB1-TNFRSF1A-SPP1']


import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt








fig, ax = plt.subplots(1,2, figsize=(10,4))
ax[0].scatter(pts[:,0], pts[:,1], c=s, s=5, cmap='Blues')
ax[0].set_title('Sender')
ax[1].scatter(pts[:,0], pts[:,1], c=r, s=5, cmap='Reds')
ax[1].set_title('Receiver')
plt.close()


## 可视化
ct.tl.communication_direction(adata, database_name='cellchat', lr_pair=('TNF', 'TNFRSF1A'), k=5)
ct.pl.plot_cell_communication(adata, database_name='cellchat', lr_pair=('TNF', 'TNFRSF1A'), plot_method='grid', background_legend=True,
    scale=0.000002, ndsize=8, grid_density=0.4, summary='sender', background='image', clustering='leiden', cmap='Alphabet',
    normalize_v = True, normalize_v_quantile=0.995)
plt.close()

