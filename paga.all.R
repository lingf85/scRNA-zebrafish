
### paga analysis for all the cell type
import numpy as np
import pandas as pd
import scanpy as sc
import os

import matplotlib.pyplot as plt
plt.switch_backend('agg')
import seaborn as sns
sns.set_style( "white" )
plt.rcParams[ "font.size" ] = 4.0
plt.rcParams[ "figure.dpi" ] = 300
plt.rcParams['axes.linewidth'] = 1
plt.rcParams['axes.labelsize'] = 6 
plt.rcParams[ "figure.figsize" ] = ( 4*0.8,4*0.8 )
plt.rcParams[ "font.serif" ] = 'Arial'


sc.settings.verbosity = 3 
sc.logging.print_versions()

data = sc.read_h5ad("pbmc.h5ad")  ### reading data 
data.obs['celltype2']=data.obs['seurat_clusters'].astype('str')+'_'+data.obs['celltype'].astype('str')
sc.pp.neighbors(data)
ident='celltype2'###…Ë÷√
data.obs[ident]=data.obs[ident].astype('category')
sc.tl.paga(data, groups=ident,model='v1.2')


#sc.settings.set_figure_params(dpi=100,dpi_save=1000,frameon=True)

data<-sc.tl.diffmap(data) ###
sc.pp.neighbors(data, n_neighbors=10, n_pcs=20,use_rep='X_diffmap')
sc.tl.paga(data, groups=ident,model='v1.2')

adata3=data

sc.pl.paga(adata3,color=[ident],save='.all.0.08.pdf',threshold=0.08,node_size_scale=0.5, node_size_power=1.0, edge_width_scale=2.0)
sc.pl.paga(adata3,color=[ident],save='.all.0.1.pdf',threshold=0.1,node_size_scale=0.5, node_size_power=1.0, edge_width_scale=2.0)

sc.tl.draw_graph(adata3,init_pos='paga' ,layout ='fa')
sc.pl.draw_graph(adata3, color=[ident], legend_loc='on data', title='',legend_fontsize=6,size=4,frameon =True,save='all.pdf')

adata3.uns['iroot'] = np.flatnonzero(adata3.obs[ident]  == '6_pMN2')[0]
sc.tl.dpt(adata3)
#sc.tl.draw_graph(adata3,init_pos='paga' ,layout ='fa')
sc.pl.draw_graph(adata3, color=[ident,'dpt_pseudotime'], legend_loc='on data', title='',legend_fontsize=6,size=4,vmax = 0.35,color_map ='Blues',save='.all.pseudo.pdf')

