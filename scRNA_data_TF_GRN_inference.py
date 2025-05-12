import os
import glob
import pickle
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
print(np.__version__) # must be 1.21.0

from dask.diagnostics import ProgressBar
from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2
from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import modules_from_adjacencies, load_motifs
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell
from matplotlib import cm
from matplotlib.colors import Normalize, to_hex

import seaborn as sns

# designate data location
DATA_FOLDER = "/home/hjlee/monocle3,pyscenic,palantir,iqcell/2025.05.belzutifan/DATA_FOLDER"
RESOURCES_FOLDER = "/home/hjlee/monocle3,pyscenic,palantir,iqcell/RESOURCES_FOLDER"
DATABASE_FOLDER = "/home/hjlee/monocle3,pyscenic,palantir,iqcell/DATABASE_FOLDER"
EXPRESSIONDATA_FOLDER = "/home/hjlee/monocle3,pyscenic,palantir,iqcell/EXPRESSIONDATA_FOLDER"

DATABASES_GLOB = os.path.join(DATABASE_FOLDER, "hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather")
MOTIF_ANNOTATIONS_FNAME = os.path.join(RESOURCES_FOLDER, "motifs-v9-nr.hgnc-m0.001-o0.0.tbl")
MM_TFS_FNAME = os.path.join(RESOURCES_FOLDER, 'hs_hgnc_tfs.txt')
SC_EXP_FNAME = os.path.join(EXPRESSIONDATA_FOLDER, "20250508_data_count.csv")
REGULONS_FNAME = os.path.join(DATA_FOLDER, "20250508_regulons_all_gene.p") # for save
MOTIFS_FNAME = os.path.join(DATA_FOLDER, "20250508_motifs_all_gene.csv") # for save

ex_matrix = pd.read_csv(SC_EXP_FNAME, header=0, index_col=0).T # for making cell (row) and gene (column)
ex_matrix.shape 
tf_names = load_tf_names(MM_TFS_FNAME)

db_fnames = glob.glob(DATABASES_GLOB)
def name(fname):
    return os.path.splitext(os.path.basename(fname))[0]
dbs = [RankingDatabase(fname=fname, name=name(fname)) for fname in db_fnames]
dbs

# STEP1 : co-expression-based TF-target interference using GRNBoost2 algorithm 
# Regulons are derived from adjacencies using arboreto package
adjacencies = grnboost2(ex_matrix, tf_names=tf_names, verbose=True)

# STEP2 : TF is added to the module & modules that have less than 20 genes are removed
modules = list(modules_from_adjacencies(adjacencies, ex_matrix))

# Calculate a list of enriched motifs and the corresponding target genes for all modules.
with ProgressBar():
    df = prune2df(dbs, modules, MOTIF_ANNOTATIONS_FNAME)

# Create regulons from this table of enriched motifs
regulons = df2regulons(df)

# Save the enriched motifs and the discovered regulons to disk.
df.to_csv(MOTIFS_FNAME)
with open(REGULONS_FNAME, "wb") as f:
    pickle.dump(regulons, f)

# The clusters can be leveraged via the dask framework:
df = prune2df(dbs, modules, MOTIF_ANNOTATIONS_FNAME)

# STEP3 : Reloading the enriched motifs and regulons from file  
df = load_motifs(MOTIFS_FNAME)
with open(REGULONS_FNAME, "rb") as f:
    regulons = pickle.load(f)

auc_mtx = aucell(ex_matrix, regulons)
print(auc_mtx)

# We characterize the different cells in a single-cell transcriptomics experiment via the enrichment of the previously discovered regulons. 
# Enrichment of a regulon is measured as the Area Under the recovery Curve (AUC) of the genes that define this regulon.
auc_mtx = aucell(ex_matrix, regulons)
auc_df = auc_mtx.copy()
auc_df.index.name = 'cell_id'
auc_df.reset_index(inplace=True)

umap_df = pd.read_csv("/home/hjlee/monocle3,pyscenic,palantir,iqcell/20250508_order.csv")
print(umap_df.head())

umap_df_sub = pd.merge(umap_df, auc_df, on='cell_id', how='inner')
umap_df_sub.to_csv('/home/hjlee/monocle3,pyscenic,palantir,iqcell/20250508_auc_matrix_meta_all_gene.csv')
umap_df_sub.set_index("cell_id", inplace=True)

unique_btypes = sorted(umap_df_sub["cell_type"].dropna().unique())
palette = sns.color_palette("tab10", n_colors=len(unique_btypes))
lut = dict(zip(unique_btypes, palette))
row_colors = umap_df_sub.reindex(auc_mtx.index)["cell_type"].map(lut)

auc_mtx = auc_mtx.loc[umap_df_sub.index]
# heatmap
selected_genes = ['NR2F1(+)','MAF(+)','PBX1(+)', 'ESRRA(+)',
                  'ETV7(+)','RELB(+)','STAT1(+)','MYC(+)','MAFF(+)','EGR1(+)','ATF3(+)','KLF6(+)','CEBPD(+)','JUN(+)','JUNB(+)','FOS(+)','CEBPB(+)','FOSB(+)']

auc_mtx_filtered = auc_mtx[selected_genes]
sns.clustermap(
    auc_mtx_filtered,
    figsize=(8, 5),
    row_colors=row_colors,
    z_score=1,
    vmin=-2,
    vmax=2,
    cmap="vlag",
    row_cluster=False,
    col_cluster=False,
    yticklabels=False
)
auc_mtx.to_csv("/home/hjlee/monocle3,pyscenic,palantir,iqcell/20250508_auc_matrix_all_gene.csv")
plt.savefig('/home/hjlee/monocle3,pyscenic,palantir,iqcell/20250508_AUCell_all_gene.png')
plt.show()
