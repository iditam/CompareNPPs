import os
import os.path
import scipy.stats as st
import scipy.spatial as sp
import pandas as pd
import yaml as yml
import math

# # construct the argument parse and parse the arguments
# ap = argparse.ArgumentParser()
# ap.add_argument("-npp", "--npp", required=True,
#                 help="npp ID")
# ap.add_argument("-s", "--start", required=True,
#                 help="first row to include")
# ap.add_argument("-e", "--end", required=True,
#                 help="last row (not included)")
# ap.add_argument("-c", "--cor", default="pearson",
#                 help="correlation method (pearson, spearman, kendall). default=pearson")
#
# args = vars(ap.parse_args())

npp_id = snakemake.wildcards.NPP
cor_method = snakemake.wildcards.CORRELATION
index_range = snakemake.params.current_range
output_file = snakemake.output.filename

with open("config.yaml", 'r') as stream:
    config_yaml = yml.load(stream)

npp_file = os.path.join(config_yaml["dataFolder"], npp_id,
                        config_yaml["nppFilePrefix"] + "_" + npp_id)
f = open(npp_file, 'r')
NPP = pd.read_table(f, delimiter="\t", index_col=0)
# print("NPP loaded")
f.close()

n = NPP.shape[0]

NPP_genes = NPP.index
NPP = NPP.values

# "METAP2" in NPP_genes
# NPP = NPP[:n, :] TBD: IS IT NEEDED?

if cor_method == "pearson":
    cor_func = st.pearsonr
elif cor_method == "spearman":
    cor_func = st.spearmanr
elif cor_method == "kendall":
    cor_func = st.kendalltau
elif cor_method == "manhattan":
    cor_func = sp.distance.cityblock
else:
    exit('Unknown correlation method: ' + cor_method)


# input:
# ix = index in upper triangle vector
# n = length of upper triangle vector
# output:  matrix indices
# Tuple [ row_index, col_index ]
def upper_triangle_rowcol(ix, n):
    row_index = math.floor(n-(1+math.sqrt(1+4*(n**2-n-2*ix)))//2)
    col_index = n-(2*n-row_index+1)*row_index//2+ix+row_index
    return([int(row_index-1), int(col_index-1)])  # adjusting to python indices starting at 0


lines = ""
for i in range(index_range[0], index_range[1]+1):
    row_index, col_index = upper_triangle_rowcol(i, n)
    if (cor_method == "manhattan"):
        cor_value = -round(cor_func(NPP[row_index, :], NPP[col_index, :]), 3)
    else:
        cor_value = round(cor_func(NPP[row_index, :], NPP[col_index, :])[0], 3)
    if (NPP_genes[row_index] < NPP_genes[col_index]):
        genes = NPP_genes[row_index] + ":" + NPP_genes[col_index]
    else:
        genes = NPP_genes[col_index] + ":" + NPP_genes[row_index]
    if (cor_method == "manhattan"):
        lines += "{},{}\n".format(cor_value, genes)
    else:
        lines += "{:.3},{}\n".format(cor_value, genes)


with open(output_file, "w+") as f:
    f.write(lines)
