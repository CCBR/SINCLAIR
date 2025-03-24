import scvelo as scv
import sys

print()
adata = scv.read(str("velocyto/" + sys.argv[1]), cache=True)

scv.pp.remove_duplicate_cells(adata)
scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata)
scv.tl.recover_dynamics(adata)
# adata.write(sys.argv[2])
scv.tl.velocity(adata, mode="dynamical", n_jobs=16)
scv.tl.velocity_graph(adata)
del adata.raw
adata.__dict__["_var"] = adata.__dict__["_var"].rename(columns={"_index": "features"})
adata.write(str("velocyto/velo_" + sys.argv[1]))
