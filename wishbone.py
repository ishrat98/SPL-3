
#!/usr/local/bin/python

import dynclipy
task = dynclipy.main()
# task = dynclipy.main(
#   ["--dataset", "/code/example.h5", "--output", "/mnt/output"],
#   "/code/definition.yml"
# )

import wishbone
import os
import sys
import json
import pandas as pd

import time
checkpoints = {}


#   ____________________________________________________________________________
#   Load data                                                               ####
task["counts"].to_csv("/counts.csv")

p = task["parameters"]

# get start cell(s)
start_cell = task["priors"]["start_id"]
if isinstance(start_cell, list):
  start_cell = np.random.choice(start_cell)

# get markers if given
if "features_id" in task["priors"]:
  markers = task["priors"]["features_id"]
else:
  markers = "~"

# get number of end states if given
if "end_n" in task["priors"]:
  branch = task["priors"]["end_n"] > 1
else:
  branch = True

checkpoints["method_afterpreproc"] = time.time()

#   ____________________________________________________________________________
#   Create trajectory                                                       ####
# normalise data
scdata = wishbone.wb.SCData.from_csv("/counts.csv", data_type='sc-seq', normalize=p["normalise"])
scdata.run_pca()
scdata.run_diffusion_map(knn=p["knn"], epsilon=p["epsilon"], n_diffusion_components=p["n_diffusion_components"], n_pca_components=p["n_pca_components"], markers=markers)

# check waypoints parameter to be lower than # cells
if p["num_waypoints"] > scdata.data.shape[0]:
  p["num_waypoints"] = scdata.data.shape[0]

# run wishbone
wb = wishbone.wb.Wishbone(scdata)
wb.run_wishbone(start_cell=start_cell, components_list=list(range(p["n_diffusion_components"])), num_waypoints=int(p["num_waypoints"]), branch=branch, k=p["k"])

checkpoints["method_aftermethod"] = time.time()

#   ____________________________________________________________________________
#   Process output & save                                                   ####
# progressions
progressions = wb.trajectory.reset_index()
progressions.columns = ["cell_id", "percentage"]
if branch:
  progressions["from"] = pd.Series(["M1", "M2", "M2"])[wb.branch - 1].tolist()
progressions["to"] = pd.Series(["M2", "M3", "M4"])[wb.branch - 1].tolist()
else:
  progressions["from"] = "M1"
progressions["to"] = "M2"

# milestone network
if branch:
  milestone_network = pd.DataFrame({"from":["M1", "M2", "M2"], "to":["M2", "M3", "M4"], "length":1, "directed":True})
else:
  milestone_network = pd.DataFrame({"from":["M1"], "to":["M2"], "length":1, "directed":True})


# pseudotime
pseudotime = wb.trajectory.reset_index()
pseudotime.columns = ["cell_id", "pseudotime"]

# dimred
dimred = wb.scdata.diffusion_eigenvectors
dimred.index.name = "cell_id"
dimred = dimred.reset_index()

# cell ids
cell_ids = pd.DataFrame({
  "cell_ids": dimred.index
})

# save
dataset = dynclipy.wrap_data(cell_ids = dimred["cell_id"])
dataset.add_trajectory(
  milestone_network = milestone_network,
  progressions = progressions
)
dataset.add_dimred(dimred = dimred)
dataset.add_timings(timings = checkpoints)
dataset.write_output(task["output"])
