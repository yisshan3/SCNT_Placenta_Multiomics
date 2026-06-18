"""
==============================================================================
Description: Compositional analysis on snRNA-seq data using scCODA to identify statistically credible changes in cell type proportions between SCNT and WT groups.
==============================================================================
"""


# Setup and Import Packages ----------------------------------------------------
import warnings
warnings.filterwarnings("ignore")

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sccoda.util import comp_ana as ana
from sccoda.util import comp_ana as mod
from sccoda.util import cell_composition_data as dat
from sccoda.util import data_visualization as viz


# Load and Prepare Metadata ----------------------------------------------------
os.chdir('/home4/ssyi/Mouse_Placenta/codes/Review/3.1.scCODA')
meta = pd.read_csv("E145_Tropho_cell_metadata_for_sccoda.csv")

# Standardize 'stage' labels (125/145/185 -> E12.5/E14.5/E18.5)
meta["stage"] = meta["stage"].astype(str).str.strip()
stage_label_map = {"125": "E12.5", "145": "E14.5", "185": "E18.5"}
meta["stage"] = meta["stage"].replace(stage_label_map)
stage_num_map = {"E12.5": 12.5, "E14.5": 14.5, "E18.5": 18.5}
meta["stage_num"] = meta["stage"].map(stage_num_map)

# Standardize 'condition' labels
meta["condition"] = meta["condition"].astype(str).str.strip()
cond_map = {"WT": "WT", "NT": "SCNT"}
meta["condition"] = meta["condition"].replace(cond_map)

print("--- Metadata Summary ---")
print(meta["condition"].value_counts(), "\n")
print(meta["stage"].value_counts(), "\n")
print(meta["sample_id"].value_counts(), "\n")


# Generate Count and Design Matrices -------------------------------------------
count_tbl = pd.crosstab(meta["sample_id"], meta["cell_type"]).astype(int)

# Generate sample-level design matrix
design = (
    meta[["sample_id", "condition", "stage", "stage_num"]]
    .drop_duplicates()
    .set_index("sample_id")
    .loc[count_tbl.index]
)


# Construct scCODA AnnData Object ----------------------------------------------
df = pd.concat([design[["condition"]], count_tbl], axis=1)
data = dat.from_pandas(df, covariate_columns=["condition"])

# Specify category order for condition (Baseline = WT, Effect = SCNT)
data.obs["condition"] = pd.Categorical(data.obs["condition"], categories=["WT", "SCNT"])


# Exploratory Data Visualization -----------------------------------------------
fig_dir = "/home4/ssyi/Mouse_Placenta/Figures"

# Stacked barplot for the levels of condition
viz.stacked_barplot(data, feature_name="condition")
plt.savefig(os.path.join(fig_dir, "E145_scCODA_stacked_barplot.pdf"), bbox_inches="tight")
# plt.show()
plt.close()

# Boxplots for relative abundance
viz.boxplots(
    data, 
    feature_name="condition",
    figsize=(8, 5),
    add_dots=False,
    args_swarmplot={"palette": ["black"]}
)
plt.savefig(os.path.join(fig_dir, "E145_scCODA_boxplots.pdf"), bbox_inches="tight")
# plt.show()
plt.close()


# Reference Cell Type Selection ------------------------------------------------
# Calculate proportions and stability metrics
prop = count_tbl.div(count_tbl.sum(axis=1), axis=0)
presence = (count_tbl > 0).mean(axis=0)         
dispersion = prop.var(axis=0)                   

# Rank candidates by presence and variance
ref_candidates = pd.DataFrame({
    "presence": presence,
    "var": dispersion,
    "mean_prop": prop.mean(axis=0),
    "zeros": (count_tbl == 0).sum(axis=0)
}).sort_values(["presence", "var"], ascending=[False, True])

print("--- Rank Reference Candidates by Stability ---")
print(ref_candidates.head(10), "\n")

# Calculate Coefficient of Variation (CV). Lower is more stable.
cv = (prop.std(axis=0) / (prop.mean(axis=0) + 1e-9)).sort_values()
print("--- Coefficient of Variation (CV) ---")
print(cv.head(10), "\n")


# Robustness Testing (Iterative MCMC) ------------------------------------------
all_results_dict = {}
cell_types = data.var.index
results_cycle = pd.DataFrame(index=cell_types, columns=["times_credible"]).fillna(0)

# Iterate each cell type as the reference to ensure robustness
for ct in cell_types:
    print(f"Running MCMC with Reference: {ct}")
    try:
        # Initialize and run inference
        model_temp = mod.CompositionalAnalysis(
          data, 
          formula="condition", 
          reference_cell_type=ct, 
          automatic_reference_absence_threshold=0.05
        )
        temp_results = model_temp.sample_hmc(
          num_results=20000,
          num_burnin=5000,
          num_leapfrog_steps=10,
          step_size=0.01,
          verbose=False
        )
        all_results_dict[ct] = temp_results
        eff_df = temp_results.effect_df
    
        safe_ct_name = str(ct).replace("/", "_").replace(" ", "_")
        csv_filename = f"/home4/ssyi/Mouse_Placenta/codes/Review/3.1.scCODA/{safe_ct_name}_ref_scCODA.csv"
        eff_df.to_csv(csv_filename)

        cred_eff = temp_results.credible_effects()
        cred_eff.index = cred_eff.index.droplevel(level=0)
        results_cycle["times_credible"] = results_cycle["times_credible"].add(cred_eff.astype("int"), fill_value=0)
      
    except Exception as e:
        print(f"ERROR: Failed for reference cell type '{ct}': {e}")
        continue

# Calculate credibility percentages
results_cycle["pct_credible"] = results_cycle["times_credible"] / len(cell_types)
results_cycle["is_credible"] = results_cycle["pct_credible"] > 0.5

print("\n--- Credible Effects Across All References ---")
print(results_cycle.sort_values("times_credible", ascending=False), "\n")


# Final MCMC Inference & Result Export -----------------------------------------
# Run final model
model_final = ana.CompositionalAnalysis(
    data,
    formula="condition",
    reference_cell_type="automatic",
    automatic_reference_absence_threshold=0.05
)

result_final = model_final.sample_hmc(
    num_results=20000,
    num_burnin=5000,
    num_leapfrog_steps=10,
    step_size=0.01,
    verbose=True
)

# Adjust FDR threshold
result_final.set_fdr(est_fdr=0.4)

# Extract and save the final effect dataframe
final_eff_df = result_final.effect_df
out_path = "/home4/ssyi/Mouse_Placenta/codes/Review/3.1.scCODA/scCODA_result.csv"
final_eff_df.to_csv(out_path)

