#%% Data Loading and Preprocessing
import scanpy as sc
import pandas as pd

adata_processed = sc.read_h5ad("data/LUAD_TPM_normalized.h5ad")
adata = sc.read_h5ad("data/TCGA_LUAD_from_source_raw.h5ad")
adata.obs = adata_processed.obs

provided_data = pd.read_csv("data/TCGA.LUAD.expression.txt", sep="\t")
jiaming_metadata = pd.read_csv("data/TCGA_LUAD_final_for_OS_model.csv")

provided_data = provided_data.drop(columns=['patient_id']).set_index('sample_id')
provided_data.index.name = None
# Remove columns whose names contain "?"
provided_data = provided_data.loc[:, ~provided_data.columns.str.contains('\?')]
# Rename the columns to the gene names (split by "|" and keep the first part)
provided_data.columns = provided_data.columns.str.split('|').str[0]

haowen_genes = set(adata.var_names)
provided_genes = set(provided_data.columns)

# Find the intersection of the two sets
common_genes = haowen_genes.intersection(provided_genes)

# Convert set to list for proper indexing
common_genes_list = list(common_genes)

adata = adata[:, common_genes_list]
provided_data = provided_data[common_genes_list]

adata.obs = adata.obs[['sample', 'patient']]
overlapping_patients = set(jiaming_metadata["patient_id"])

# Find the intersection of the two sets (adata.obs["patient"] and overlapping_patients)
common_patients = set(adata.obs["patient"]).intersection(overlapping_patients)
# Filter the adata to only include the common patients
adata = adata[adata.obs["patient"].isin(common_patients)]
# add the jiaming_metadata to the adata.obs
adata.obs = adata.obs.merge(jiaming_metadata, left_on="patient", right_on="patient_id", how="left")

#%% DEG Analysis
import os
import pickle as pkl
import pandas as pd

from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

# Replace this with the path to directory where you would like results to be saved
OUTPUT_PATH = "deseq2_corrected_age_gender_race"
os.makedirs(OUTPUT_PATH, exist_ok=True)

counts_df = pd.DataFrame(adata.X.todense(), index=list(adata.obs["sample"]), columns=adata.var_names)
# Make sure all values mod 1 are 0 if not raise an error
if not (counts_df % 1 == 0).all().all():
    raise ValueError("Counts matrix contains non-integer values")
counts_df = counts_df.round(0).astype(int)

metadata_df = pd.DataFrame(adata.obs)
metadata_df.index = list(adata.obs["sample"])

# Convert gender to binary
all_genders = metadata_df["gender"].unique()
gender_to_int = {gender: i for i, gender in enumerate(all_genders)}
metadata_df["gender"] = metadata_df["gender"].map(gender_to_int)

# Convert race to int
all_races = metadata_df["race"].unique()
race_to_int = {race: i for i, race in enumerate(all_races)}
metadata_df["race"] = metadata_df["race"].map(race_to_int)

# metadata_df["condition"] = metadata_df["OS"]
metadata_df = metadata_df[["OS", "age", "gender", "race"]]

inference = DefaultInference(n_cpus=8)
dds = DeseqDataSet(
    counts=counts_df,
    metadata=metadata_df,
    design="~OS+age+gender+race",
    refit_cooks=True,
    inference=inference,
)

dds.deseq2()
ds = DeseqStats(dds, contrast=["OS", 0, 1], inference=inference)
ds.summary()

#%%
ds.summary(lfc_null=0.1, alt_hypothesis="greaterAbs")
ds.plot_MA(s=20)

#%%
ds.results_df.to_csv(os.path.join(OUTPUT_PATH, "deseq2_results_full.csv"))
deseq2_results_df = ds.results_df
# Sort by padj
deseq2_results_df = deseq2_results_df.sort_values(by="padj")
deseq2_results_df_filtered = deseq2_results_df[deseq2_results_df["padj"] < 0.05]
deseq2_results_df_filtered.to_csv(os.path.join(OUTPUT_PATH, "deseq2_results_padj_0.05.csv"))