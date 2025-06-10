from scipy.sparse import csr_matrix
import h5py
import pandas as pd

# Step 1: Load the sparse expression matrix
with h5py.File("combined_data.h5ad", "r") as f:
    data = f["X"]["data"][:]
    indices = f["X"]["indices"][:]
    indptr = f["X"]["indptr"][:]
    X_shape = (len(indptr) - 1, indices.max() + 1)
    X_csr = csr_matrix((data, indices, indptr), shape=X_shape)
    X_dense = X_csr.toarray()

    gene_names = [g.decode("utf-8") for g in f["var"]["_index"][:X_shape[1]]]
    survival_labels = f["obs"]["survival_3yr_label"][:]
    obs_index = [i.decode("utf-8") for i in f["obs"]["_index"][:]]

# Step 2: Match dimensions
n_valid_rows = len(survival_labels)
X_dense = X_dense[:n_valid_rows, :]
gene_names = gene_names[:X_dense.shape[1]]

# Step 3: Build DataFrame
df = pd.DataFrame(X_dense, columns=gene_names)
df["label"] = survival_labels
df = df[(df["label"] != -9223372036854775808) & (~pd.isna(df["label"]))].reset_index(drop=True)

# Step 4: Reformat into prompt/response
def make_prompt(row):
    features = "; ".join([f"{col}: {row[col]:.4f}" for col in row.index if col != "label"])
    return f"Patient Record:\n{features}\nWhat is the 3-year survival outcome?\nAnswer:"

def make_response(label):
    return "Yes" if str(int(label)) == "1" else "No"

# Build Hugging Face style dataset
llm_ready_df = pd.DataFrame()
llm_ready_df["prompt"] = df.apply(make_prompt, axis=1)
llm_ready_df["response"] = df["label"].apply(make_response)

# Step 5: Save to JSONL
llm_ready_df.to_json("llm_supervised_data.jsonl", orient="records", lines=True)

print("✅ LLM-ready dataset saved to 'llm_supervised_data.jsonl'")
