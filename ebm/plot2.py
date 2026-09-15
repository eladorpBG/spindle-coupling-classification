import os
import pickle
import numpy as np

# --- 1. Load the Model ---
BASE_DIR = "ebm_results_adults"
MODEL_PATH = os.path.join(BASE_DIR, "ebm_model.pkl")

print(f"Loading model from: {MODEL_PATH}")
with open(MODEL_PATH, "rb") as f:
    model = pickle.load(f)

expl = model.explain_global()
output_dir = os.path.join(BASE_DIR, "plots_png_sorted")
os.makedirs(output_dir, exist_ok=True)

# --- 2. Extract and Sort Importances ---
# The summary (global importance) is stored at index -1
summary_data = expl.data() 

names = summary_data['names']   # List of feature names
scores = summary_data['scores'] # List of importance scores

# Zip them together so we can sort them as pairs
# We also need to map names back to their original index for visualize()
name_to_index = {name: i for i, name in enumerate(model.term_names_)}

# Create a list of tuples: (score, name, original_index)
ranked_features = []
for name, score in zip(names, scores):
    if name in name_to_index:
        original_idx = name_to_index[name]
        ranked_features.append((score, name, original_idx))

# Sort by score descending (Highest score first)
ranked_features.sort(key=lambda x: x[0], reverse=True)

# --- 3. Generate Plots ---
print(f"Generating {len(ranked_features)} plots...")

for rank, (score, name, original_idx) in enumerate(ranked_features):
    # if " x " in name or " & " in name:
    #     print(f"Skipping interaction term: {name}")
    #     continue
    # Sanitize filename
    safe_name = name.replace(" ", "_").replace("/", "-")
    
    # Filename: 01_FeatureName.png
    # rank+1 ensures it starts at 01
    filename = f"{str(rank + 1).zfill(2)}_{safe_name}.png"
    save_path = os.path.join(output_dir, filename)
    
    # Visualize using the ORIGINAL index
    fig = expl.visualize(original_idx)
    
    # Optional: Styling
    fig.update_layout(width=800, height=600, template="plotly_white")
    
    # Save
    fig.write_image(save_path)
    
    print(f"Saved Rank {rank+1}: {filename} (Score: {score:.4f})")