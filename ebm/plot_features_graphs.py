import os
import pickle
import matplotlib.pyplot as plt
import numpy as np


def plot_ebm_shape_safe(ebm_global_obj, feature_name, save_path):
    try:
        feat_idx = ebm_global_obj.term_names_.index(feature_name)
    except ValueError:
        print(f"Skipping {feature_name} (likely an interaction term)")
        return

    # Use explain_global().data() to get the correct plotting data
    explanation = ebm_global_obj.explain_global()
    feat_data = explanation.data(feat_idx)
    
    x_vals = feat_data['names']  # Bin edges
    y_vals = feat_data['scores']  # Scores for each bin
    
    plt.figure(figsize=(8, 5))
    
    # Check if categorical (strings) or continuous (numbers)
    is_categorical = False
    if len(x_vals) > 0:
        # If the first element is a string that CANNOT be converted to float
        try:
            float(x_vals[0])
        except (ValueError, TypeError):
            is_categorical = True

    if is_categorical:
        plt.bar(x_vals, y_vals, color='darkblue', alpha=0.7)
    else:
        # Convert to numpy arrays
        x_vals = np.array(x_vals, dtype=float)
        y_vals = np.array(y_vals, dtype=float)
        
        # Handle dimension mismatch: x_vals (bin edges) is typically 1 longer than y_vals (scores)
        # For step plot with where='post', we use the bin starts
        if len(x_vals) == len(y_vals) + 1:
            x_vals = x_vals[:-1]
        
        plt.step(x_vals, y_vals, color='darkblue', linewidth=2, where='post')

    plt.axhline(0, color='gray', linestyle='--')
    plt.title(f"EBM Shape: {feature_name}")
    plt.ylabel("Score (Log-Odds contribution)")
    plt.xlabel(feature_name)
    plt.grid(True, alpha=0.3)
    
    plt.savefig(save_path)
    plt.close()
    print(f"Saved plot for {feature_name}")

# --- 1. Load the Model ---
BASE_DIR = "ebm_results_children"
MODEL_PATH = os.path.join(BASE_DIR, "ebm_model.pkl")
print(f"Loading model from: {MODEL_PATH}")
with open(MODEL_PATH, "rb") as f:
    global_ebm = pickle.load(f)
print("Model loaded successfully.\n")

# --- 2. Plot all feature shapes ---
OUTPUT_DIR = os.path.join(BASE_DIR, "feature_plots")
os.makedirs(OUTPUT_DIR, exist_ok=True)
print(f"Plotting feature shapes to directory: {OUTPUT_DIR}\n")
for i, feature_name in enumerate(global_ebm.term_names_):
    if " x " in feature_name or " & " in feature_name:
        print(f"Skipping interaction term: {feature_name}")
        continue
    print(f"Plotting feature {i+1}/{len(global_ebm.term_names_)}: {feature_name}")
    save_path = os.path.join(OUTPUT_DIR, f"feature_{i+1}_{feature_name.replace(' ', '_')}.png")
    plot_ebm_shape_safe(global_ebm, feature_name, save_path=save_path)