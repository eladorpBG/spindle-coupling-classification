import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from interpret.glassbox import ExplainableBoostingClassifier
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.model_selection import LeaveOneGroupOut
import os
import pickle

from plot_features_graphs import plot_ebm_shape_safe

# --- Configuration ---
DATA_PATH = "~/spindles/children_data_2025.csv"
OUTPUT_DIR = "ebm_results_children"
SEED = 42
SKIP_LOSO_CV = False  # Set to True to skip LOPO-CV and only do global interpretability
SKIP_INTERPRETABILITY = False  # Set to True to skip interpretability plots

# Ensure output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

# --- 1. Load and Clean Data ---
print(f"Loading data from {DATA_PATH}...")
data = pd.read_csv(DATA_PATH)

# Remove "O" channels
data = data[~data['channel'].str.contains("O")].copy()

# Clean channel names (remove -M1, -M2)
data['channel'] = data['channel'].str.replace(r'-M[12]$', '', regex=True)

# Impute NA in 'refrPeriod' based on 'coupling_label' mean (mimicking R's na.aggregate)
data['refrPeriod'] = data['refrPeriod'].fillna(
    data.groupby('coupling_label')['refrPeriod'].transform('mean')
)

# Drop unused columns
cols_to_drop = ['numCycles', 'channel']
data = data.drop(columns=[c for c in cols_to_drop if c in data.columns])

# Define predictors and target
# Exclude metadata columns
metadata_cols = ["patient_id", "spindle_idx", "detSample", "startSample", "endSample", "coupling_label"]
feature_cols = [c for c in data.columns if c not in metadata_cols]

X = data[feature_cols]
y = data['coupling_label']
groups = data['patient_id']

print(f"Data prepared. N={len(data)}, Features={len(feature_cols)}")

# --- 2. LOPO-CV for Performance Metrics ---
print("\n--- Starting LOPO-CV ---")

lopo_auc = []
lopo_prauc = []
unique_patients = data['patient_id'].unique()
    
if not SKIP_LOSO_CV:
    # We can use SKLearn's LeaveOneGroupOut or manual loop. Manual loop allows easy printing per patient.
    for patient in unique_patients:
        train_mask = data['patient_id'] != patient
        val_mask = data['patient_id'] == patient
    
        X_train, y_train = X[train_mask], y[train_mask]
        X_val, y_val = X[val_mask], y[val_mask]
        
        # Train EBM
        ebm = ExplainableBoostingClassifier(random_state=SEED, n_jobs=1)
        ebm.fit(X_train, y_train)
        
        # Predict Probability (for class 1)
        probs = ebm.predict_proba(X_val)[:, 1]
        
        # Metrics
        if len(np.unique(y_val)) > 1:
            auc = roc_auc_score(y_val, probs)
            pr_auc = average_precision_score(y_val, probs)
            
            lopo_auc.append(auc)
            lopo_prauc.append(pr_auc)
            print(f"Patient: {patient} - AUC: {auc:.3f}")
        else:
            print(f"Patient: {patient} - Skipped (Only one class in validation set)")

print(f"Average AUC: {np.mean(lopo_auc):.3f}")
print(f"Average PR-AUC: {np.mean(lopo_prauc):.3f}")

# --- 3. Global Interpretability ---
print("\n--- Training Global Model ---")
global_ebm = ExplainableBoostingClassifier(random_state=SEED, n_jobs=-1)
global_ebm.fit(X, y)

# A. Feature Importance
print("Calculating importance...")
ebm_global = global_ebm.explain_global() # <--- This is the key fix
importances = global_ebm.term_importances()
feature_names = global_ebm.term_names_

imp_df = pd.DataFrame({'feature': feature_names, 'importance': importances})
imp_df = imp_df.sort_values('importance', ascending=False)
top_features = imp_df.head(4)['feature'].tolist()

print("Top Features:", top_features)

# B. Plotting (The Fixed Version)
print("Generating Plots...")

# Run the plotting loop
for feat in top_features:
    # Skip interaction terms for 2D plots (they contain ' x ' or ' & ')
    if " x " not in feat and " & " not in feat:
        safe_name = feat.replace(" ", "_")
        plot_ebm_shape_safe(ebm_global, feat, os.path.join(OUTPUT_DIR, f"shape_{safe_name}.png"))

print(f"\nCompleted successfully. Results in: {os.path.abspath(OUTPUT_DIR)}")

# --- 4. Save the Global Model ---
print("Saving model to disk...")
model_path = os.path.join(OUTPUT_DIR, "ebm_model.pkl")

with open(model_path, 'wb') as f:
    pickle.dump(global_ebm, f)

print(f"Model saved to: {model_path}")
print(f"To load it later, use: model = pickle.load(open('{model_path}', 'rb'))")