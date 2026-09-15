"""
Fix for EBM Interaction Visualization Issues

Problem: The interpret library's visualization shows all interactions as constant
(single color), but the data actually has proper variation. This is a visualization
bug, not a data problem.

This script provides a workaround to create proper visualizations.
"""

import pickle
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def verify_interaction_data(model, interaction_name):
    """Verify that an interaction has varying values"""
    idx = model.term_names_.index(interaction_name)
    scores = model.term_scores_[idx]
    
    print(f"Interaction: {interaction_name}")
    print(f"  Shape: {scores.shape}")
    print(f"  Range: [{np.min(scores):.6f}, {np.max(scores):.6f}]")
    print(f"  Mean: {np.mean(scores):.6f}")
    print(f"  Std: {np.std(scores):.6f}")
    print(f"  Unique values: {len(np.unique(scores))}")
    print(f"  Non-zero entries: {np.count_nonzero(scores)}/{scores.size}")
    print()

def plot_interaction(model, interaction_name, save_path=None):
    """
    Create a proper heatmap visualization for an interaction term
    
    Args:
        model: Trained EBM model
        interaction_name: Name of interaction (e.g., "fano & corrCoef")
        save_path: Optional path to save figure
    """
    # Find the interaction
    idx = model.term_names_.index(interaction_name)
    scores = model.term_scores_[idx]
    feature_indices = model.term_features_[idx]
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, 9))
    
    # Create heatmap with proper colorscale
    im = ax.imshow(scores.T, aspect='auto', cmap='RdBu_r', origin='lower')
    
    # Labels
    ax.set_title(f'Interaction Effect: {interaction_name}', 
                 fontsize=16, fontweight='bold', pad=20)
    ax.set_xlabel(model.feature_names_in_[feature_indices[0]], fontsize=13)
    ax.set_ylabel(model.feature_names_in_[feature_indices[1]], fontsize=13)
    
    # Add colorbar with label
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Effect on log-odds', fontsize=12, rotation=270, labelpad=20)
    
    # Add statistics text
    stats_text = (f'Range: [{np.min(scores):.3f}, {np.max(scores):.3f}]\\n'
                  f'Mean: {np.mean(scores):.3f}\\n'
                  f'Std: {np.std(scores):.3f}')
    ax.text(1.15, 0.5, stats_text, transform=ax.transAxes,
            fontsize=10, verticalalignment='center',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Saved to {save_path}")
    
    return fig

def plot_all_interactions(model, save_dir=None):
    """Plot all interaction terms from a model"""
    import os
    
    interaction_indices = [i for i, name in enumerate(model.term_names_) 
                          if " & " in name]
    
    print(f"Found {len(interaction_indices)} interactions to plot\\n")
    
    for idx in interaction_indices:
        interaction_name = model.term_names_[idx]
        print(f"Plotting: {interaction_name}")
        
        fig = plot_interaction(model, interaction_name)
        
        if save_dir:
            os.makedirs(save_dir, exist_ok=True)
            safe_name = interaction_name.replace(" & ", "_and_").replace(" ", "_")
            save_path = os.path.join(save_dir, f"interaction_{safe_name}.png")
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            print(f"  Saved to {save_path}")
            plt.close(fig)
        else:
            plt.show()
        
        print()

# Example usage
if __name__ == "__main__":
    # Load models
    with open("ebm_results_adults/ebm_model.pkl", "rb") as f:
        adults_model = pickle.load(f)
    
    print("=" * 70)
    print("VERIFICATION: Interactions Have Varying Data")
    print("=" * 70)
    print()
    
    # Verify a few interactions have variationverify_interaction_data(adults_model, "fano & corrCoef")
    verify_interaction_data(adults_model, "peakLoc & raisingSlope")
    verify_interaction_data(adults_model, "duration & fano")
    
    print("=" * 70)
    print("CREATING CUSTOM VISUALIZATIONS")
    print("=" * 70)
    print()
    
    # Plot the problematic interaction from the user's screenshot
    plot_interaction(adults_model, "fano & corrCoef", 
                    save_path="interaction_fano_corrCoef_fixed.png")
    plt.show()
    
    # Optionally plot all interactions
    # Uncomment to generate all interaction plots:
    # plot_all_interactions(adults_model, save_dir="interaction_plots_adults")
