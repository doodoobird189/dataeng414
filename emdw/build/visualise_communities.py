import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D

def load_data():
    """Load predictions, ground truth, and adjacency matrix"""
    # Load predictions
    pred_nodes = []
    pred_communities = []
    with open('community_predictions.txt') as f:
        for line in f:
            node, comm = map(int, line.strip().split())
            pred_nodes.append(node)
            pred_communities.append(comm)
    
    # Load ground truth
    true_nodes = []
    true_communities = []
    with open('community_truth.txt') as f:
        for line in f:
            parts = list(map(float, line.strip().split()))
            true_nodes.append(len(true_nodes))  # node index
            true_communities.append(np.argmax(parts) + 1)  # convert one-hot to label
    
    # Load adjacency matrix
    adj_matrix = np.loadtxt('adjacency_matrix.txt')
    return adj_matrix, np.array(pred_communities), np.array(true_communities)

def visualize_comparison(adj_matrix, pred_comms, true_comms):
    """Visualize predicted vs ground truth communities"""
    G = nx.from_numpy_array(adj_matrix)
    pos = nx.spring_layout(G, seed=42)  # Consistent layout
    
    # Create combined plot
    plt.figure(figsize=(14, 6))
    
    # Plot 1: Predicted communities
    plt.subplot(121)
    plot_communities(G, pos, pred_comms, "Predicted Communities")
    
    # Plot 2: Ground truth
    plt.subplot(122)
    plot_communities(G, pos, true_comms, "Ground Truth Communities")
    
    plt.tight_layout()
    plt.savefig('community_comparison.png', dpi=300, bbox_inches='tight')
    plt.show()

def plot_communities(G, pos, communities, title):
    """Helper function to plot communities"""
    unique_comms = np.unique(communities)
    colors = plt.cm.tab10(np.linspace(0, 1, len(unique_comms)))
    
    # Map communities to colors
    node_colors = [colors[np.where(unique_comms == c)[0][0]] for c in communities]
    
    nx.draw_networkx_nodes(
        G, pos,
        node_color=node_colors,
        node_size=200,
        alpha=0.9
    )
    nx.draw_networkx_edges(
        G, pos,
        width=0.5,
        alpha=0.2,
        edge_color='gray'
    )
    
    # Create custom legend
    legend_elements = [Line2D([0], [0], marker='o', color='w', 
                      markerfacecolor=colors[i], markersize=10,
                      label=f'Community {comm}')
                     for i, comm in enumerate(unique_comms)]
    plt.legend(handles=legend_elements)
    plt.title(title)
    plt.axis('off')

if __name__ == "__main__":
    adj_matrix, pred_comms, true_comms = load_data()
    visualize_comparison(adj_matrix, pred_comms, true_comms)