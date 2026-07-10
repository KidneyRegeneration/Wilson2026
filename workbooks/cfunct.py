import numpy as np
import scanpy as sc
import anndata as ad
import pandas as pd
from scipy.spatial.distance import pdist, squareform, euclidean, cdist
from scipy.stats import zscore, pearsonr
import matplotlib.pyplot as plt
import networkx as nx
import plotly.graph_objects as go
import seaborn as sns
from collections import defaultdict


def get_nz_median(arr):
    """Returns non-zero, non-NaN median of an array."""
    arr = arr[~np.isnan(arr)]
    if not np.any(arr):
        return 1
    return np.median(arr[np.nonzero(arr)])


def normalize_data(mat):
    """Normalizes matrix by dividing each column by the non-zero median of that column."""
    nzm = np.apply_along_axis(get_nz_median, 0, mat)
    print(nzm)
    nzm[nzm == 0] = 1
    return mat / nzm


def normalize_with_reference_set(adata1, adata2):
    """Normalizes two AnnData objects together, returns split normalized sets."""
    adata = sc.concat([adata1, adata2], axis=0)
    print('normalizing')
    adata.X = normalize_data(adata.X)
    adata1_normalized = adata[:len(adata1)].copy()
    adata2_normalized = adata[len(adata1):].copy()
    return adata1_normalized, adata2_normalized


def generate_colors(n, colormap='viridis'):
    cm = plt.get_cmap(colormap)
    return [cm(1. * i / n) for i in range(n)]


def make_hist(counts):
    n_bins = int(len(counts) / 30)
    x_min, x_max = np.min(counts), np.max(counts)
    ys, bin_edges = np.histogram(counts, bins=n_bins, range=(x_min, x_max))
    if x_min == x_max:
        x_min -= 0.05
        x_max += 0.05
    bin_width = (x_max - x_min) / n_bins
    xs = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    return xs, ys / bin_width, (x_min + x_max) / 2.0


def cluster_names(data: ad.AnnData):
    """Returns ordered unique cluster names from obs['cat_cluster'] of an AnnData object."""
    unique_values, indices = np.unique(data.obs['cat_cluster'], return_index=True)
    return unique_values[np.argsort(indices)]


def remove_small_clusters(adata, min_size):
    """Separates clusters into those above and below a size threshold."""
    min_size = min_size
    clusters = cluster_names(adata)
    large_enough, removed = [], []
    for cluster in clusters:
        num_obs = (adata.obs['cat_cluster'] == cluster).sum()
        print(f"{cluster} has {num_obs} observations")
        if num_obs >= min_size:
            large_enough.append(cluster)
        else:
            print(f"{cluster} Removed")
            removed.append(cluster)
    return large_enough, removed


def add_sizes_to_cluster_names(adata):
    """Returns cluster names with sizes appended."""
    clusters = cluster_names(adata)
    return [f"{c}_{(adata.obs['cat_cluster'] == c).sum()}" for c in clusters]


def add_sizes_to_cluster_names2(adata):
    """Returns cluster names with sizes and an array of sizes."""
    clusters = cluster_names(adata)
    clusters_with_sizes = []
    cluster_sizes = []
    for cluster in clusters:
        size = (adata.obs['cat_cluster'] == cluster).sum()
        clusters_with_sizes.append(f"{cluster}_{size}")
        cluster_sizes.append(size)
    return clusters_with_sizes, cluster_sizes


def make_sankey_ordered(labels, sources, targets, values, title):
    """Creates an ordered Sankey diagram."""
    link_count = defaultdict(lambda: {'source': 0, 'target': 0})
    for src, tgt in zip(sources, targets):
        link_count[src]['source'] += 1
        link_count[tgt]['target'] += 1

    def sorting_criteria(link):
        src, tgt = link[0], link[1]
        src_c, tgt_c = link_count[src], link_count[tgt]
        return (
            src_c['source'] == 1 and tgt_c['target'] == 1,
            src_c['source'] == 2 and tgt_c['target'] == 0,
            src_c['source'] == 2,
            src_c['source'] > 2 and tgt_c['target'] == 0
        )

    sorted_links = sorted(zip(sources, targets, values), key=sorting_criteria, reverse=True)
    sorted_sources, sorted_targets, sorted_values = zip(*sorted_links)

    fig = go.Figure(go.Sankey(
        node={"label": labels, 'pad': 10, 'thickness': 15},
        link={"source": sorted_sources, "target": sorted_targets, "value": sorted_values}
    ))

    fig.update_layout(title_text=title, font_size=13, height=900)
    fig.write_html("sankey_diagram_interset.html")


def display_network(labels, sources, targets):
    G = nx.DiGraph()
    G.add_nodes_from(labels)
    G.add_edges_from((labels[s], labels[t]) for s, t in zip(sources, targets))
    pos = nx.spring_layout(G, k=0.5)

    plt.figure(figsize=(12, 8))
    nx.draw(G, pos, with_labels=True, node_size=900, font_size=14,
            width=3, edge_color="gray", alpha=0.9, arrowsize=20)
    plt.title("Network Diagram")
    plt.show()


def display_network2(labels, sources, targets):
    G = nx.DiGraph()
    G.add_nodes_from(labels)
    G.add_edges_from((labels[s], labels[t]) for s, t in zip(sources, targets))
    pos = nx.spring_layout(G, k=0.35)
    colors = plt.cm.rainbow(np.linspace(0, 1, len(G.nodes())))

    plt.figure(figsize=(12, 8))
    nx.draw(G, pos, with_labels=True, node_size=900, font_size=14,
            width=3, edge_color="gray", alpha=0.9, arrowsize=20, node_color=colors)
    plt.title("Network Diagram")
    plt.show()


def display_network_colornodes(labels, sources, targets, type1, type2):
    G = nx.DiGraph()
    G.add_nodes_from(labels)
    G.add_edges_from((labels[s], labels[t]) for s, t in zip(sources, targets))
    pos = nx.spring_layout(G, k=0.5)

    node_colors = [
        'salmon' if l in type1 else
        'lightblue' if l in type2 else
        (0.8, 0.8, 0.8) for l in labels
    ]

    plt.figure(figsize=(12, 8))
    nx.draw(G, pos, with_labels=True, node_size=900, node_color=node_colors,
            font_size=14, width=3, edge_color="gray", alpha=0.9, arrowsize=20)
    plt.title("Network Diagram")
    plt.show()


def histograms_distances(distrs_d, clusters_x, clusters_y):
    """Makes histograms of distances between clusters in two datasets."""
    rows_per_figure = 4
    colors = generate_colors(len(clusters_x), 'tab20')

    for fig_num in range(-(-len(clusters_y) // rows_per_figure)):
        fig, axs = plt.subplots(rows_per_figure, 1, figsize=(20, 11), constrained_layout=True)
        axs = np.array(axs).reshape(-1)

        start = fig_num * rows_per_figure
        end = min(start + rows_per_figure, len(clusters_y))

        for i, y_idx in enumerate(range(start, end)):
            axs[i].set_ylabel(clusters_y[y_idx], fontsize=15)
            for x_idx, clust_x in enumerate(clusters_x):
                if len(clusters_x) == distrs_d.shape[1]:
                    xs, ys, x_annotate = make_hist(distrs_d[:, x_idx, y_idx])
                    axs[i].plot(xs, ys, linewidth=2, color=colors[x_idx])
                    axs[i].annotate(clust_x, (x_annotate, max(ys)), color=colors[x_idx], fontsize=13)
                else:
                    raise ValueError("Length of clusters_x does not match the shape of distances.")

        for j in range(i + 1, rows_per_figure):
            axs[j].axis('off')

    plt.savefig("Histo.pdf", format="pdf", bbox_inches="tight")
    plt.show()


def qs_directed_interset(clusts_from, clusts_into, means, stds):
    """Computes q_ij metrics for directed inter-set comparisons."""
    rows = []
    for i, clust in enumerate(clusts_from):
        clust_means = means[:, i]
        clust_means = clust_means[np.nonzero(clust_means)]
        clust_means = np.sort(clust_means)
        coord_nn = np.where(means[:, i] == clust_means[0])[0][0]

        rows.append({'q_ij': 0, '<d_ij>': 0, 'std_<d_ij>': 0, 'oc': clust, 'c_j': clusts_into[coord_nn]})

        for m in clust_means[1:]:
            coord = np.where(means[:, i] == m)[0][0]
            dij = means[coord, i] - means[coord_nn, i]
            stdij = np.sqrt(stds[coord_nn, i] ** 2 + stds[coord, i] ** 2)
            rows.append({
                'q_ij': dij / stdij,
                '<d_ij>': dij,
                'std_<d_ij>': stdij,
                'oc': clust,
                'c_j': clusts_into[coord],
            })

    return pd.DataFrame(rows)


def alignment_table_interset_oneway(cut_off, dists, clusters, clusters_ref):
    """Displays alignment significance table for inter-set cluster comparison."""
    qs_df = qs_directed_interset(
        clusters, clusters_ref, np.mean(dists, axis=0), np.std(dists, axis=0)
    )
    df_sorted = qs_df.sort_values(by='q_ij')
    aligned = df_sorted[df_sorted['q_ij'] <= cut_off]
    indices = aligned.index
    combined_indices = np.unique(np.concatenate([indices, indices + 1, indices + 2]))
    combined_rows = df_sorted.loc[combined_indices]

    for clust in clusters:
        print(clust)
        print(combined_rows[combined_rows['oc'] == clust].drop(columns=['oc']))
        print()


def alignment_diagram_interset_oneway(cut_off, distances, clusters, reference_clusters, remove_this, normalize):
    """Creates Sankey diagram to visualize significant inter-set cluster alignments."""
    means = np.mean(distances, axis=0)
    sources, targets, values = [], [], []

    qs_df = qs_directed_interset(
        clusters.copy(), reference_clusters.copy(),
        means, np.std(distances, axis=0)
    )

    for c, cluster in enumerate(clusters):
        filtered = qs_df[qs_df['oc'] == cluster][qs_df['q_ij'] <= cut_off]
        for _, row in filtered.iterrows():
            j = np.where(reference_clusters == row['c_j'])[0][0]
            sources.append(c)
            targets.append(len(clusters) + j)
            values.append(1 / means[j, c])

    labels = list(clusters) + list(reference_clusters)
    np.save('labels_1to2.npy', labels)
    np.save('targets_1to2.npy', targets)
    np.save('sources_1to2.npy', sources)
    np.save('values_1to2.npy', values)

    if remove_this.size == 0:
        title = "CAT. "
    else:
        title = f"CAT. Removed clusters with fewer than 50 cells: {remove_this}"

    if normalize:
        title += " Non-zero median normalized"

    make_sankey_ordered(labels, sources, targets, values, title)


from collections import defaultdict
import plotly.graph_objects as go
import numpy as np # Assuming numpy is available for array operations

def make_sankey_ordered2(labels, sources, targets, values, title, node_weights=None):
    """
    Makes an ordered Sankey diagram where each connection is colored by its source node.
    """
    
    # Create a dictionary to count links for each node
    link_count = defaultdict(lambda: {'source': 0, 'target': 0})
    for src, tgt in zip(sources, targets):
        link_count[src]['source'] += 1
        link_count[tgt]['target'] += 1
    
    # Define sorting criteria
    def sorting_criteria(link):
        src, tgt = link[0], link[1]
        src_count, tgt_count = link_count[src], link_count[tgt]
        is_SL = src_count['source'] == 1 and tgt_count['target'] == 1
        is_S_two_T_not_target = src_count['source'] == 2 and tgt_count['target'] == 0
        is_S_two_T_can_be_target = src_count['source'] == 2
        is_S_more_than_two_T_not_target = src_count['source'] > 2 and tgt_count['target'] == 0
        return (is_SL, is_S_two_T_not_target, is_S_two_T_can_be_target, is_S_more_than_two_T_not_target)

    # Sort links based on criteria
    sorted_links = sorted(zip(sources, targets, values), key=sorting_criteria, reverse=True)
    
    # Unzip the sorted links
    sorted_sources, sorted_targets, sorted_values = zip(*sorted_links)

    # --- New: Color connections by source node ---
    # Define a color palette. You can expand this list for more distinct colors.
    # These are semi-transparent for better visualization of overlaps.
    color_palette = [
        'rgba(31, 119, 180, 0.8)',  # Blue
        'rgba(255, 127, 14, 0.8)',  # Orange
        'rgba(44, 160, 44, 0.8)',   # Green
        'rgba(214, 39, 40, 0.8)',   # Red
        'rgba(148, 103, 189, 0.8)', # Purple
        'rgba(140, 86, 75, 0.8)',   # Brown
        'rgba(227, 119, 194, 0.8)', # Pink
        'rgba(127, 127, 127, 0.8)', # Gray
        'rgba(188, 189, 34, 0.8)',  # Olive
        'rgba(23, 190, 207, 0.8)'   # Cyan
    ]

    # Create a mapping from unique source node indices to colors
    # The 'sources' list contains integer indices corresponding to the 'labels' list.
    unique_source_indices = np.unique(sources)
    source_index_to_color = {
        idx: color_palette[i % len(color_palette)]
        for i, idx in enumerate(unique_source_indices)
    }

    # Generate the list of colors for each link based on its source node
    link_colors = [source_index_to_color[src] for src in sorted_sources]
    
    node_thickness = 15 # Default thickness if weights aren't provided
    
    if node_weights is not None and len(node_weights) == len(labels):
        # Calculate a scale factor to prevent very small nodes from disappearing
        # and very large nodes from dominating. A simple linear scaling works best.
        max_weight = max(node_weights)
        min_weight = min(node_weights) if min(node_weights) > 0 else 1
        
        # Scale the weights to a reasonable thickness range (e.g., 10 to 40)
        # Using a fixed scaling factor based on the max is common, or a min/max normalization.
        # Let's use a simple scaling relative to the max, but ensure a minimum size.
        
        scaled_weights = [(w / max_weight) * 30 + 10 for w in node_weights] # Scales to 10-40 range
        node_thickness = scaled_weights

    fig = go.Figure(go.Sankey(
        node={
            "label": labels,
            'pad': 10,
            'thickness': node_thickness,
            # You could also color nodes by their type/group if needed
            # "color": [source_index_to_color[i] if i < len(labels)/2 else 'lightgray' for i in range(len(labels))]
        },
        link={
            "source": sorted_sources,
            "target": sorted_targets,
            "value": sorted_values,
            "color": link_colors # Apply the generated colors here
        }
    ))

    fig.update_layout(title_text=title, font_size=13, height=900)
    fig.write_html("sankey_diagram_interset.html")

    # fig.show() # Keep commented out as per original