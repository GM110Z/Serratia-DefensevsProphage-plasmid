#!/usr/bin/env python3
import sys
import pandas as pd
from collections import defaultdict, deque

# Get arguments
args = sys.argv
mash_dist_file = args[1]  # e.g., distances.txt from mash dist -t
mash_threshold = float(args[2])  # e.g., 0.01 for near-identical Serratia strains
output_file = args[3]  # Output: genome-to-cluster table

# Load MASH distance matrix
# Use header=0 to treat first row as column names, index_col=0 for row names
df = pd.read_csv(mash_dist_file, sep="\t", header=0, index_col=0, low_memory=False)

# Convert matrix to pairwise comparisons
pairs = df.stack().reset_index()
pairs.columns = ["genome1", "genome2", "distance"]

# Ensure distance is float
pairs["distance"] = pairs["distance"].astype(float)

# Create undirected graph based on threshold
edges = defaultdict(set)
all_genomes = set()
for _, row in pairs.iterrows():
    g1 = row["genome1"].rsplit("/", 1)[-1].replace(".fasta", "")
    g2 = row["genome2"].rsplit("/", 1)[-1].replace(".fasta", "")
    all_genomes.update([g1, g2])
    if row["distance"] <= mash_threshold:
        edges[g1].add(g2)
        edges[g2].add(g1)

# Cluster genomes using connected components
visited = set()
cluster_id = 0
clusters = {}
for genome in sorted(all_genomes):
    if genome in visited:
        continue
    cluster_id += 1
    queue = deque([genome])
    while queue:
        g = queue.popleft()
        if g in visited:
            continue
        visited.add(g)
        clusters[g] = cluster_id
        queue.extend(edges[g] - visited)

# Fill in singleton genomes
for genome in all_genomes:
    if genome not in clusters:
        cluster_id += 1
        clusters[genome] = cluster_id

# Write output
out_df = pd.DataFrame({
    "genome": sorted(clusters.keys()),
    "cluster_id": [clusters[g] for g in sorted(clusters.keys())]
})
out_df.to_csv(output_file, sep="\t", index=False)
