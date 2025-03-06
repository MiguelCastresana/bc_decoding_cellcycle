import yaml
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from adjustText import adjust_text  # optional; install via pip if needed

# Load configuration from external file
with open("~/src/repository/repository_final/GRN/config_pentapartitenetwork.yaml", "r") as f:
    config = yaml.safe_load(f)

# Load the data table from Excel
table_path = config["table_path"]
df = pd.read_excel(table_path)

# Automatically select all unique drug names for the given database
selected_drugs = sorted(df[df["database"] == config["database"]]["Drug"].unique())
print(f"Automatically selected drugs: {selected_drugs}")

# Filter the DataFrame based on the chosen database and selected drugs
df_filtered = df[(df["database"] == config["database"]) & (df["Drug"].isin(selected_drugs))].copy()

# Define the columns to include in the network (from config)
columns_to_include = config["columns_to_include"]

# Define custom ordering for specific columns if provided
custom_order = config.get("custom_order", {})

def get_ordered_values(series, col_name):
    unique_vals = series.unique()
    if col_name in custom_order:
        # Use the custom order from config, only include values present in the data
        return [val for val in custom_order[col_name] if val in unique_vals]
    else:
        return sorted(unique_vals)

# ------------------------------------------------------------------
# Create the Graph
# ------------------------------------------------------------------
G = nx.Graph()
# Define shapes and colors for each column
shape_color_map = {
    "Drug": ("o", "lightcoral"),
    "Gene": ("s", "lightskyblue"),
    "Target_Pathway": ("^", "lightgreen"),
    "BC_subtype": ("p", "mediumpurple"),
    "cell_cycle_phase": ("D", "orange")
}

pos = {}
column_items = {}
max_count = 0

# Determine ordered unique values for each column and compute maximum count
for col in columns_to_include:
    ordered_vals = get_ordered_values(df_filtered[col], col)
    column_items[col] = ordered_vals
    if len(ordered_vals) > max_count:
        max_count = len(ordered_vals)

# Add nodes for each column in a row-based layout (horizontally centered)
for col in columns_to_include:
    items = column_items[col]
    n_items = len(items)
    offset = (max_count - n_items) / 2  # centering offset
    y_coord = config["row_y_positions"].get(col, 0)
    for i, val in enumerate(items):
        node_id = f"{col}_{val}"
        x_coord = i + offset
        shape, color = shape_color_map[col]
        G.add_node(node_id, label=val, shape=shape, color=color)
        pos[node_id] = (x_coord, y_coord)

# Add edges between consecutive columns for each row in the filtered data
for _, row in df_filtered.iterrows():
    for i in range(len(columns_to_include) - 1):
        col1 = columns_to_include[i]
        col2 = columns_to_include[i + 1]
        node1 = f"{col1}_{row[col1]}"
        node2 = f"{col2}_{row[col2]}"
        if node1 in G and node2 in G:
            G.add_edge(node1, node2)

fig, ax = plt.subplots(figsize=(12, 6))

# Draw nodes by shape on ax
for shape in set(nx.get_node_attributes(G, "shape").values()):
    nodes = [n for n in G.nodes if G.nodes[n]["shape"] == shape]
    colors = [G.nodes[n]["color"] for n in nodes]
    nx.draw_networkx_nodes(G, pos, nodelist=nodes, node_shape=shape, 
                           node_color=colors, node_size=800, ax=ax)

# Draw edges on ax
nx.draw_networkx_edges(G, pos, edgelist=G.edges(), width=1.0, alpha=0.5, ax=ax)

# Draw labels at the node centers and store text objects
texts = []
for n, (x, y) in pos.items():
    t = ax.text(x, y, G.nodes[n]["label"], fontsize=8,
                ha="center", va="center", zorder=5)
    texts.append(t)

# Adjust text positions to reduce overlap
adjust_text(texts, arrowprops=dict(arrowstyle="-", color="gray", lw=0.5),
            force_text=0.1, expand_text=(1,1))

# Reset all label positions to be exactly at their node centers
for n, t in zip(G.nodes, texts):
    t.set_position(pos[n])

ax.axis("off")

out_path = config["output_path"]
fig.savefig(out_path, format="pdf", bbox_inches="tight")
plt.show()
print(f"Saved figure to {out_path}")





