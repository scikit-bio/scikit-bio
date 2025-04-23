import pandas as pd
import skbio
from skbio import TreeNode
import cProfile


def treenode_to_dataframe(tree_node):
   rows_for_dataframe = []
   tree_node.assign_ids()
   for node in tree_node.traverse(include_self=True):
       # if node.nazme == 'None':


       rows_for_dataframe.append(
           {'parent': node.parent.id if node.parent is not None else -1, "node": node.id, 'name': node.name,
            'length': node.length,"is_tip":node.is_tip()})
       # add in name of the node and edge length of the node
   df = pd.DataFrame(rows_for_dataframe)
   return df



def dataframe_to_treenode(dataframe):
   # Find root node
   root_row = dataframe[dataframe["parent"] == -1].iloc[0]
   root = TreeNode(name=root_row["name"]) if root_row["name"] != 'None' else TreeNode()
   root.id = root_row["node"]


   # Convert DataFrame to dictionary for fast lookup
   node_dict = {row["node"]: row for _, row in dataframe.iterrows()}


   # Create dictionary to store TreeNodes
   tree_nodes = {root.id: root}


   # Build tree structure
   for node_id, row in node_dict.items():
       if node_id == root.id:
           continue  # Skip root node


       parent_id = row["parent"]
       node_name = row["name"]
       node_length = row["length"]


       # Create TreeNode
       node = TreeNode(name=node_name if node_name != 'None' else None,
                       length=float(node_length) if pd.notna(node_length) else None)
       node.id = node_id
       tree_nodes[node_id] = node


       # Append to parent (faster than repeated DataFrame lookups)
       tree_nodes[parent_id].append(node)


   return root

def tip_to_root_conversion(dataframe, tip_node_names):
    # Ensure 'name' column is treated as string type and convert 'None' strings to actual None
    dataframe['name'] = dataframe['name'].replace('None', None)

    required_columns = {"node", "parent", "name", "length"}
    if not required_columns.issubset(dataframe.columns):
        raise ValueError(f"DataFrame is missing required columns: {required_columns - set(dataframe.columns)}")

    # Ensure unique names for nodes with names
    named_nodes = dataframe[dataframe['name'].notna()]  # nodes with names
    unnamed_nodes = dataframe[dataframe['name'].isna()]  # nodes without names

    # Ensure unique names for nodes with names
    name_counts = named_nodes['name'].value_counts()
    name_map = {}
    for name, count in name_counts.items():
        if count > 1:
            duplicates = named_nodes[named_nodes['name'] == name].index
            for i, idx in enumerate(duplicates, 1):
                name_map[idx] = f"{name}_{i}"

    # Apply the unique name mapping to named nodes only
    dataframe.loc[named_nodes.index, 'name'] = named_nodes.index.map(lambda i: name_map.get(i, dataframe.at[i, 'name']))

    # Build lookup tables
    node_dict = dataframe.set_index('node').to_dict(orient='index')  # Change indexing to 'node'
    id_to_name = dataframe.set_index('node')['name'].to_dict()
    name_to_id = {v: k for k, v in id_to_name.items()}  # Reverse mapping for names to node IDs

    tree_nodes = {}
    seen_pairs = set()

    def get_or_create_node(node_id, length=None):
        # First check if node already exists
        if node_id in tree_nodes:
            return tree_nodes[node_id]

        # Fetch node name from id_to_name, and if there's no name, leave it as None
        node_name = id_to_name.get(node_id, None)

        # Create node with no name if it's missing
        if node_name is None:
            node = TreeNode(length=float(length) if pd.notna(length) else None)  # Do not assign any name
        else:
            node = TreeNode(name=node_name, length=float(length) if pd.notna(length) else None)

        # Store the node in the dictionary for future reference
        tree_nodes[node_id] = node
        return node

    roots = set()

    for tip_name in tip_node_names:
        if tip_name not in name_to_id:
            print(f"Warning: Tip node '{tip_name}' not found in DataFrame")
            continue

        current_id = name_to_id[tip_name]  # Convert tip name to node ID

        while True:
            if current_id not in node_dict:
                print(f"Warning: Node ID '{current_id}' not found in the tree.")
                break

            row = node_dict[current_id]  # Corrected lookup
            parent_id = row["parent"]
            node_length = row["length"]

            node = get_or_create_node(current_id, node_length)

            if parent_id == -1 or pd.isna(parent_id):
                roots.add(node)
                break

            parent_node = get_or_create_node(parent_id)

            if (parent_id, current_id) not in seen_pairs:
                parent_node.append(node)
                seen_pairs.add((parent_id, current_id))

            current_id = parent_id  # Move up to the parent node

    if len(roots) > 1:
        root = TreeNode(name="root")
        root.extend(roots)
    else:
        root = roots.pop()

    return root

