import sys
import re
from ete3 import Tree

def extract_newick_trees(file_content):
    trees = {}
    content = ''.join(file_content)
    matches = re.findall(r'(dS|dN) tree:(.*?)(?=\n(?:d|w)|$)', content, re.DOTALL)

    for key, tree in matches:
        tree = tree.strip()
        trees[key] = tree

    return trees


def get_tab_split_table(tree_str):
    tree = Tree(tree_str)
    table = {}

    def traverse(node):
        if node.is_leaf():
            table[node.name]=node.dist
        else:
            child_names = []
            for child in node.get_children():
                traverse(child)
                child_names.extend(child.name.split('_') if not child.is_leaf() else [child.name])
            node.name = '_'.join(child_names)
            table[node.name]=node.dist

    traverse(tree)
    return table
input_file = sys.argv[1]
output_file = sys.argv[2]
test_id = sys.argv[3]
# Read the input file content
file_content = open(input_file, 'r').readlines()

# Extract the Newick trees and save them in a dictionary
newick_trees = extract_newick_trees(file_content)

# Process each Newick tree and create the tables
tables = {key: get_tab_split_table(tree) for key, tree in newick_trees.items()}

# Combine the tables and write the result into a tab-separated file
with open(output_file, 'w') as outfile:
    # Write the header
    outfile.write("nodename\tdN\tdS\tid\n")
    for nodename in tables['dN']:
        dN = tables['dN'][nodename]
        dS = tables['dS'][nodename]
        outfile.write(f"{nodename}\t{dN}\t{dS}\t{test_id}\n")
