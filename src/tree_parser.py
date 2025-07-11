from Bio import Phylo

def load_tree(path):
    """
    Load a phylogenetic tree from a file in Newick format.
    """
    tree = Phylo.read(path, 'newick')
    # Ensure all clades have names
    counter = 0
    for clade in tree.find_clades():
        if clade.name is None:
            clade.name = f'Node{counter}'
            counter += 1
    return tree
