from io import StringIO
from Bio import Phylo
from ete3 import Tree, TreeStyle, NodeStyle, TextFace, faces, add_face_to_node
import io

def biophylo_to_ete3(tree):
    """
    transform Bio.Phylo tree to ete3 Tree object.
    """
    buf = StringIO()
    Phylo.write(tree, buf, "newick")
    nwk = buf.getvalue().strip()
    t = Tree(nwk, format=1)
    return t

def draw_svg(tree, final, changes, out_svg):
    """
    use ete3 to visualize phylogenetic tree in SVG:
    1. tree: Bio.Phylo object
    2. final: dict of {node_name: list of state-sets}
    3. changes: dict of {(parent, child): mask_list}
    4. out_svg: output SVG path
    """

    # read tree from Bio.Phylo and convert to ete3 Tree
    t = biophylo_to_ete3(tree)

    # set up tree style
    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.show_branch_length = False
    ts.show_branch_support = False

    # layout function to customize node appearance
    def layout(node):
        name = node.name or ""
        # label = node.name or ""
        st = final.get(name, [])
        label = "".join(f"[{''.join(sorted(s))}]" for s in st)
        tf = TextFace(label, fsize=8)
        faces.add_face_to_node(tf, node, column=0, position="branch-right")

        # color the node based on its state change
        if node.up:
            p = node.up.name or ""
            mask = changes.get((p, name), [])
            nstyle = NodeStyle()
            nstyle["vt_line_width"] = 2
            nstyle["hz_line_width"] = 2
            nstyle["fgcolor"] = "black"
            nstyle["bgcolor"] = "white"
            if any(mask):
                nstyle["vt_line_color"] = "red"
                nstyle["hz_line_color"] = "red"
            else:
                nstyle["vt_line_color"] = "gray"
                nstyle["hz_line_color"] = "gray"
            node.set_style(nstyle)

    ts.layout_fn = layout

    # render the tree to SVG
    t.render(out_svg, tree_style=ts, w=800, units="px")
    print(f"Tree visualization saved to {out_svg}")

def draw_tree_withname(tree, out_svg):
    """
    use ete3 to visualize Bio.Phylo tree as SVG with node names only.
    """
    # read tree from Bio.Phylo and convert to ete3 Tree
    t = biophylo_to_ete3(tree)

    # style the tree
    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.show_branch_length = False
    ts.show_branch_support = False

    # only show node names
    def layout(node):
        name = node.name or ""
        if name:
            tf = TextFace(name, fsize=10)
            faces.add_face_to_node(tf, node, column=0, position="branch-right")

    ts.layout_fn = layout

    # render the tree to SVG
    t.render(out_svg, tree_style=ts, w=800, units="px")
    print(f"Tree with node names saved to {out_svg}")

