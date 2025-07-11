import argparse
from tree_parser import load_tree
from matrix_parser import parse_tnt_matrix
from fitch import down_pass, up_pass, detect_changes, detect_semiapomorphy, drop_last_character
from visualize import draw_svg, draw_tree_withname

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--tree', required=True, help='Newick tree file')
    p.add_argument('--matrix', required=True, help='TNT matrix file')
    p.add_argument('--out-text', required=True, help='output text file')
    p.add_argument('--out-svg', required=True, help='output SVG file')
    p.add_argument('--out-orisvg', required=True, help='output original SVG file with node names')
    args = p.parse_args()

    tree = load_tree(args.tree)
    tip_states = parse_tnt_matrix(args.matrix)

    # down pass - up pass - detect changes - detect semi-apomorphy - drop last character
    down, full = down_pass(tree, tip_states)
    final = up_pass(tree, down, full)
    changes = detect_changes(tree, final, tip_states)
    changes = detect_semiapomorphy(changes, final)
    changes = drop_last_character(changes)

    # write output to text file
    with open(args.out_text, 'w') as f:
        # progres report
        f.write("Node\t↓n\t↓B\t↑N\n")
        for clade in tree.find_clades(order='level'):
            n = [''.join(sorted(s)) for s in down[clade.name]['n']]
            B = [''.join(sorted(s)) for s in down[clade.name]['B']]
            N = [''.join(sorted(s)) for s in final[clade.name]]
            f.write(f"{clade.name}\t{','.join(n)}\t{','.join(B)}\t{','.join(N)}\n")

        f.write("\n# Character changes (only list changed characters)\n")
        for (u, v), info in changes.items():
            # filter out unchanged characters
            changed = [i for i, m in enumerate(info['mask']) if m]
            if not changed:
                continue
            # use mask to get the indices of changed characters
            pos_str = ",".join(str(i) for i in changed)
            # compare old and new states
            change_strs = []
            for i in changed:
                #old, new = info['status'][i]
                new = info['status'][i]
                #change_strs.append(f"{old}->{new}")
                change_strs.append(f"{new}")
            status_str = ",".join(change_strs)
            # write the change information
            f.write(f"{u} -> {v}:\tcharacter[{pos_str}]\tstates[{status_str}]\n")
        # write apomorphies, semi-apomorphies and homoplasies
        f.write("\n# Apomorphies\n")
        for (u, v), info in changes.items():
            # filter out type marked as apomorphy
            apos = [i for i, t in enumerate(info['type']) if t == 'apomorphy']
            if not apos:
                continue
            # assemble the indices and new states
            pos_str = ",".join(str(i) for i in apos)
            state_str = ",".join(info['status'][i] for i in apos)
            f.write(f"{u} -> {v}:\tcharacter[{pos_str}]\tstates[{state_str}]\n")
        f.write("\n# Semi-apomorphies\n")
        for (u, v), info in changes.items():
            apos = [i for i, t in enumerate(info['type']) if t == 'semi-apomorphy']
            if not apos:
                continue
            pos_str = ",".join(str(i) for i in apos)
            state_str = ",".join(info['status'][i] for i in apos)
            f.write(f"{u} -> {v}:\tcharacter[{pos_str}]\tstates[{state_str}]\n")
        f.write("\n# Homoplasy\n")
        for (u, v), info in changes.items():
            apos = [i for i, t in enumerate(info['type']) if t == 'homoplasy']
            if not apos:
                continue
            pos_str = ",".join(str(i) for i in apos)
            state_str = ",".join(info['status'][i] for i in apos)
            f.write(f"{u} -> {v}:\tcharacter[{pos_str}]\tstates[{state_str}]\n")
    # SVG visualization
    draw_svg(tree, final, changes, args.out_svg)
    draw_tree_withname(tree, args.out_orisvg)

if __name__ == '__main__':
    main()
