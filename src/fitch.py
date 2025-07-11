from copy import deepcopy
from collections import Counter, defaultdict

def down_pass(tree, tip_states):
    """
    from tip to root:
    """
    all_states = set()

    # collect all states from tip
    for states in tip_states.values():
        for c in states:
            # only add non-unknown states
            if c != '?':
                all_states.add(c)
    k = len(next(iter(tip_states.values())))
    full_set = [set(all_states) for _ in range(k)]
    result = {}

    def recurse_down(clade):
        # if clade is terminal, initialize n and B
        if clade.is_terminal():
            seq = tip_states[clade.name]
            n_list, B_list = [], []
            for i, c in enumerate(seq):
                if c == '?':
                    # unknown state, use full set
                    n_list.append(full_set[i])
                    B_list.append(full_set[i])
                else:
                    n_list.append({c})
                    B_list.append({c})
            result[clade.name] = {'n': n_list, 'B': B_list}
            return
        for child in clade.clades:
            recurse_down(child)
        L = result[clade.clades[0].name]['n']
        R = result[clade.clades[1].name]['n']
        n = []
        B = []
        for i in range(k):
            # if both children have the same state, use intersection for n
            # also use union for B
            if L[i] & R[i]:
                n.append(L[i] & R[i])
                B.append(L[i] | R[i])
            # if they differ, use union for n and full set for B
            else:
                n.append(L[i] | R[i])
                B.append(full_set[i])
        result[clade.name] = {'n': n, 'B': B}

    recurse_down(tree.root)
    return result, full_set

def up_pass(tree, down, full_sets):
    """
    from root to tip:
    """
    final = {}

    k = len(full_sets)

    def recurse_up(clade, parent_N=None, parent_name=None):
        down_n = down[clade.name]['n']
        down_B = down[clade.name]['B']

        # if parent_N is None, this is the root node, so N = n
        if parent_N is None:
            N = deepcopy(down_n)
        else:
            # if clade is terminal, use n directly
            if clade.is_terminal():
                N = deepcopy(down_n)
            # if clade is not terminal, check if parent_N is subset of down_n
            else:
                N = []
                mask = []
                for i in range(k):
                    # if parent_N[i] is a subset of down_n[i], use parent_N[i]
                    if parent_N[i].issubset(down_n[i]):
                        Ni = parent_N[i]
                    # otherwise, use the union of n and (intersection with down_B[i])
                    else:
                        Ni = down_n[i] | (parent_N[i] & down_B[i])
                    N.append(Ni)

        final[clade.name] = N
        for ch in clade.clades:
            recurse_up(ch, N, clade.name)

    recurse_up(tree.root, None, None)
    return final


def detect_changes(tree, final, tip_states):
    """
    traversal the whole tree, comparing each parent → child N,
    if there is no intersection, mark the sites as changed,
    at the same time record the new state of child N,
    record apomorphy and homoplasy

    parameters
      tree       – Bio.Phylo.Tree
      final      – dict[node_name] = List[Set], N for every node
      tip_states – dict[tip_name] = List[str], states for every tip

    return
      changes: dict[
        (parent_name, child_name) → {
          'mask':      List[bool],   # character change number
          'status':    List[None/str], # states number
          'type':      List[None/str]  # "apomorphy" or "homoplasy"
        }
      ]
    """
    # collect all tips in the tree
    all_tips = {t.name for t in tree.get_terminals()}

    # First step: collect all changes
    raw = {}
    counter = Counter()
    for parent in tree.find_clades(order="level"):
        # parent.name is the parent node, pN is its N
        pN = final[parent.name]
        for child in parent.clades:
            # child.name is the child node, cN is its N
            cN = final[child.name]
            # obtain the set of states in the subtree rooted at child
            subtree = {t.name for t in child.get_terminals()}
            mask = []
            status = []
            for i, (ps, cs) in enumerate(zip(pN, cN)):
                # if parent and child N have no intersection, record the change
                if not (ps & cs):
                    new = ''.join(sorted(cs))
                    mask.append(True)
                    # record the new state of child N
                    status.append(new)
                    # count the occurrence of (i, new) in the whole tree
                    counter[(i, new)] += 1
                else:
                    mask.append(False)
                    status.append(None)
            raw[(parent.name, child.name)] = {
                'mask': mask,
                'status': status,
                'subtree': subtree
            }

    # Second step: classify changes
    changes = {}
    # for every branch, check the status and classify
    for br, info in raw.items():
        types = []
        # count how many times each (index, new_state) appears outsides the subtree
        outside = all_tips - info['subtree']
        for i, new in enumerate(info['status']):
            if not info['mask'][i]:
                types.append(None)
                continue

            # if the states present in the tree more than once even it's not recorded as changes, then it's homoplasy (defined in WinClada)
            if counter[(i, new)] > 1:
                types.append("homoplasy")
            else:
                # if the state is unique, it's an apomorphy, or it's homoplasy
                found = any(tip_states[t][i] == new for t in outside)
                types.append("homoplasy" if found else "apomorphy")

        changes[br] = {
            'mask':   info['mask'],
            'status': info['status'],
            'type':   types
        }

    return changes

def detect_semiapomorphy(changes, final):
    """
    For the existing changes marked as 'homoplasy', further classify them
    If they have unique last-position (character) Fitch states, mark them
    as 'semi-apomorphy'.

    Logic
      1) for each (position, new_state) that appears multiple times in homoplasy branches,
         collect their last-position Fitch states.
      2) If these Fitch states are more than one, only mark those branches where 
         the last-position state occurs exactly once as 'semi-apomorphy'

    parameters:
      changes: dict[
        (parent_name, child_name) → {
          'mask':   List[bool],
          'status': List[str or None],
          'type':   List[str or None]  # 'apomorphy'|'homoplasy'
        }
      ]
      final: dict[node_name] = List[Set]  # collect the N for every node
    
    return:
      changes。
    """
    # last position is the last character in the final sequences
    sample = next(iter(final.values()))
    last_pos = len(sample) - 1

    # classify homoplasy branches by their (index, new_state) node, states
    groups = defaultdict(list)
    for branch, info in changes.items():
        for idx, t in enumerate(info['type']):
            if t == "homoplasy":
                new_state = info['status'][idx]
                groups[(idx, new_state)].append(branch)

    # for branches whose state occurs more than twice including twice, collect
    for (idx, new_state), branches in groups.items():
        if len(branches) < 2:
            continue

        # collect the states in the final character
        branch_last = {}
        for u, v in branches:
            cs = final[v][last_pos]
            s = ''.join(sorted(cs))
            branch_last[(u, v)] = s

        # count how many times each last character state appears
        counts = Counter(branch_last.values())
        if len(counts) <= 1:
            # if there is only one state, no semi-apomorphy
            continue

        # if the last character state occurs only once, mark it as semi-apomorphy
        for branch, last_state in branch_last.items():
            if counts[last_state] == 1:
                changes[branch]['type'][idx] = "semi-apomorphy"

    return changes

def drop_last_character(changes):
    """
    from the changes, remove the last character from all records.
    
    parameters:
      changes: dict[
        (parent_name, child_name) → {
          'mask':   List[bool],
          'status': List[str or None],
          'type':   List[str or None]
        }
      ]
    
    return:
      a new changes
    """
    trimmed = {}
    for branch, info in changes.items():
        trimmed[branch] = {
            'mask':   info['mask'][:-1],
            'status': info['status'][:-1],
            'type':   info['type'][:-1]
        }
    return trimmed
