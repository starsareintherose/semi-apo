module fitch;

import std.stdio;
import std.string;
import std.conv;
import std.array;
import std.algorithm;
import std.range;
import tree_parser;

/**
 * Efficient character state representation using bit masks
 */
struct StateMask {
    uint mask;  // Bit mask for states (supports up to 32 states)
    
    this(uint mask) {
        this.mask = mask;
    }
    
    this(string[] states, string[] allStates) {
        this.mask = 0;
        foreach (state; states) {
            foreach (i, s; allStates) {
                if (s == state) {
                    this.mask |= (1 << i);
                    break;
                }
            }
        }
    }
    
    string[] toStringArray(string[] allStates) {
        string[] result;
        foreach (i, state; allStates) {
            if (mask & (1 << i)) {
                result ~= state;
            }
        }
        return result;
    }
    
    StateMask intersect(StateMask other) {
        return StateMask(this.mask & other.mask);
    }
    
    StateMask union_(StateMask other) {
        return StateMask(this.mask | other.mask);
    }
    
    bool isEmpty() {
        return mask == 0;
    }
    
    bool isSubsetOf(StateMask other) {
        return (this.mask & other.mask) == this.mask;
    }
    
    uint popcount() {
        uint count = 0;
        uint m = mask;
        while (m) {
            count += m & 1;
            m >>= 1;
        }
        return count;
    }
}

/**
 * Data structure for down pass results using efficient masks
 */
struct MaskDownPassResult {
    StateMask[] n;  // Character state masks
    StateMask[] B;  // Character state masks
}

/**
 * Perform down pass from tips to root using mask optimization
 */
MaskDownPassResult[string] downPassOptimized(Tree tree, string[][string] tipStates) {
    string[] allStates;
    
    // Collect all states from tips
    foreach (states; tipStates.values) {
        foreach (c; states) {
            // Only add non-unknown states
            if (c != "?") {
                if (!allStates.canFind(c)) {
                    allStates ~= c;
                }
            }
        }
    }
    
    // Get sequence length
    int k = cast(int)tipStates.values.front.length;
    
    // Create full set masks for each position
    StateMask fullMask = StateMask(0);
    foreach (i; 0..allStates.length) {
        fullMask.mask |= (1 << i);
    }
    
    MaskDownPassResult[string] result;
    
    void recurseDown(TreeNode clade) {
        // If clade is terminal, initialize n and B
        if (clade.isTerminal()) {
            if (clade.name !in tipStates) {
                throw new Exception("Terminal node '" ~ clade.name ~ "' not found in tip states");
            }
            string[] seq = tipStates[clade.name];
            StateMask[] nList, BList;
            
            foreach (i, c; seq) {
                if (c == "?") {
                    // Unknown state, use full set
                    nList ~= fullMask;
                    BList ~= fullMask;
                } else {
                    StateMask singleState = StateMask([c], allStates);
                    nList ~= singleState;
                    BList ~= singleState;
                }
            }
            result[clade.name] = MaskDownPassResult(nList, BList);
            return;
        }
        
        // Process children first
        foreach (child; clade.children) {
            recurseDown(child);
        }
        
        // Combine results from children (assuming binary tree)
        StateMask[] L = result[clade.children[0].name].n;
        StateMask[] R = result[clade.children[1].name].n;
        StateMask[] n, B;
        
        foreach (i; 0..k) {
            // Get intersection and union
            StateMask intersection = L[i].intersect(R[i]);
            StateMask unionSet = L[i].union_(R[i]);
            
            if (!intersection.isEmpty()) {
                n ~= intersection;
                B ~= unionSet;
            } else {
                n ~= unionSet;
                B ~= fullMask;
            }
        }
        
        result[clade.name] = MaskDownPassResult(n, B);
    }
    
    recurseDown(tree.root);
    return result;
}

/**
 * Data structure for down pass results
 */
struct DownPassResult {
    string[][] n;  // Sets represented as sorted string arrays
    string[][] B;  // Sets represented as sorted string arrays
}

/**
 * Perform down pass from tips to root
 */
DownPassResult[string] downPass(Tree tree, string[][string] tipStates) {
    // Use optimized mask-based implementation internally
    auto maskResult = downPassOptimized(tree, tipStates);
    
    // Convert back to string arrays for compatibility
    DownPassResult[string] result;
    string[] allStates;
    
    // Collect all states
    foreach (states; tipStates.values) {
        foreach (c; states) {
            if (c != "?" && !allStates.canFind(c)) {
                allStates ~= c;
            }
        }
    }
    
    foreach (nodeName, maskData; maskResult) {
        string[][] nArrays, BArrays;
        
        foreach (nMask; maskData.n) {
            nArrays ~= nMask.toStringArray(allStates);
        }
        foreach (BMask; maskData.B) {
            BArrays ~= BMask.toStringArray(allStates);
        }
        
        result[nodeName] = DownPassResult(nArrays, BArrays);
    }
    
    return result;
}

/**
 * Perform up pass from root to tips
 * Returns: mapping from node name to array of string arrays (each representing a set)
 */
string[][][string] upPass(Tree tree, DownPassResult[string] down, string[][] fullSets) {
    string[][][string] final_;
    int k = cast(int)fullSets.length;
    
    void recurseUp(TreeNode clade, string[][] parentN = null) {
        string[][] downN = down[clade.name].n;
        string[][] downB = down[clade.name].B;
        string[][] N;
        
        // If parentN is null, this is the root node, so N = n
        if (parentN is null) {
            foreach (set; downN) {
                N ~= set.dup;
            }
        } else {
            // If clade is terminal, use n directly
            if (clade.isTerminal()) {
                foreach (set; downN) {
                    N ~= set.dup;
                }
            } else {
                // Check if parentN is subset of downN
                foreach (i; 0..k) {
                    string[] Ni;
                    if (isSubset(parentN[i], downN[i])) {
                        Ni = parentN[i].dup;
                    } else {
                        string[] intersection = setIntersection(parentN[i], downB[i]);
                        Ni = setUnion(downN[i], intersection);
                    }
                    N ~= Ni;
                }
            }
        }
        
        final_[clade.name] = N;
        
        foreach (child; clade.children) {
            recurseUp(child, N);
        }
    }
    
    recurseUp(tree.root, null);
    return final_;
}

/**
 * Detect changes between parent and child nodes
 */
struct ChangeInfo {
    bool[] mask;
    string[] status;
    string[] type;
    string[] subtreeTaxa;
}

ChangeInfo[string] detectChanges(Tree tree, string[][][string] final_, string[][string] tipStates) {
    // Collect all tips in the tree
    string[] allTips;
    foreach (tip; tree.getTerminals()) {
        allTips ~= tip.name;
    }
    
    // Count occurrences of state changes
    int[string] counter;
    ChangeInfo[string] rawChanges;
    
    // First pass: collect all changes
    foreach (parent; tree.findClades()) {
        if (parent.name !in final_) continue;
        string[][] pN = final_[parent.name];
        
        foreach (child; parent.children) {
            if (child.name !in final_) continue;
            string[][] cN = final_[child.name];
            string branchKey = parent.name ~ "->" ~ child.name;
            
            // Get subtree taxa
            string[] subtree;
            foreach (tip; child.getTerminals()) {
                subtree ~= tip.name;
            }
            
            bool[] mask;
            string[] status;
            
            foreach (i; 0..pN.length) {
                string[] intersection = setIntersection(pN[i], cN[i]);
                if (intersection.length == 0) {
                    // Change detected
                    string[] sortedChars = cN[i].dup.sort.array;
                    string newState = sortedChars.join("");
                    mask ~= true;
                    status ~= newState;
                    
                    // Count occurrences
                    string key = to!string(i) ~ ":" ~ newState;
                    counter[key]++;
                } else {
                    mask ~= false;
                    status ~= "";
                }
            }
            
            rawChanges[branchKey] = ChangeInfo(mask, status, [], subtree);
        }
    }
    
    // Second pass: classify changes as apomorphy/homoplasy
    ChangeInfo[string] changes;
    foreach (branchKey, info; rawChanges) {
        string[] types;
        string[] outside = setDifference(allTips, info.subtreeTaxa);
        
        foreach (i, newState; info.status) {
            if (!info.mask[i]) {
                types ~= "";
                continue;
            }
            
            string key = to!string(i) ~ ":" ~ newState;
            if (counter[key] > 1) {
                types ~= "homoplasy";
            } else {
                // Check if state appears in outside taxa
                bool found = false;
                foreach (taxon; outside) {
                    if (i < tipStates[taxon].length && tipStates[taxon][i] == newState) {
                        found = true;
                        break;
                    }
                }
                types ~= found ? "homoplasy" : "apomorphy";
            }
        }
        
        changes[branchKey] = ChangeInfo(info.mask, info.status, types, info.subtreeTaxa);
    }
    
    return changes;
}

/**
 * Detect semi-apomorphies by analyzing last position states
 */
ChangeInfo[string] detectSemiapomorphy(ChangeInfo[string] changes, string[][][string] final_) {
    // Get last position index
    auto sample = final_.values.front;
    int lastPos = cast(int)sample.length - 1;
    
    // Group homoplasy branches by (index, new_state)
    string[][string] groups;
    
    foreach (branchKey, info; changes) {
        foreach (idx, typeStr; info.type) {
            if (typeStr == "homoplasy") {
                string newState = info.status[idx];
                string groupKey = to!string(idx) ~ ":" ~ newState;
                groups[groupKey] ~= branchKey;
            }
        }
    }
    
    // Process groups with multiple branches
    foreach (groupKey, branches; groups) {
        if (branches.length < 2) continue;
        
        // Collect last position states for these branches
        string[string] branchLastStates;
        foreach (branchKey; branches) {
            auto parts = branchKey.split("->");
            string childName = parts[1];
            string[] sortedChars = final_[childName][lastPos].dup.sort.array;
            string lastState = sortedChars.join("");
            branchLastStates[branchKey] = lastState;
        }
        
        // Count occurrences of each last state
        int[string] counts;
        foreach (lastState; branchLastStates.values) {
            counts[lastState]++;
        }
        
        if (counts.length <= 1) continue;
        
        // Mark branches with unique last states as semi-apomorphy
        foreach (branchKey, lastState; branchLastStates) {
            if (counts[lastState] == 1) {
                auto parts = groupKey.split(":");
                int idx = to!int(parts[0]);
                changes[branchKey].type[idx] = "semi-apomorphy";
            }
        }
    }
    
    return changes;
}

/**
 * Drop last character from all changes
 */
ChangeInfo[string] dropLastCharacter(ChangeInfo[string] changes) {
    ChangeInfo[string] trimmed;
    
    foreach (branchKey, info; changes) {
        ChangeInfo newInfo;
        newInfo.mask = info.mask[0..$-1];
        newInfo.status = info.status[0..$-1];
        newInfo.type = info.type[0..$-1];
        newInfo.subtreeTaxa = info.subtreeTaxa;
        trimmed[branchKey] = newInfo;
    }
    
    return trimmed;
}

/**
 * Helper functions for set operations
 */
string[] setIntersection(string[] a, string[] b) {
    string[] result;
    foreach (item; a) {
        if (b.canFind(item) && !result.canFind(item)) {
            result ~= item;
        }
    }
    return result.sort.array;
}

string[] setUnion(string[] a, string[] b) {
    string[] result = a.dup;
    foreach (item; b) {
        if (!result.canFind(item)) {
            result ~= item;
        }
    }
    return result.sort.array;
}

string[] setDifference(string[] a, string[] b) {
    string[] result;
    foreach (item; a) {
        if (!b.canFind(item)) {
            result ~= item;
        }
    }
    return result;
}

bool isSubset(string[] subset, string[] superset) {
    foreach (item; subset) {
        if (!superset.canFind(item)) {
            return false;
        }
    }
    return true;
}