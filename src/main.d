module main;

import std.stdio;
import std.getopt;
import std.conv;
import std.array;
import std.algorithm;
import std.string;
import std.process;
import tree_parser;
import matrix_parser;
import fitch;

void main(string[] args) {
    string treeFile;
    string matrixFile;
    string outText;
    string outSvg;
    string outOriSvg;
    
    auto helpInformation = getopt(args,
        "tree", "Newick tree file", &treeFile,
        "matrix", "TNT matrix file", &matrixFile,
        "out-text", "output text file", &outText,
    );
    
    if (helpInformation.helpWanted || treeFile.length == 0 || matrixFile.length == 0 || 
        outText.length == 0) {
        defaultGetoptPrinter("Semi-apomorphy detection in phylogenetic trees", helpInformation.options);
        return;
    }
    
    // Load tree and matrix data
    Tree tree = loadTree(treeFile);
    string[][string] tipStates = parseTntMatrix(matrixFile);
    
    // Perform the analysis pipeline
    auto downResult = downPass(tree, tipStates);
    
    // Extract fullSet for upPass (get all unique states)
    string[] allStates;
    foreach (states; tipStates.values) {
        foreach (c; states) {
            if (c != "?" && !allStates.canFind(c)) {
                allStates ~= c;
            }
        }
    }
    
    int k = cast(int)tipStates.values.front.length;
    string[][] fullSet;
    foreach (i; 0..k) {
        fullSet ~= allStates.dup;
    }
    
    auto final_ = upPass(tree, downResult, fullSet);
    auto changes = detectChanges(tree, final_, tipStates);
    changes = detectSemiapomorphy(changes, final_);
    changes = dropLastCharacter(changes);
    
    // Write output to text file
    auto outFile = File(outText, "w");
    scope(exit) outFile.close();
    
    // Progress report
    outFile.writeln("Node\t↓n\t↓B\t↑N");
    foreach (clade; tree.findClades()) {
        string[] nStrs, BStrs, NStrs;
        
        if (clade.name in downResult) {
            foreach (nSet; downResult[clade.name].n) {
                nStrs ~= nSet.dup.sort.array.join("");
            }
            foreach (BSet; downResult[clade.name].B) {
                BStrs ~= BSet.dup.sort.array.join("");
            }
        }
        
        if (clade.name in final_) {
            foreach (NSet; final_[clade.name]) {
                NStrs ~= NSet.dup.sort.array.join("");
            }
        }
        
        outFile.writefln("%s\t%s\t%s\t%s", 
                        clade.name, 
                        nStrs.join(","), 
                        BStrs.join(","), 
                        NStrs.join(","));
    }
    
    outFile.writeln("\n# Character changes (only list changed characters)");
    foreach (branchKey, info; changes) {
        // Filter out unchanged characters
        int[] changed;
        foreach (i, mask; info.mask) {
            if (mask) changed ~= cast(int)i;
        }
        
        if (changed.length == 0) continue;
        
        string posStr = changed.map!(to!string).join(",");
        string[] statusStrs;
        foreach (i; changed) {
            statusStrs ~= info.status[i];
        }
        string statusStr = statusStrs.join(",");
        
        outFile.writefln("%s:\tcharacter[%s]\tstates[%s]", branchKey, posStr, statusStr);
    }
    
    // Write apomorphies, semi-apomorphies and homoplasies
    writeChangesByType(outFile, changes, "apomorphy", "# Apomorphies");
    writeChangesByType(outFile, changes, "semi-apomorphy", "# Semi-apomorphies");
    writeChangesByType(outFile, changes, "homoplasy", "# Homoplasy");
    
}

void writeChangesByType(File outFile, ChangeInfo[string] changes, string targetType, string header) {
    outFile.writeln("\n" ~ header);
    
    foreach (branchKey, info; changes) {
        int[] indices;
        foreach (i, typeStr; info.type) {
            if (typeStr == targetType) {
                indices ~= cast(int)i;
            }
        }
        
        if (indices.length == 0) continue;
        
        string posStr = indices.map!(to!string).join(",");
        string[] stateStrs;
        foreach (i; indices) {
            stateStrs ~= info.status[i];
        }
        string stateStr = stateStrs.join(",");
        
        outFile.writefln("%s:\tcharacter[%s]\tstates[%s]", branchKey, posStr, stateStr);
    }
}


