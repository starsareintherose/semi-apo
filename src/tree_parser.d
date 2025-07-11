module tree_parser;

import std.stdio;
import std.string;
import std.conv;
import std.array;
import std.algorithm;
import std.ascii;
import std.file;

/**
 * Represents a node in a phylogenetic tree
 */
class TreeNode {
    string name;
    TreeNode[] children;
    TreeNode parent;
    double branchLength = 0.0;
    
    this(string name = null) {
        this.name = name;
    }
    
    /**
     * Check if this node is a terminal (leaf) node
     */
    bool isTerminal() {
        return children.length == 0;
    }
    
    /**
     * Get all terminal nodes in the subtree rooted at this node
     */
    TreeNode[] getTerminals() {
        TreeNode[] terminals;
        if (isTerminal()) {
            terminals ~= this;
        } else {
            foreach (child; children) {
                terminals ~= child.getTerminals();
            }
        }
        return terminals;
    }
    
    /**
     * Find all nodes in the tree using level-order traversal
     */
    TreeNode[] findClades(string order = "level") {
        TreeNode[] result;
        TreeNode[] queue = [this];
        
        while (queue.length > 0) {
            TreeNode current = queue[0];
            queue = queue[1..$];
            result ~= current;
            
            foreach (child; current.children) {
                queue ~= child;
            }
        }
        
        return result;
    }
}

/**
 * Represents a phylogenetic tree
 */
class Tree {
    TreeNode root;
    
    this(TreeNode root) {
        this.root = root;
    }
    
    /**
     * Find all nodes in the tree
     */
    TreeNode[] findClades(string order = "level") {
        return root.findClades(order);
    }
    
    /**
     * Get all terminal nodes in the tree
     */
    TreeNode[] getTerminals() {
        return root.getTerminals();
    }
}

/**
 * Simple Newick format parser
 */
class NewickParser {
    private string input;
    private size_t pos;
    private int nodeCounter;
    
    this(string newick) {
        this.input = newick.strip();
        this.pos = 0;
        this.nodeCounter = 0;
    }
    
    Tree parse() {
        TreeNode root = parseNode();
        ensureNodesHaveNames(root);
        return new Tree(root);
    }
    
    private TreeNode parseNode() {
        skipWhitespace();
        TreeNode node = new TreeNode();
        
        if (pos < input.length && input[pos] == '(') {
            // Internal node with children
            pos++; // skip '('
            skipWhitespace();
            
            // Parse first child
            node.children ~= parseNode();
            node.children[$-1].parent = node;
            
            // Parse additional children
            while (pos < input.length && input[pos] == ',') {
                pos++; // skip ','
                skipWhitespace();
                node.children ~= parseNode();
                node.children[$-1].parent = node;
            }
            
            skipWhitespace();
            if (pos < input.length && input[pos] == ')') {
                pos++; // skip ')'
            }
            skipWhitespace();
        }
        
        // Parse node name if present
        if (pos < input.length && (isAlphaNum(input[pos]) || input[pos] == '_')) {
            string name = "";
            while (pos < input.length && (isAlphaNum(input[pos]) || input[pos] == '_')) {
                name ~= input[pos];
                pos++;
            }
            node.name = name;
        }
        
        // Parse branch length if present
        if (pos < input.length && input[pos] == ':') {
            pos++; // skip ':'
            string lengthStr = "";
            while (pos < input.length && (isDigit(input[pos]) || input[pos] == '.')) {
                lengthStr ~= input[pos];
                pos++;
            }
            if (lengthStr.length > 0) {
                try {
                    node.branchLength = to!double(lengthStr);
                } catch (ConvException) {
                    node.branchLength = 0.0;
                }
            }
        }
        
        return node;
    }
    
    private void skipWhitespace() {
        while (pos < input.length && isWhite(input[pos])) {
            pos++;
        }
    }
    
    private void ensureNodesHaveNames(TreeNode node) {
        if (node.name is null || node.name.length == 0) {
            node.name = "Node" ~ to!string(nodeCounter++);
        }
        
        foreach (child; node.children) {
            ensureNodesHaveNames(child);
        }
    }
}

/**
 * Load a phylogenetic tree from a file in Newick format
 */
Tree loadTree(string path) {
    string content = std.file.readText(path).strip();
    
    // Remove trailing semicolon if present
    if (content.endsWith(";")) {
        content = content[0..$-1];
    }
    
    NewickParser parser = new NewickParser(content);
    return parser.parse();
}