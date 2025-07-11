module matrix_parser;

import std.stdio;
import std.string;
import std.conv;
import std.regex;
import std.array;

/**
 * Parse TNT matrix file
 * 
 * 1. Skip comments starting with ' and any non-numeric leading lines.
 * 2. The first numeric line is parsed as nchar (number of characters) and ntaxa (number of taxa).
 * 3. Read the following taxa lines, stopping at lines starting with ';' or when ntaxa records are read.
 * 4. Each line's first column is the taxon, the second column is the status sequence.
 * 5. Finally, check if the number of taxa read and the length of each sequence match the header; if not, print WARNING.
 */
string[][string] parseTntMatrix(string path) {
    string[][string] data;
    int nchar = -1;
    int ntaxa = -1;
    bool haveHeader = false;
    
    auto file = File(path, "r");
    scope(exit) file.close();
    
    foreach (rawLine; file.byLine()) {
        string line = rawLine.idup.strip();
        // Skip empty lines and comments
        if (line.length == 0 || line.startsWith("'")) {
            continue;
        }
        // Numeric header line
        if (!haveHeader) {
            auto headerMatch = matchFirst(line, regex(r"^(\d+)\s+(\d+)$"));
            if (headerMatch) {
                nchar = to!int(headerMatch[1]);
                ntaxa = to!int(headerMatch[2]);
                haveHeader = true;
                continue;
            }
        }
        
        // If there is a semicolon, it indicates the end of the matrix
        if (line.startsWith(";")) {
            break;
        }
        
        // Taxa line
        auto parts = line.split();
        if (parts.length >= 2) {
            string taxon = parts[0];
            string seq = parts[1];
            string[] charArray;
            foreach (c; seq) {
                charArray ~= [c];
            }
            data[taxon] = charArray;
            
            // If ntaxa is specified, stop reading when we reach it
            if (ntaxa != -1 && data.length >= ntaxa) {
                break;
            }
        }
    }
    
    // Check if the number of taxa and sequence length match the header
    if (ntaxa != -1 && data.length != ntaxa) {
        writefln("WARNING: the number of taxa parsed (%d) does not match the header (%d).", 
                 data.length, ntaxa);
    }
    
    if (nchar != -1) {
        foreach (taxon, seq; data) {
            if (seq.length != nchar) {
                writefln("WARNING: the sequence length for taxon %s (%d) does not match the declared length (%d).", 
                         taxon, seq.length, nchar);
            }
        }
    }
    
    return data;
}
