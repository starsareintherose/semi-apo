import re

def parse_tnt_matrix(path):
    """
    parser the tnt matrix file.
    1. Skip comments starting with ' and any non-numeric leading lines.
    2. The first numeric line is parsed as nchar (number of characters) and ntaxa (number of taxa).
    3. Read the following taxa lines, stopping at lines starting with ';' or when ntaxa records are read.
    4. Each line's first column is the taxon, the second column is the status sequence.
    5. Finally, check if the number of taxa read and the length of each sequence match the header; if not, print WARNING.
    """
    data = {}
    nchar = None
    ntaxa = None
    have_header = False

    with open(path) as f:
        for raw in f:
            line = raw.strip()
            # skip empty lines and comments
            if not line or line.startswith("'"):
                continue

            # numeric header line
            if not have_header:
                m = re.fullmatch(r'(\d+)\s+(\d+)', line)
                if m:
                    nchar, ntaxa = map(int, m.groups())
                    have_header = True
                continue

            # if there is a semicolon, it indicates the end of the matrix
            if line.startswith(';'):
                break

            # taxa line
            parts = re.split(r'\s+', line)
            if len(parts) >= 2:
                taxon, seq = parts[0], parts[1]
                data[taxon] = list(seq)
                # if ntaxa is specified, stop reading when we reach it
                if ntaxa is not None and len(data) >= ntaxa:
                    break

    # check if the number of taxa and sequence length match the header
    if ntaxa is not None and len(data) != ntaxa:
        print(f"WARNING: the number of taxa parsed ({len(data)}) does not match the header ({ntaxa}).")
    if nchar is not None:
        for taxon, seq in data.items():
            if len(seq) != nchar:
                print(f"WARNING: the sequence length for taxon {taxon} ({len(seq)}) does not match the declared length ({nchar}).")

    return data
