#!/usr/bin/env python3

def test_line(line):
    """Returns true for valid lines. Not comments or blank lines."""
    if not line.strip():
        return False  # if the line is blank
    if line.startswith("#"):
        return False  # comment line
    return True

def parse_gff_gft(file_path):
    """
    Parse a GFF/GTF file and return a list of features.

    Parameters:
    file_path (str): Path to the GFF or GTF file.

    Returns:
    list: A list of tuples, each containing the fields of a feature.
    """
    features = []
    with open(file_path, 'r', encoding='utf-8', errors='ignore') as file:
        for line in file:
            if not test_line(line):
                continue
            
            # Split the line into fields
            fields = line.strip().split('\t')
            
            # Extract necessary fields
            seqname = fields[0]
            source = fields[1]
            feature = fields[2]
            start = int(fields[3])
            end = int(fields[4])
            score = fields[5]
            strand = fields[6]
            frame = fields[7]
            attribute = fields[8].rstrip(';')
            
            # Append feature information as a tuple to the list
            features.append((seqname, source, feature, start, end, score, strand, frame, attribute))
    
    return features
