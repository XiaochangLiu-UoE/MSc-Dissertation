import csv
Result = "Z:\\project\\residue.txt"
hydrophobic_aliphatic = {'ALA', 'VAL', 'ILE', 'LEU', 'MET', 'GLY', 'PRO', 'CYS'}
hydrophobic_aromatic = {'PHE', 'TYR', 'TRP'}
positive = {'ARG', 'HIS', 'LYS'}
negative = {'ASP', 'GLU'}
polar = {'SER', 'THR', 'ASN', 'GLN'}
metal = {'ZN', 'MG', 'FE', 'CU', 'AL'}
hdw = {'ALA', 'VAL', 'LEU', 'ILE', 'PRO', 'PHE', 'CYS', 'MET', 'GLY'}
hbd = {'ARG', 'ASN', 'GLN', 'HIS', 'LYS', 'SER', 'THR', 'TRP', 'TYR'}
hba = {'ASN', 'ASP', 'GLN', 'GLU', 'HIS', 'SER', 'THR', 'TYR'}
hydrophobicity = {"ALA": 0.31, "ARG": -1.01, "ASN": -0.60, "ASP": -0.77, "CYS": 1.54, "GLN": -0.22, "GLU": -0.64,
                  "HIS": 0.13, "ILE": 1.80, "LEU": 1.70, "LYS": -0.99, "MET": 1.23, "PHE": 1.79, "PRO": 0.72,
                  "SER": -0.04, "THR": 0.26, "TRP": 2.25, "TYR": 0.96, "VAL": 1.22, "GLY": 0}
caldict = {}
with open(Result, "r") as o:
    for line in o:
        line = line.rstrip("\n")
        words = line.split()
        protein = words[0]
        if protein not in caldict:
            caldict[protein] = list(0 for i in range(11))
        AA = words[1].split('-')[0]
        if AA in hydrophobic_aliphatic:
            caldict[protein][0] += 1
        elif AA in hydrophobic_aromatic:
            caldict[protein][1] += 1
        elif AA in polar:
            caldict[protein][2] += 1
        elif AA in positive:
            caldict[protein][3] += 1
        elif AA in negative:
            caldict[protein][4] += 1
        elif AA in metal:
            caldict[protein][5] += 1
        else:
            caldict[protein][6] += 1
        if AA in hdw:
            caldict[protein][7] += 1
        if AA in hbd:
            caldict[protein][8] += 1
        if AA in hba:
            caldict[protein][9] += 1
        try:
            caldict[protein][10] += hydrophobicity[AA]
        except KeyError:
            caldict[protein][10] += 0

with open('Z:\\project\\Metadata\\BindingSitesAnalysis.csv', 'w+', newline='') as w:
    fieldnames = ('Targets', 'Hydrophobic_aliphatic', 'Hydrophobic_aromatic', 'Polar Uncharged', 'Positive', 'Negative',
                  'Metal', 'Other Cofactors',  'HDW', 'HBD', 'HBA', 'Hydrophobicity')
    writer = csv.writer(w)
    writer.writerow(fieldnames)
    for protein in caldict:
        row = [protein]
        row.extend(caldict[protein])
        writer.writerow(row)
