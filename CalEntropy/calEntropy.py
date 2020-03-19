from CalEntropy.Atoms import *

# (cal/mol-1, 300K)
__Amino_acid_conformational_entropy = {"LYS": -189, "ARG": -188, "GLN": -173, "MET": -146, "GLU": -146,
                                           "ILE": -76, "LEU": -71, "ASN": -103, "THR": -108, "VAL": -43,
                                           "TYR": -113, "SER": -111,  "HIS": -95, "ASP": -78, "CYS": -85,
                                           "TRP": -99, "PHE": -62, "ALA": -0, "PRO": -6, "GLY": 0}


def get_pdb_atoms(pdbfile):
    atoms = []
    with open(pdbfile, 'r') as pdb:
        lines = tuple(atom.strip() for atom in pdb.readlines()[0: -2])
    for atom in lines:
        identifier = atom[0: 6].strip()
        if identifier == 'ATOM':
            atom_type = atom[12: 16].strip()
            residue = atom[17: 20].strip()
            index = atom[22: 26].strip()
            atoms.append(ResidueAtom(atom_type=atom_type,
                                     coordinate=tuple(map(float, (atom[30: 38].strip(),
                                                                  atom[38: 46].strip(),
                                                                  atom[46: 54].strip()))),
                                     index=index, residue=residue))
    return tuple(atoms)


def get_sdf_atoms(sdffile, single=True):
    atom_types = []
    atoms = []
    if single:
        with open(sdffile, 'r') as sdf:
            lines = tuple(sdf.readlines())
    else:
        lines = sdffile
    split_line = lines[3].strip().split()
    natoms = int(split_line[0])
    nbonds = int(split_line[1])
    if natoms > len(lines):
        count = len(split_line[0])
        while natoms > len(lines):
            natoms = int(split_line[0][:count])
            count -= 1
        nbonds = int(split_line[0][count + 1:])
    atom_lines = lines[4: 4 + natoms]
    property_lines = lines[5 + natoms + nbonds:]

    for i, line in enumerate(property_lines):
        index = i
        while line == ">  <TYPE_INFO>\n":
            next_line = property_lines[index + 1]
            if next_line == "\n":
                break
            atom_types.extend(next_line.strip().split())
            index += 1

    for i, atom in enumerate(atom_lines):
        split_atom = atom.strip().split()
        atom_type = atom_types[i]
        coordinate = tuple(map(float, split_atom[0: 3]))
        atoms.append(LigandAtom(atom_type=atom_type, coordinate=coordinate))

    return tuple(atoms)


def cal_lw_fs(ligand, protein):
    ligand_protein_contacts = {}
    total_flexible_entropy = 0
    total_protein_water_loss = 0
    total_ligand_water_loss = 0
    deduct_water = 1
    protein_contacted_water = {}
    contacted_atom = set()

    for ligand_atom in ligand:
        ligand_atom_coordinate = ligand_atom.get_coordinate()
        ligand_protein_contacts[ligand_atom] = set()
        for residue_atom in protein:
            residue_atom_coordinate = residue_atom.get_coordinate()
            dist = ligand_atom_coordinate.cal_dist(residue_atom_coordinate)
            if dist > 3.5:
                pass
            else:
                ligand_protein_contacts[ligand_atom].add(residue_atom)

    for latm in ligand_protein_contacts:

        # calculate ligand water loss
        if len(ligand_protein_contacts[latm]) >= latm.get_water_count():
            ligand_water_loss = latm.get_water_count()
            total_ligand_water_loss += ligand_water_loss
        else:
            total_ligand_water_loss += len(ligand_protein_contacts[latm])

        for atm in ligand_protein_contacts[latm]:
            # calculate flexibility entropy
            flexible_entropy = 0
            index, residue = atm.get_residue()
            if index not in contacted_atom:
                contacted_atom.add(index)
                flexible_entropy = __Amino_acid_conformational_entropy[residue]
            total_flexible_entropy += flexible_entropy

            # calculate protein water loss
            if atm not in protein_contacted_water:
                protein_contacted_water[atm] = atm.get_water_count()
            original_water = protein_contacted_water[atm]
            if original_water == 0:
                pass
            else:
                left_water = original_water - deduct_water
                total_protein_water_loss += deduct_water
                protein_contacted_water.update({atm: left_water})
    total_water_loss = total_ligand_water_loss + total_protein_water_loss
    water_entropy = total_water_loss * 10 * 300 * 0.001
    return water_entropy, total_flexible_entropy
