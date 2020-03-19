#!/usr/bin/python
import os
import csv
import multiprocessing as mp
import time
import math


class Atom(object):
    def __init__(self, atom_type, coordinate):
        self.__type = atom_type
        self.__coordinate = Point(coordinate)

    def get_coordinate(self):
        return self.__coordinate

    def get_type(self):
        return self.__type


class Point(object):
    def __init__(self, coordinate):
        self.x = coordinate[0]
        self.y = coordinate[1]
        self.z = coordinate[2]

    def cal_dist(self, point):
        __distance = math.sqrt((self.x - point.x) ** 2 + (self.y - point.y) ** 2 + (self.z - point.z) ** 2)
        return __distance


class ResidueAtom(Atom):
    def __init__(self, atom_type, coordinate, index, residue):
        __residue_number_of_water = {"ARG": {"NE": 1, "NH1": 2, "NH2": 2},
                                     "ASN": {"ND2": 2, "OD1": 2},
                                     "ASP": {"OD1": 2, "OD2": 2},
                                     "GLN": {"NE2": 2, "OE1": 2},
                                     "GLU": {"OE1": 2, "OE2": 2},
                                     "HIS": {"ND1": 1, "NE2": 1},
                                     "LYS": {"NZ": 3},
                                     "SER": {"OG": 2},
                                     "THR": {"OG1": 2},
                                     "TRP": {"NE1": 1},
                                     "TYR": {"OH": 1}}
        Atom.__init__(self, atom_type, coordinate)
        self.__residue = residue
        self.__index = index
        try:
            self.__number_of_water = __residue_number_of_water[residue][atom_type]
        except KeyError:
            self.__number_of_water = 0

    def get_residue(self):
        return self.__index, self.__residue

    def get_water_count(self):
        return self.__number_of_water


class LigandAtom(Atom):
    def __init__(self, atom_type, coordinate):
        # LIDAEUS default setting
        __ligand_number_of_water = {"N.4": 1, "N.3": 2, "N.2": 3, "N.1": 4, "N.ar": 1, "N.pl3": 1, "N.am": 1,
                                    "O.3": 1, "O.2": 2, "O.co2": 2}
        Atom.__init__(self, atom_type, coordinate)
        try:
            self.__number_of_water = __ligand_number_of_water[atom_type]
        except KeyError:
            self.__number_of_water = 0

    def get_water_count(self):
        return self.__number_of_water


# (cal/mol-1, 300K)
__Amino_acid_conformational_probability = {"LYS": -189, "ARG": -188, "GLN": -173, "MET": -146, "GLU": -146,
                                           "ILE": -76, "LEU": -71, "ASN": -103, "THR": -108, "VAL": -43,
                                           "TYR": -113, "SER": -111, "HIS": -95, "ASP": -78, "CYS": -85,
                                           "TRP": -99, "PHE": -62, "ALA": 0, "PRO": -6, "GLY": 0}


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
                                     coordinate=tuple(map(float, (
                                         atom[30: 38].strip(),
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
        dists = set()
        ligand_atom_coordinate = ligand_atom.get_coordinate()
        for residue_atom in protein:
            residue_atom_coordinate = residue_atom.get_coordinate()
            dist = ligand_atom_coordinate.cal_dist(residue_atom_coordinate)
            if dist > 3.5:
                pass
            else:
                dists.add(residue_atom)
        ligand_protein_contacts[ligand_atom] = dists

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
                flexible_entropy = __Amino_acid_conformational_probability[residue]
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


def multi_cal_entropy(job, patms):
    entropies = []
    for mol in job:
        try:
            latms = get_sdf_atoms(mol[1], single=False)
            entropies.append((mol[0], cal_lw_fs(latms, patms)))
        except:
            print(mol[0])
    return tuple(entropies)


def multi_collect(entropies):
    global Multi_results
    Multi_results.extend(entropies)


def split_jobs(jobs, part=3):
    ratio = float(1) / part
    job = []
    for each_part in range(0, part):
        job.append(
            jobs[int(math.ceil(ratio * each_part * len(jobs))): int(math.ceil(ratio * (each_part + 1) * len(jobs)))])
    return tuple(job)


if __name__ == "__main__":
    path = "/datastore/home/s1829725/project/Targets/"
    targets = os.listdir(path)
    protein = "/protein_nowater.pdb"
    vs_results = "/outputmod.sdf"
    output = "/datastore/home/s1829725/project/CalEntropy/water_side_chain/"
    outputfiles = os.listdir(output)

    for target in targets:
        if target + ".csv" in outputfiles:
            print(target + ".csv" + " done")
            pass
        else:
            Multi_results = []
            multi_sdf = []
            protein_atoms = get_pdb_atoms(path + target + protein)
            with open(path + target + vs_results, "r") as SDF:
                sdfs = SDF.readlines()
            sdffile = []
            info = []
            for i, line in enumerate(sdfs):
                sdffile.append(line)
                if line == "> <Name>\n":
                    name = sdfs[i + 1].strip()
                    info.append(name)
                if line == ">  <SCORE_INFO>\n":
                    original_energies = sdfs[i + 1].strip().split()[0: 4]
                    info.extend(original_energies)
                if line == "$$$$\n":
                    multi_sdf.append((tuple(info), sdffile))
                    sdffile = []
                    info = []
            mols = multi_sdf
            # split jobs
            start = time.time()
            p = mp.Pool(mp.cpu_count())
            jobs = split_jobs(mols, part=mp.cpu_count())
            for job in jobs:
                p.apply_async(func=multi_cal_entropy, args=(job, protein_atoms), callback=multi_collect)
            p.close()
            p.join()
            print("%s %s" % (target, time.time() - start))
            field_names = ("Name", "Enthalpy", "vDW", "HBD", "HBA",
                           "Lost Water entropy", "Side chain conformational entropy")
            with open(output + target + ".csv", "wb") as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow(field_names)
                for entropy in Multi_results:
                    writer.writerow((entropy[0][0], entropy[0][1], entropy[0][2], entropy[0][3], entropy[0][4],
                                     entropy[1][0], entropy[1][1] * 0.01))
            # debug
            """for mol in mols:
                latms = get_sdf_atoms(mol[1], single=False)
                _, __ = cal_lw_fs(latms, protein_atoms)
                print("%s %s" % (_, __))"""