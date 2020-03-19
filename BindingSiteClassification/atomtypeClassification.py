import os
import csv
import numpy as np
import plotly.figure_factory as ff
from plotly.offline import plot
from BindingSiteClassification.Descriptor import *


class AtomTypeDescriptor(object):
    def __init__(self, carbons, oxygens, nitrogens, target):
        coordinate = []
        self.__target = target
        self.__carbon = Descriptor(carbons)
        self.__oxygen = Descriptor(oxygens)
        self.__nitrogen = Descriptor(nitrogens)
        coordinate.extend(self.__carbon.get_coordinate())
        coordinate.extend(self.__oxygen.get_coordinate())
        coordinate.extend(self.__nitrogen.get_coordinate())
        self.__coordinate = tuple(coordinate)

    def get_coordinate(self):
        return self.__coordinate

    def get_target(self):
        return self.__target


def read_residues(pdbfile):
    residues = dict()
    with open(pdbfile, 'r') as pdb:
        pdb_file = tuple(atom.strip() for atom in pdb.readlines()[1: -2])
    for atom in pdb_file:
        identifier = atom[0: 6].strip()
        if identifier == 'ATOM':
            index = atom[22: 26].strip()
            atom_type = atom[13]
            if index not in residues:
                residues[index] = []
            residues[index].append((atom_type, tuple(map(float, (atom[30: 38].strip(),
                                                                 atom[38: 46].strip(),
                                                                 atom[46: 54].strip())))))
    return residues


def read_aalist(file):
    aalist = dict()
    with open(file, 'r') as aa_list:
        for line in aa_list:
            words = line.split()
            if words[0] not in aalist:
                aalist[words[0]] = []
            aalist[words[0]].append(words[1])
    return aalist


def cal_similarity(targetlist):
    result = dict()
    for i in range(len(targetlist)):
        suj_similarity = dict()
        ref = targetlist[i].get_coordinate()
        for j in range(len(targetlist)):
            can = targetlist[j].get_coordinate()
            dist = math.sqrt(sum(((ref[k] - can[k]) ** 2 for k in range(36)))) / 36
            similarity = float(1) / (1 + dist)
            suj_similarity[targetlist[j].get_target()] = similarity
        result[targetlist[i].get_target()] = suj_similarity
    return result


def read_csv(file, heading=True):
    if heading:
        with open(file, 'r') as csv:
            lines = tuple(tuple(line.strip().split(',')) for line in csv.readlines()[1:])
        return lines
    else:
        with open(file, 'r') as csv:
            lines = tuple(tuple(line.strip().split(',')) for line in csv.readlines())
        return lines


AA_List = read_aalist('Z:\\project\\residue.txt')
BP_list = []
SM = []
path = 'Z:\\project\\Targets\\'
for target in AA_List:
    print(target)
    Carbons = []
    Oxygens = []
    Nitrogens = []
    Atoms = []
    files = os.listdir(path + target)
    if 'protein_nowater_nocofactor.pdb' in files:
        Residues = read_residues(path + target + '\\protein_nowater_nocofactor.pdb')
    elif 'protein_nowater_2.pdb' in files:
        Residues = read_residues(path + target + '\\protein_nowater_2.pdb')
    else:
        Residues = read_residues(path + target + '\\protein_nowater.pdb')
    for KeyResidue in AA_List[target]:
        try:
            resn = KeyResidue.split('-')[1]
            for atom in Residues[resn]:
                # Atoms.append(atom[1])
                if atom[0] == 'C':
                    Carbons.append(atom[1])
                elif atom[0] == 'O':
                    Oxygens.append(atom[1])
                elif atom[0] == 'N':
                    Nitrogens.append(atom[1])
        except KeyError:
            print(KeyResidue)
    BP = AtomTypeDescriptor(carbons=Carbons, oxygens=Oxygens, nitrogens=Nitrogens, target=target)
    # print(BP.get_coordinate())
    SM.append(BP.get_coordinate())
    BP_list.append(BP)

SM = np.array(SM)
fig = ff.create_dendrogram(X=SM, labels=list(AA_List.keys()), color_threshold=45)
fig['layout'].update({'width': 1240})
fig['layout'].update({'height': 500})

"""csv_file = 'Z:\\project\\Metadata\\BindingSitesAnalysis.csv'
csv_lines = read_csv(csv_file)

ranking = tuple(fig['layout']['xaxis']['ticktext'])

AA_dict = {}
for line in csv_lines:
    AA_dict[line[0]] = line[1:]
rows = []
for aa in ranking:
    rows.append(AA_dict[aa])
rows.reverse()
with open('Z:\\project\\Metadata\\AA_clustering_atomtype.csv', 'w+', newline='') as w:
    w = csv.writer(w)
    w.writerows(rows)"""

plot(fig)

"""Result = cal_similarity(BP_list)
x_axis = [' ']
x_axis.extend(tuple(AA_List.keys()))
similarity_matrix = tuple(tuple(Result[receptor].values()) for receptor in Result)
with open('Z:\\dissertation\\fixed_AtomType.csv', 'w+', newline='') as csv_file:
    writer = csv.writer(csv_file)
    writer.writerow(x_axis)
    for i in range(len(similarity_matrix)):
        row = [x_axis[i + 1]]
        row.extend(similarity_matrix[i])
        writer.writerow(row)"""
