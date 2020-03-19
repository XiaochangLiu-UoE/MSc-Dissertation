import os
import csv
from BindingSiteClassification.Descriptor import *


def cal_similarity(targetlist):
    result = dict()
    for i in range(len(targetlist)):
        suj_similarity = dict()
        ref = targetlist[i].get_coordinate()
        for j in range(len(targetlist)):
            can = targetlist[j].get_coordinate()
            dist = math.sqrt(sum(((ref[k] - can[k]) ** 2 for k in range(12)))) / 12
            similarity = float(1) / (1 + dist)
            suj_similarity[targetlist[j].getTarget()] = similarity
        result[targetlist[i].getTarget()] = suj_similarity
    return result


def read_residues(pdbfile):
    residues = dict()
    with open(pdbfile, 'r') as pdb:
        pdb_file = tuple(atom.strip() for atom in pdb.readlines()[1: -2])
    for atom in pdb_file:
        identifier = atom[0: 6].strip()
        if identifier == 'ATOM':
            index = atom[22: 26].strip()
            if index not in residues:
                residues[index] = []
            residues[index].append(tuple(map(float, (atom[30: 38].strip(),
                                                     atom[38: 46].strip(),
                                                     atom[46: 54].strip()))))
    return residues


class SimpleDescriptor(Descriptor):
    def __init__(self, atoms, target):
        Descriptor.__init__(self, atoms)
        self.target = target

    def get_target(self):
        return self.target


AA_List = dict()
path = 'Z:\\project\\Targets\\'
with open('Z:\\project\\residue.txt', 'r') as aa_list:
    for line in aa_list:
        words = line.split()
        if words[0] not in AA_List:
            AA_List[words[0]] = []
        AA_List[words[0]].append(words[1])

BP_list = []
for target in AA_List:
    print(target)
    files = os.listdir(path + target)
    KeyResidues = []
    if 'protein_nowater_nocofactor.pdb' in files:
        Residues = read_residues(path + target + '\\protein_nowater_nocofactor.pdb')
    elif 'protein_nowater_2.pdb' in files:
        Residues = read_residues(path + target + '\\protein_nowater_2.pdb')
    else:
        Residues = read_residues(path + target + '\\protein_nowater.pdb')
    for KeyResidue in AA_List[target]:
        try:
            resn = KeyResidue.split('-')[1]
            KeyResidues.extend(Residues[resn])
        except KeyError:
            print(KeyResidue)
    BP = SimpleDescriptor(atoms=KeyResidues, target=target)
    print(BP.get_coordinate())
    BP_list.append(BP)

Result = cal_similarity(BP_list)

x_axis = [' ']
x_axis.extend(tuple(AA_List.keys()))
similarity_matrix = tuple(tuple(Result[receptor].values()) for receptor in Result)

with open('Z:\\dissertation\\test.csv', 'w+', newline='') as csv_file:
    writer = csv.writer(csv_file)
    writer.writerow(x_axis)
    for i in range(len(similarity_matrix)):
        row = [x_axis[i + 1]]
        row.extend(similarity_matrix[i])
        writer.writerow(row)
