import csv
import multiprocessing as mp
import time
from CalEntropy.calEntropy import *
from Rescoring.different_weights_PIP_ENERGY import split_jobs


"""def cal_entropy(sdf, protein):
    ligand = get_sdf_atoms(sdf, single=False)
    water_entropy, side_chain_entropy = cal_lw_fs(ligand=ligand, protein=protein)
    return water_entropy, side_chain_entropy"""


def collect(entropy):
    global Results
    Results.append(entropy)


def multi_cal_entropy(job, patms):
    entropies = []
    for mol in job:
        try:
            latms = get_sdf_atoms(mol[1], single=False)
            entropies.append((mol[0], cal_lw_fs(latms, patms)))
        except IndexError:
            print(mol[0])
    return tuple(entropies)


def multi_collect(entropies):
    global Multi_results
    Multi_results.extend(entropies)


if __name__ == "__main__":
    Results = []
    Multi_results = []
    multi_sdf = []
    protein_atoms = get_pdb_atoms("ERBB2\\protein_nowater.pdb")

    with open("ERBB2\\output.sdf", "r") as SDF:
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
    # set the number of molecules
    mols = multi_sdf
    # split jobs
    start = time.time()
    p = mp.Pool(mp.cpu_count())
    jobs = split_jobs(mols, part=mp.cpu_count())
    for queue in jobs:
        p.apply_async(func=multi_cal_entropy, args=(queue, protein_atoms), callback=multi_collect)
    p.close()
    p.join()
    print("%s" % (time.time() - start))
    field_names = ("Name", "Enthalpy", "vDW", "HBD", "HBA", "Lost Water entropy", "Side chain conformational entropy")
    with open("ERBB2\\Entropy.csv", "w+", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(field_names)
        for entropy in Multi_results:
            writer.writerow((entropy[0][0], entropy[0][1], entropy[0][2], entropy[0][3], entropy[0][4], entropy[1][0],
                             entropy[1][1] * 0.01))
    # one by one
    """start = time.time()
    p = mp.Pool(mp.cpu_count())
    for sdf in mols:
        p.apply_async(func=cal_entropy, args=(sdf, protein_atoms), callback=collect)
    p.close()
    p.join()
    print("%s" % (time.time() - start))

    Results = collections.Counter(Results)
    Multi_results = collections.Counter(Multi_results)
    if Results == Multi_results:
        print("Identical")"""

    # single process
    """start = time.time()
    for sdf in mols:
        ligand_atoms = get_sdf_atoms(sdf, single=False)
        water, side = cal_lw_fs(ligand_atoms, protein_atoms)
    print("%s" % (time.time() - start))"""
