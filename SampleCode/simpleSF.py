#!/usr/bin/python
import csv
import multiprocessing as mp
import time
import os
from Descriptor.calEntropy import *


def split_jobs(jobs, part=3):
    ratio = float(1) / part
    job = []
    for each_part in range(0, part):
        job.append(
            jobs[int(math.ceil(ratio * each_part * len(jobs))): int(math.ceil(ratio * (each_part + 1) * len(jobs)))])
    return tuple(job)


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
