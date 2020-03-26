import math


def cal_coordinate(distlist):
    mean_dist = sum(tuple(distlist.values()))/len(distlist)
    var_dist = sum(((dist - mean_dist) ** 2 for dist in tuple(distlist.values()))) / (len(distlist) - 1)
    sd = math.sqrt(var_dist)
    skew_dist = sum(((dist - mean_dist) / sd) ** 3 for dist in tuple(distlist.values())) / len(distlist)
    return mean_dist, var_dist, skew_dist


def pack_coordinates(distlist1, distlist2, distlist3, distlist4):
    coordinate = []
    for dl in (distlist1, distlist2, distlist3, distlist4):
        coordinate.extend(cal_coordinate(distlist=dl))
    return tuple(coordinate)


def cal_dist(atom, ref):
    dist = math.sqrt((atom[0] - ref[0]) ** 2 + (atom[1] - ref[1]) ** 2 + (atom[2] - ref[2]) ** 2)
    return dist


def cal_cog(atoms):
    x = 0
    y = 0
    z = 0
    for atom in atoms:
        x = x + atom[0]
        y = y + atom[1]
        z = z + atom[2]
    cog = (x / len(atoms), y / len(atoms), z / len(atoms))
    return cog


def accu_dist(atoms, ref):
    dists = dict()
    for atom in atoms:
        dist = cal_dist(atom=atom, ref=ref)
        dists[atom] = dist
    return dists


class Descriptor(object):
    def __init__(self, atoms):
        if len(atoms) != 0:
            self.__cog = cal_cog(atoms=atoms)
            self.__cog_dists = accu_dist(atoms=atoms, ref=self.__cog)
            self.__cog_fur = max(self.__cog_dists, key=lambda k: k[1])
            self.__cog_near = min(self.__cog_dists, key=lambda k: k[1])
            self.__cog_near_dists = accu_dist(atoms=atoms, ref=self.__cog_near)
            self.__cog_fur_dists = accu_dist(atoms=atoms, ref=self.__cog_fur)
            self.__cog_fur_fur = max(self.__cog_fur_dists, key=lambda k: k[1])
            self.__cog_fur_fur_dists = accu_dist(atoms=atoms, ref=self.__cog_fur_fur)
            self.__coordinate = pack_coordinates(self.__cog_dists, self.__cog_near_dists,
                                                 self.__cog_fur_dists, self.__cog_fur_fur_dists)
        else:
            self.__cog = None
            self.__cog_near = None
            self.__cog_fur = None
            self.__cog_fur_fur = None
            self.__cog_dists = None
            self.__cog_near_dists = None
            self.__cog_fur_dists = None
            self.__cog_fur_fur_dists = None
            self.__coordinate = tuple(99 for i in range(12))

    def get_coordinate(self):
        return self.__coordinate

