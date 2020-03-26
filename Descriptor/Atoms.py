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
