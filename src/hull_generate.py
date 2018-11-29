#!/usr/bin/env python3
# coding: utf-8

import matrix_generate
import numpy as np
import pandas as pd
import os
import re


class HullGenerate():

    def __init__(self, atp, directory, delta_angle, initial_distance, delta_r,
                 total_layers):
        self.n_pat = re.compile(r'(.*[0-9]+).+')
        # dataFile = open(files).read().splitlines()
        dataFile = os.listdir(directory)
        self.molecules = [x.replace('.gro', '')
                          for x in dataFile if x.endswith('gro')]
        self.molecules = np.array(sorted(self.molecules, key=self.__sort_replace))
        dataFile = [directory + '/' + fileName for fileName in dataFile]

        groFiles = [x for x in dataFile if x.endswith('gro')]
        groFiles.sort(key=self.__sort_replace)
        topFiles = [x for x in dataFile if x.endswith('top')]
        topFiles.sort(key=self.__sort_replace)
        itpFiles = [x for x in dataFile if x.endswith('nb.itp')]
        itpFiles.sort(key=self.__sort_replace)

        matrices = [matrix_generate.MatrixGenerate(gro, top, itp)
                    for gro, top, itp in zip(groFiles, topFiles, itpFiles)]

        self.cCoulomb = []
        self.cLJ = []

        delta_r /= 10
        initial_distance /= 10
        self.coulombMatrix = []
        self.ljMatrix = []
        self.points = []
        for mols, matrix in zip(self.molecules[:5], matrices[:5]):
            print('Calculating ' + mols)
            matrix.hullGenerate(atp, delta_angle, initial_distance,
                                total_layers, delta_r)

            points, coulombAtom, ljAtom = matrix.getHullList()
            # self.output += "\n" + txt_val_c + txt_val_lj
            # print(points[0])
            self.points.append([','.join(map(str, p)) for p in points])
            self.coulombMatrix.append(coulombAtom)
            self.ljMatrix.append(ljAtom)
            # break
        # self.coulombMatrix = np.array(self.coulombMatrix)
        # print(self.coulombMatrix.shape)
        # self.ljMatrix = np.array(self.ljMatrix)
        # print(self.ljMatrix.shape)

    def saveGrid(self, output):
        np.savetxt(output + '.mol', self.molecules, fmt='%s')
        dfCoulomb = pd.DataFrame(self.coulombMatrix, index=self.molecules[:5])
        dfLj = pd.DataFrame(self.ljMatrix, index=self.molecules[:5])
        dfPoints = pd.DataFrame(self.points, index=self.molecules[:5])
        dfCoulomb.to_csv(output + '_C.csv', sep=';', header=False)
        dfLj.to_csv(output + '_LJ.csv', sep=';', header=False)
        dfPoints.to_csv(output + '_Points.csv', sep=';', header=False)

    def __sort_replace(self, x):
        return self.n_pat.search(x).group(1)
