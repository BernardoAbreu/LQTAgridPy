#!/usr/bin/env python3
# coding: utf-8

import matrix_generate
import numpy as np
import pandas as pd
import os
import re


class GridGenerate():

    def __init__(self, coordinates, dimensions, atp, directory, step):
        n_pat = re.compile(r'(.*[0-9]+).+')
        # dataFile = open(files).read().splitlines()
        dataFile = os.listdir(directory)
        self.molecules = [x.replace('.gro', '')
                          for x in dataFile if x.endswith('gro')]

        def __sort_replace(x):
            return n_pat.search(x).group(1)

        self.molecules.sort(__sort_replace)
        dataFile = [directory + '/' + fileName for fileName in dataFile]

        groFiles = [x for x in dataFile if x.endswith('gro')]
        groFiles.sort(key=__sort_replace)
        topFiles = [x for x in dataFile if x.endswith('top')]
        topFiles.sort(key=__sort_replace)
        itpFiles = [x for x in dataFile if x.endswith('nb.itp')]
        itpFiles.sort(key=__sort_replace)

        matrices = []

        minimos = np.array([999999.0, 999999.0, 999999.0])
        maximos = np.array([-999999.0, -999999.0, -999999.0])
        # for fileGro, fileItp, fileTop in utils.pairwise(dataFile):
        # for fileGro, fileTop, fileItp in utils.triplewise(dataFile):
        for i in range(len(groFiles)):
            matrix = matrix_generate.MatrixGenerate(groFiles[i],
                                                    topFiles[i],
                                                    itpFiles[i])
            minimos = np.minimum(minimos, matrix.minimos)
            maximos = np.maximum(maximos, matrix.maximos)
            matrices.append(matrix)

        x0, y0, z0 = coordinates if coordinates else (maximos.astype(int) - 5)

        dim_x, dim_y, dim_z = dimensions if dimensions else \
            ((maximos - minimos).astype(int) + 10)

        self.cCoulomb = []
        self.cLJ = []
        # esse loop gera o cabeçalho do arquivo de saida
        for atp_l in atp:
            for i in np.arange(x0, dim_x + x0 + step, step):
                for j in np.arange(y0, dim_y + y0 + step, step):
                    for k in np.arange(z0, dim_z + z0 + step, step):
                        self.cCoulomb.append(
                            "%.2f_%.2f_%.2f_%s_C:" % (i, j, k, atp_l))
                        self.cLJ.append("%.2f_%.2f_%.2f_%s_LJ:" %
                                        (i, j, k, atp_l))
        self.output = ' \t'.join(self.cCoulomb + self.cLJ)

        self.coulombMatrix = []
        self.ljMatrix = []
        for matrix in matrices:
            matrix.gridGenerate(dim_x, dim_y, dim_z, atp, x0, y0, z0, step)
            # valuesCoulomb = matrix.getMatrix("C")
            # valuesLj = matrix.getMatrix("L")
            txt_val_c, txt_val_lj, coulombMatrix, ljMatrix = matrix.getMatrix()
            self.output += "\n" + txt_val_c + txt_val_lj
            self.coulombMatrix.append(coulombMatrix)
            self.ljMatrix.append(ljMatrix)

    def saveGrid(self, output):
        arq = open(output + '.txt', "w")
        arq.write(self.output)
        arq.close()
        dfCoulomb = pd.DataFrame(self.coulombMatrix,
                                 columns=self.cCoulomb,
                                 index=self.molecules)
        dfLj = pd.DataFrame(self.ljMatrix,
                            columns=self.cLJ,
                            index=self.molecules)
        df = dfCoulomb.join(dfLj)
        df.to_csv(output + '.csv', sep=';')
