#!/usr/bin/env
# coding: utf-8

import math
import os
import numpy as np
from scipy.spatial import ConvexHull
import generate_points


def isclose(a, b, rel_tol=1e-9, abs_tol=0.0):
    return abs(a - b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


class MatrixGenerate():

    def __init__(self, fileGro, fileTop, fileItp):
        # print('hello')
        self.setX(fileGro)
        self.atomsTypes(fileTop)
        self.loadConstants(fileItp)
        self.loadAP()
        self.determineConstants()

    def setX(self, fileName):
        print(fileName)
        start_coords = 3

        self.c6 = []
        self.c12 = []

        coords = []

        with open(fileName, 'r') as f:
            while f.readline():
                self.numberElements = int(f.readline())
                for i in range(self.numberElements):
                    coords.append(np.array(f.readline().split()[start_coords:],
                                           dtype=np.float16))
                f.readline()

        # cria-se matriz com 3 linhas (uma para cada coordenada x, y e z)
        # e as colunas como coordenadas
        self.X = np.array(coords) * 10

        # numero de linhas da matriz gerada
        self.m = self.X.shape[0]

        self.minimos = np.amin(self.X, axis=0)
        self.maximos = np.amax(self.X, axis=0)

    def atomsTypes(self, fileName):
        # print(fileName)
        with open(fileName) as f:
            line = f.readline()
            while '[ atoms ]' not in line:
                line = f.readline()
            next(f)

            self.types = []
            self.cargas = np.empty(self.numberElements)
            # print(self.numberElements)
            for i in range(self.numberElements):
                line = f.readline()
                atoms = line.split()
                # print(atoms)
                self.types.append(atoms[1])
                self.cargas[i] = atoms[6]

    def loadConstants(self, fileName):
        self.typeConstants = []
        sigma = []
        epsilon = []

        with open(fileName) as f:
            line = f.readline()
            while '[ atomtypes ]' not in line:
                line = f.readline()
            next(f)
            for line in f:
                atoms = line.split()
                self.typeConstants.append(atoms[0])
                sigma.append(atoms[4])
                epsilon.append(atoms[5])

        sigma = np.array(sigma, dtype=float)
        epsilon = np.array(epsilon, dtype=float)

        self.constantc6 = (epsilon * (sigma ** 6)) * 4.
        self.constantc12 = (epsilon * (sigma ** 12)) * 4.

    def loadAP(self):
        filename = os.path.dirname(__file__) + '/defaultsFiles/AtomProva.atp'
        # filename = 'defaultsFiles/AtomProva.atp'
        self.ap = {}
        with open(filename, 'r') as f:
            next(f)
            next(f)
            for line in f:
                group, charge, c6, c12 = line.split()
                self.ap[group] = {'carga': float(charge), 'c6': float(c6),
                                  'c12': float(c12)}

    def determineConstants(self):
        for c_type in self.types:
            index = self.typeConstants.index(c_type)
            self.c6.append(self.constantc6[index])
            self.c12.append(self.constantc12[index])
        self.c6 = np.array(self.c6)
        self.c12 = np.array(self.c12)

    def gridGenerate(self, dimX, dimY, dimZ, atp, x0, y0, z0, step):
        self.DimX = dimX
        self.DimY = dimY
        self.DimZ = dimZ
        self.natp = len(atp)

        f = 138.935485  # conversion factor

        nframes = self.m / self.numberElements

        self.gridCoulomb = {}
        self.gridLJ = {}

        # esse loop roda sobre o número de sondas escolhidas
        for h in range(self.natp):
            # carrega-se as respectivas constantes
            q1 = self.ap[atp[h]]['carga']  # self.cargasap[elem]
            c6a = self.ap[atp[h]]['c6']  # self.c6ap[elem]
            c12a = self.ap[atp[h]]['c12']  # self.c12ap[elem]
            Vlj = 0
            Vc = 0
            npontos = 0

            r1 = np.zeros(3)
            self.gridCoulomb[atp[h]] = {}
            self.gridLJ[atp[h]] = {}

            # aqui começa o loop que gera as coordenadas cartesianas
            # acho que você pode gerar os pontos com o fecho convexo e
            # substituir esses 3 loops por um loop sobre os pontos gerados
            for i in np.arange(x0, self.DimX + x0 + step, step):
                r1[0] = i + x0
                self.gridCoulomb[atp[h]][i] = {}
                self.gridLJ[atp[h]][i] = {}
                for j in np.arange(y0, self.DimY + y0 + step, step):
                    r1[1] = j + y0
                    self.gridCoulomb[atp[h]][i][j] = {}
                    self.gridLJ[atp[h]][i][j] = {}
                    for k in np.arange(z0, self.DimZ + z0 + step, step):
                        r1[2] = k + z0
                        Vlj = 0
                        Vc = 0
                        npontos += 1
                        self.gridCoulomb[atp[h]][i][j][k] = {}
                        self.gridLJ[atp[h]][i][j][k] = {}

                        # geradas as coordenadas cartesianas começa o loop
                        # sobre os átomos do PAC para calcular os descriotres
                        # com base na distância entre a sonda e os átomos
                        for l in range(self.m):
                            r = np.linalg.norm(r1 - self.X[l]) / 10.
                            index = l % self.numberElements
                            c6ij = math.sqrt(c6a * self.c6[index])
                            c12ij = math.sqrt(c12a * self.c12[index])

                            if r != 0:
                                Vlj += (c12ij / (r ** 12)) - (c6ij / (r ** 6))
                                Vc += f * float(q1) * \
                                    float(self.cargas[index]) / r
                            else:
                                Vlj = float("inf")
                                Vc = float("inf")

                        self.gridCoulomb[atp[h]][i][j][k] = Vc / nframes
                        self.gridLJ[atp[h]][i][j][k] = Vlj / math.sqrt(nframes)

    def hullGenerate(self, atp, step, initial_distance, total_layers, delta_r):
        self.natp = len(atp)

        f = 138.935485  # conversion factor

        nframes = self.m / self.numberElements

        # esse loop roda sobre o número de sondas escolhidas
        for h in range(self.natp):
            # carrega-se as respectivas constantes
            q1 = self.ap[atp[h]]['carga']
            c6a = self.ap[atp[h]]['c6']
            c12a = self.ap[atp[h]]['c12']

            hull = ConvexHull(self.X)
            self.points = generate_points.generate_points(hull, step,
                                                          initial_distance,
                                                          total_layers,
                                                          delta_r)

            self.hullCoulombList = np.empty(self.points.shape[0])
            self.hullLJList = np.empty(self.points.shape[0])

            c6ij = np.sqrt(c6a * self.c6)
            c12ij = np.sqrt(c12a * self.c12)
            for p_index, point in enumerate(self.points):
                r = np.linalg.norm(point - self.X, axis=1) / 10.
                r6 = r ** 6
                Vlj = ((c12ij /
                        (r6 ** 2).reshape(int(nframes), self.numberElements)) -
                       (c6ij / r6.reshape(int(nframes), self.numberElements))
                       ).ravel().sum()

                Vc = f * q1 * (self.cargas /
                               r.reshape(int(nframes), self.numberElements)
                               ).ravel().sum()

                self.hullCoulombList[p_index] = Vc / nframes
                self.hullLJList[p_index] = Vlj / math.sqrt(nframes)

    def getMatrix(self):
        textValuesCoulomb = ""
        textValuesLj = ""
        coulombMatrix = []
        ljMatrix = []
        count = 0
        for h in self.gridCoulomb:
            for i in self.gridCoulomb[h]:
                for j in self.gridCoulomb[h][i]:
                    for k in self.gridCoulomb[h][i][j]:
                        textValuesCoulomb += "%g\t" % \
                            (self.gridCoulomb[h][i][j][k])
                        textValuesLj += "%g\t" % (self.gridLJ[h][i][j][k])
                        coulombMatrix.append(self.gridCoulomb[h][i][j][k])
                        ljMatrix.append(self.gridLJ[h][i][j][k])
                        count += 1
        return textValuesCoulomb, textValuesLj, coulombMatrix, ljMatrix

# TODO
    def getHullList(self):
        return self.points, self.hullCoulombList, self.hullLJList
