#!/usr/bin/env python3
# coding: utf-8

import argparse
import hull_generate


def run(mols, atom, angle, initial, radius, output):
    '''LQTAgridPy is a python version of LQTAgrid, a practical application of
    4D QSAR analysis methodology developed at University of Campinas.

    More: https://github.com/rougeth/LQTAgridPy
    '''
    hull = hull_generate.HullGenerate(atom, mols, angle, initial, radius, 7)
    hull.saveGrid(output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--mols', required=True,
                        metavar='<path>',
                        help='files path, gro and itp.'
                        )
    parser.add_argument('--atom', '-a', required=True,
                        metavar='[atom]',
                        action='append',
                        help='Probe.'
                        )
    parser.add_argument('--initial', '-i', required=True, type=float,
                        metavar='<x>',
                        help='Initial distance from molecule hull.'
                        )
    parser.add_argument('--angle', '-d', required=True, type=int,
                        help='Angle step.'
                        )
    parser.add_argument('--radius', '-r', required=True, type=float,
                        metavar='<x>',
                        help='Radius step.'
                        )
    parser.add_argument('--output', '-o', required=True,
                        metavar='<path_output>',
                        help='Output matrix file.'
                        )
    args = parser.parse_args()
    mols = args.mols
    atom = args.atom
    angle = args.angle
    initial = args.initial
    radius = args.radius
    output = args.output
    run(mols, atom, angle, initial, radius, output)
