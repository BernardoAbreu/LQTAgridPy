#!/usr/bin/env python3
# coding: utf-8

import click
import grid_generate
import hull_generate
# from . import grid_generate


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--mols',
              metavar='<path>',
              type=click.Path(exists=True),
              required=True,
              help='files path, gro and itp.'
              )
# @click.option('--coordinates', '-c',
#               metavar='<x> <y> <z>',
#               type=int,
#               nargs=3,
#               required=False,
#               help='Coordinates of the box.'
#               )
# @click.option('--dimensions', '-d',
#               metavar='<x> <y> <z>',
#               type=int,
#               nargs=3,
#               required=False,
#               help='Dimensions of the box.'
#               )
@click.option('--atom', '-a',
              metavar='[atom]',
              multiple=True,
              required=True,
              help='Probe.'
              )
# @click.option('--step', '-s',
#               metavar='<x>',
#               type=float,
#               nargs=1,
#               required=True,
#               help='Steps for navegation on matrix.'
#               )
@click.option('--initial', '-i',
              metavar='<x>',
              type=float,
              nargs=1,
              required=True,
              help='Initial distance from molecule hull.'
              )
@click.option('--angle', '-d',
              metavar='<x>',
              type=int,
              nargs=1,
              required=True,
              help='Angle step.'
              )
@click.option('--radius', '-r',
              metavar='<x>',
              type=float,
              nargs=1,
              required=True,
              help='Radius step.'
              )
@click.option('--output', '-o',
              metavar='<path_output>',
              type=click.Path(),
              required=True,
              help='Output matrix file.'
              )
def run(mols, atom, angle, initial, radius, output):
    '''LQTAgridPy is a python version of LQTAgrid, a practical application of
    4D QSAR analysis methodology developed at University of Campinas.

    More: https://github.com/rougeth/LQTAgridPy
    '''
    # print('hello')
    # grid = grid_generate.GridGenerate(
    #     coordinates,
    #     dimensions,
    #     atom,
    #     mols,
    #     step
    # )
    hull = hull_generate.HullGenerate(atom, mols, angle, initial, radius, 7)
    hull.saveGrid(output)
    # grid.saveGrid(output)


if __name__ == '__main__':
    run()
