#!/usr/bin/env
# coding: utf-8

from itertools import izip
import math

def Distance(r1, r2):
    d = math.sqrt(math.pow((r1[0] - r2[0]), 2)
            + math.pow((r1[1] - r2[1]) ,2) + math.pow((r1[2] - r2[2]), 2))
    return d

def pairwise(iterable):
    a = iter(iterable)
    return izip(a, a)
    