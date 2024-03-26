#!/usr/bin/env python
# coding: utf-8

def make_xyz_str(symbols, coords):
    xyz_str = ''
    for s, c in zip(symbols, coords):
        xyz_str = xyz_str + f'{s}  {c[0]: .10f}  {c[1]: .10f}  {c[2]: .10f}\n'
    return xyz_str