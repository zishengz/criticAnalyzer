#!/usr/bin/env python3 -u
# --------------------------------------------------------
# CP info extractor and printer for CP search stdout
#   using critic2, periodic systems
# 
# Copyright (C) 2021 Zisheng Zhang
# Email: zisheng@chem.ucla.edu
# --------------------------------------------------------
# Dependencies: numpy, ase
# Usage:
# $ python cp_printer.py cps.crit.out
# --------------------------------------------------------

from ase.io import read, write
import numpy as np
import sys


def get_cp_flags(out_raw):
    '''Get the line # flag for key printout blocks'''
    flag_cp_list_unique = None
    flag_bond_list_unique = None
    flag_cp_groups = None
    flag_cp_list_all = None
    flag_connect = None
    for i in range(len(out_raw)):
        if flag_cp_list_unique is None:
            if "* Critical point list, final report (non-equivalent cps)" in out_raw[i]:
                flag_cp_list_unique = i
        if flag_bond_list_unique is None:
            if "* Analysis of system bonds" in out_raw[i]:
                flag_bond_list_unique = i
        if flag_cp_groups is None:
            if "# (x symbols are the non-equivalent representative atoms)" in out_raw[i]:
                flag_cp_groups = i
        if flag_cp_list_all is None:
            if "* Complete CP list, bcp and rcp connectivity table" in out_raw[i]:
                flag_cp_list_all = i
        if flag_connect is None:
            if "* Attractor connectivity matrix" in out_raw[i]:
                flag_connect = i
    return flag_cp_list_unique, flag_bond_list_unique, flag_cp_groups, flag_cp_list_all, flag_connect


def get_cp_list_unique(out_raw, flag):
    '''CP list, non-equivalent'''
    spacegroup_unique = []
    cp_type_unique = []
    rho_unique = []
    grad_unique = []
    lap_unique = []
    for i in range(flag+4, len(out_raw)):
        l = out_raw[i]
        if l[0] != '\n':
            l = l.split()
            spacegroup_unique.append(l[1])
            cp_type_unique.append(l[-9][0])
            rho_unique.append(eval(l[-3]))
            grad_unique.append(eval(l[-2]))
            lap_unique.append(eval(l[-1]))
        else:
            break
    return spacegroup_unique, cp_type_unique, rho_unique, grad_unique, lap_unique
    #spcgrp, cptype, rho, grad, lap


def get_bond_list_unique(out_raw, flag):
    '''Bond list, non-equivalent'''
    bond_r1 = []
    bond_r2 = []
    angle_r1Br2 = []
    for i in range(flag+7, len(out_raw)):
        l = out_raw[i]
        if l[0] != '\n':
            l = l.split()
            bond_r1.append(eval(l[-4]))
            bond_r2.append(eval(l[-3]))
            angle_r1Br2.append(eval(l[-1]))
        else:
            break
    return bond_r1, bond_r2, angle_r1Br2


def get_cp_list_all(out_raw, flag):
    '''CP list, including equivalent ones'''
    ncp = []
    bcp = []
    for i in range(flag+3, len(out_raw)):
        l = out_raw[i]
        if l[0] != '\n':
            if l.split()[2] == 'n':
                ncp.append(l)
            elif l.split()[2] == 'b':
                bcp.append(l)
            else:
                break
        else:
            break
    return ncp, bcp


def get_nat_index(ncp, bcp, coord_direct):
    '''Get the indeces according to POSCAR'''
    ncp_index = []
    for i in range(len(ncp)):
        ncp_coord = np.array([eval(j) for j in ncp[i].split()[-3:]])
        exist_match = False
        for j in range(len(coord_direct)):
            if np.abs(coord_direct[j]-ncp_coord).sum() < 0.01 or np.abs(coord_direct[j]-ncp_coord+np.array([0, 0, 1])).sum() < 0.01 or np.abs(coord_direct[j]-ncp_coord+np.array([0, 1, 0])).sum() < 0.01 or np.abs(coord_direct[j]-ncp_coord+np.array([1, 0, 0])).sum() < 0.01:
                ncp_index.append(j)
                exist_match = True
                break
        if not exist_match:
            ncp_index.append(f'X {i}')
    bcp_index = []
    bond_pair_index = []
    for i in range(len(bcp)):
        l = bcp[i]
        bcp_index.append(eval(l.split()[0]))
        crit_ind1 = eval(l.split('(')[0].split()[-1])
        crit_ind2 = eval(l.split('(')[1].split()[-1])
        bond_pair_index.append([
            crit_ind1, crit_ind2
            #ncp_index[crit_ind1-1], ncp_index[crit_ind2-1]
        ])
    return ncp_index, bcp_index, bond_pair_index


def get_cp_groups(out_raw, flag):
    '''Get th groups of equivalent CPs from symmetry'''
    groups = []
    for i in range(flag+2, len(out_raw)):
        l = out_raw[i]
        if l[0] == '\n':
            break
        else:
            if l[0] == 'x':
                tmp = [eval(l.split()[1])]
                groups.append(tmp)
            else:
                tmp.append(eval(l.split()[0]))
    return groups


def print_bcp_list(elem_list, ncp_index, bcp_index, bond_pair_index, groups, spcgrp, cptype, rho, grad, lap, bond_r1, bond_r2, angle_r1Br2):
    '''Print formated CP info'''
    print(f'   B O N D       Rho     |Grad|      Lap      r1     r2    r1/r1  r1-B-r2')
    ncp_num = cptype.count('n')
    for b in range(len(bcp_index)):
        bp = bond_pair_index[b]
        bp = [ncp_index[j-1] for j in bp]
        mygroup = None
        for i in range(len(groups)):
            if bcp_index[b] in groups[i]:
                mygroup = i
                break
        if mygroup is not None:
            if type(bp[0]) is str:
                elem1 = bp[0]
            else:
                elem1 = elem_list[bp[0]]+' '+str(bp[0])
            if type(bp[1]) is str:
                elem2 = bp[1]
            else:
                elem2 = elem_list[bp[1]]+' '+str(bp[1])
            # remove bad bcp
            if lap[mygroup] > 1000 or lap[mygroup] < -1000:
                continue
            print(
                f'{elem1:<5} - {elem2:<5}  {rho[mygroup]:.6f}  {grad[mygroup]:.4f}  {lap[mygroup]:10.6f}  {bond_r1[mygroup-ncp_num]:.3f}  {bond_r2[mygroup-ncp_num]:.3f}  {bond_r1[mygroup-ncp_num]/bond_r2[mygroup-ncp_num]:.3f}  {angle_r1Br2[mygroup-ncp_num]:.3f}')


fname_inp = sys.argv[1]
crystal = read('POSCAR')
elem_list = crystal.get_chemical_symbols()
coord_direct = crystal.get_scaled_positions()
out_raw = open(fname_inp).readlines()


# print bcp info
flag_cp_list_unique, flag_bond_list_unique, flag_cp_groups, flag_cp_list_all, flag_connect = get_cp_flags(
    out_raw)
spcgrp, cptype, rho, grad, lap = get_cp_list_unique(
    out_raw, flag_cp_list_unique)
bond_r1, bond_r2, angle_r1Br2 = get_bond_list_unique(out_raw, flag_bond_list_unique)

groups = get_cp_groups(out_raw, flag_cp_groups)
ncp, bcp = get_cp_list_all(out_raw, flag_cp_list_all)
ncp_index, bcp_index, bond_pair_index = get_nat_index(ncp, bcp, coord_direct)
print_bcp_list(elem_list, ncp_index, bcp_index, bond_pair_index,
               groups, spcgrp, cptype, rho, grad, lap, bond_r1, bond_r2, angle_r1Br2)
