# -*- coding:utf-8 -*-
import sys
import math
import numpy as np
import random

#lattice constant 10 and strain settings
LAT_CONST = 2.855324
POISSON = 0.293
STRAIN_X = 0.000
STRAIN_Y = -POISSON * STRAIN_X
STRAIN_Z = -POISSON * STRAIN_X

DIST_X = LAT_CONST
DIST_Y = LAT_CONST
DIST_Z = LAT_CONST

# the size of the whole simulation cell
# N_X, N_Y and N_Z represents the length of the cell (unit: lattice constant)
num = 20  # determins the size of the plane 100
N_X = 1100  # length of 100 direction
N_X_D = 190  # location of the detector
N_Y = num
N_Z = num
# NXNYNZ means the number of cells, and 2 means 2 atoms per unit cell. This calculates the total number of atoms in the entire cell
n_atoms = N_X * N_Y * N_Z * 2
n_arrange = N_X_D * N_Y * N_Z * 2  # the numbers of atoms that can be modified
# atoms located farther than NXD can't be modified

# Set here is the density rate of vacancies/precipitates
RATE = 0.000  # vacancy density rate #if RATE = 1, material will be nothing
n_vacancy = math.floor(n_arrange * RATE)
cut_off = 5.6

# in order to get the whole simulation box to be non-cyclic, we must put vacuum zone at the both end in x-direction.
#the thickness of the vaccum is set as the variable "indent", seen below:
indent = 3  # indent should not be zero

print(n_vacancy)

# Set here are the fundamental arrangement of atoms; as for bcc crystal, put atoms at (0,0,0), (1/2,1/2,1/2), and the corresponding location of each unit cell 
vert = [0., 0., 0.]
surf1 = [0.5*LAT_CONST, 0.5*LAT_CONST, 0.5*LAT_CONST]
#surf2 = [0.5*LAT_CONST,0.5*LAT_CONST,0.]
#surf3 = [0.5*LAT_CONST,0.,0.5*LAT_CONST]

size_x = (DIST_X * (N_X+indent*2)) * (1.0 + STRAIN_X)
size_y = (DIST_Y * N_Y) * (1.0 + STRAIN_Y)
size_z = (DIST_Z * N_Z) * (1.0 + STRAIN_Z)


# the function below creates the list of atoms in perfect crystal
# Each atom has the list of:
    # ID, coordinates in order of x, y, z, groupID, present/absent binary, distance over the cyclic boundaries
# from the left to the right, respectively.
# group ID
# 1:Most of the atoms located at the left of the detector
# 2:atoms at the source
# 3:atoms at the detector
# 4:Most of the atoms located at the right of the detector
# 5:atoms at the right end of the simulation box
# 6:Cu atoms

# atom_number x y z grouping_num TF inv_y inv_z
perfect = np.zeros((n_atoms, 8))


def make_perfect():

    count = 0
    # make perfect lattice
    for ix in range(indent, N_X+indent):
        for iy in range(N_Y):
            for iz in range(N_Z):
                if ix == indent:
                    perfect[count] = [count+1, vert[0]+ix*LAT_CONST, vert[1]+iy*LAT_CONST, vert[2] +
                                      iz*LAT_CONST, 2, 1, vert[1]+(iy-N_Y)*LAT_CONST, vert[2]+(iz-N_Z)*LAT_CONST]
                    count += 1
                    perfect[count] = [count+1, surf1[0]+ix*LAT_CONST, surf1[1]+iy*LAT_CONST, surf1[2] +
                                      iz*LAT_CONST, 2, 1, surf1[1]+(iy-N_Y)*LAT_CONST, surf1[2]+(iz-N_Z)*LAT_CONST]
                    count += 1
                elif ix == N_X_D+indent-1:
                    perfect[count] = [count+1, vert[0]+ix*LAT_CONST, vert[1]+iy*LAT_CONST, vert[2] +
                                      iz*LAT_CONST, 3, 1, vert[1]+(iy-N_Y)*LAT_CONST, vert[2]+(iz-N_Z)*LAT_CONST]
                    count += 1
                    perfect[count] = [count+1, surf1[0]+ix*LAT_CONST, surf1[1]+iy*LAT_CONST, surf1[2] +
                                      iz*LAT_CONST, 3, 1, surf1[1]+(iy-N_Y)*LAT_CONST, surf1[2]+(iz-N_Z)*LAT_CONST]
                    count += 1
                elif ix < N_X_D+indent-1:
                    perfect[count] = [count+1, vert[0]+ix*LAT_CONST, vert[1]+iy*LAT_CONST, vert[2] +
                                      iz*LAT_CONST, 1, 1, vert[1]+(iy-N_Y)*LAT_CONST, vert[2]+(iz-N_Z)*LAT_CONST]
                    count += 1
                    perfect[count] = [count+1, surf1[0]+ix*LAT_CONST, surf1[1]+iy*LAT_CONST, surf1[2] +
                                      iz*LAT_CONST, 1, 1, surf1[1]+(iy-N_Y)*LAT_CONST, surf1[2]+(iz-N_Z)*LAT_CONST]
                    count += 1
                elif ix == N_X+indent-1:
                    perfect[count] = [count+1, vert[0]+ix*LAT_CONST, vert[1]+iy*LAT_CONST, vert[2] +
                                      iz*LAT_CONST, 5, 1, vert[1]+(iy-N_Y)*LAT_CONST, vert[2]+(iz-N_Z)*LAT_CONST]
                    count += 1
                    perfect[count] = [count+1, surf1[0]+ix*LAT_CONST, surf1[1]+iy*LAT_CONST, surf1[2] +
                                      iz*LAT_CONST, 5, 1, surf1[1]+(iy-N_Y)*LAT_CONST, surf1[2]+(iz-N_Z)*LAT_CONST]
                    count += 1
                else:
                    perfect[count] = [count+1, vert[0]+ix*LAT_CONST, vert[1]+iy*LAT_CONST, vert[2] +
                                      iz*LAT_CONST, 4, 1, vert[1]+(iy-N_Y)*LAT_CONST, vert[2]+(iz-N_Z)*LAT_CONST]
                    count += 1
                    perfect[count] = [count+1, surf1[0]+ix*LAT_CONST, surf1[1]+iy*LAT_CONST, surf1[2] +
                                      iz*LAT_CONST, 4, 1, surf1[1]+(iy-N_Y)*LAT_CONST, surf1[2]+(iz-N_Z)*LAT_CONST]
                    count += 1


def center(x):
    if x % 2 == 0:
        center = x/2
    else:
        center = (x-1)/2
    return center


# the function below introduces vacancies or precipitates
# in case of vacancies, this will delete the atoms
# in case of precipitates, this will change the group ID of some atoms into 6
def arrange_lat(R, Vac_Cu):
    del_num = 0
    klist = [1]  # absolutely not be void in this position
    for j in range(n_atoms):  # any big value is ok
        if del_num >= n_vacancy:
            break
        p = random.randint(1, n_arrange)
        if R == 0:
            if perfect[p, 4] == 1 and perfect[p, 1] >= (cut_off+indent*LAT_CONST) and perfect[p, 1] <= ((N_X_D+indent)*LAT_CONST - (cut_off)):
                for k in range(len(klist)):
                    if perfect[klist[k], 5] == 0 and Vac_Cu == 1:          # if Vac [k,5] ==0
                        if neighbor(p, klist[k], R, cut_off) == 1:
                            break
                    if perfect[klist[k], 4] == 6 and Vac_Cu == 2:          # if Cu [k,4] ==5
                        if neighbor(p, klist[k], R, cut_off) == 1:
                            break
                    if k == len(klist)-1:
                        if Vac_Cu == 1:
                            perfect[p, 5] = 0           # if Vac[p,5] ==0
                            del_num += 1
                            klist.append(p)
                            print(len(klist)-1)
                            # print(del_num)
                        if Vac_Cu == 2:
                            perfect[p, 4] = 6           # if Cu [p,4] ==5
                            del_num += 1
                            # print(del_num)
                            klist.append(p)
                            print(len(klist)-1)
        else:  # means??
            #delta = 0.05
            delta = 0.5
            if perfect[p, 4] == 1 and perfect[p, 1] >= (N_X_D*(0.5-delta)+indent)*LAT_CONST+R+cut_off and perfect[p, 1] <= (N_X_D*(0.5+delta)+indent)*LAT_CONST-(R+cut_off):
                for k in range(len(klist)):
                    if perfect[klist[k], 5] == 0 and Vac_Cu == 1:          # if Vac [k,5] ==0
                        if neighbor(p, klist[k], R, cut_off) == 1:
                            break
                    if perfect[klist[k], 4] == 6 and Vac_Cu == 2:          # if Cu [k,4] ==5
                        if neighbor(p, klist[k], R, cut_off) == 1:
                            break
                    if k == len(klist)-1:
                        for v in range(n_arrange):
                            if neighbor(p, v, R, 0) == 1:
                                if Vac_Cu == 1:
                                    perfect[v, 5] = 0  # if Vac [v,5] ==0
                                    del_num += 1
                                    klist.append(v)
                                    # print(klist[del_num-1])
                                    # print(del_num)
                                    print(len(klist)-1)
                                if Vac_Cu == 2:
                                    perfect[v, 4] = 6  # if Cu [v,4] ==5
                                    del_num += 1
                                    # print(del_num)
                                    klist.append(v)
                                    print(len(klist)-1)

    name(Vac_Cu+4)


def name(type_num):
    f = open("Fe.lat", "w")
    f.write("# Fe\n")
    f.write("\n")
    #f.write("%d atoms\n" % perfect_num)
    f.write("\n")
    f.write("%d atom types\n" % type_num)
    f.write("0 %12.6f xlo xhi\n" % size_x)
    f.write("0 %12.6f ylo yhi\n" % size_y)
    f.write("0 %12.6f zlo zhi\n" % size_z)
    f.write("\n")
    f.write("Atoms\n")
    f.write("\n")
    lat = 0
    for i in range(n_atoms):
        if perfect[i, 5] == 1:
            lat += 1
            f.write("%12d %6d %12.6f %12.6f %12.6f\n" % (
                lat, perfect[i, 4], perfect[i, 1], perfect[i, 2], perfect[i, 3]))
    f.close()

    with open("Fe.lat") as f:
        l = f.readlines()

    l.insert(2, "%d atoms" % lat)
    with open("Fe.lat", mode='w') as f:
        f.writelines(l)



def input_empty():
    f = open("in.lammps", "a")
    for jx in range(indent):
        for jy in range(num):
            for jz in range(num):
                f.write("create_atoms %d single %12.6f %12.6f %12.6f\n" % (
                    4, vert[0]+jx*LAT_CONST, vert[1]+jy*LAT_CONST, vert[2]+jz*LAT_CONST))
                f.write("create_atoms %d single %12.6f %12.6f %12.6f\n" % (
                    4, surf1[0]+jx*LAT_CONST, surf1[1]+jy*LAT_CONST, surf1[2]+jz*LAT_CONST))

    for jx in range(N_X+indent, N_X+2*indent):
        for jy in range(num):
            for jz in range(num):
                f.write("create_atoms %d single %12.6f %12.6f %12.6f\n" % (
                    4, vert[0]+jx*LAT_CONST, vert[1]+jy*LAT_CONST, vert[2]+jz*LAT_CONST))
                f.write("create_atoms %d single %12.6f %12.6f %12.6f\n" % (
                    4, surf1[0]+jx*LAT_CONST, surf1[1]+jy*LAT_CONST, surf1[2]+jz*LAT_CONST))
    f.write("write_restart restart.equil\n")
    f.close()



def inverse_lat():
    f = open("inverse.lat", "w")
    lat = 0
    for i in range(perfect_num):
        if perfect[i, 5] == 0:  # if Cu [i,4] ==5
            lat += 1
            f.write("%12d %6d %12.6f %12.6f %12.6f\n" % (
                lat, perfect[i, 4], perfect[i, 1], perfect[i, 2], perfect[i, 3]))
    f.close()


# This function determines whether the two atoms are closer than the cutoff length
def neighbor(p, k, R, cut):
    if distance(p, k, 1)**2 + distance(p, k, 2)**2 + distance(p, k, 3)**2 <= (cut+R)**2:
        return 1
    else:
        return 0

# This function culculates the distance


def distance(p, k, a):
    if a == 1:
        return abs(perfect[p, 1]-perfect[k, 1])
    elif a == 2:
        b1 = abs(perfect[p, 2]-perfect[k, 2])
        b2 = abs(perfect[p, 2]-perfect[k, 6])
        b3 = abs(perfect[p, 6]-perfect[k, 2])
        # theorically b4 is not needed
        return min(b1, b2, b3)

    else:
        b1 = abs(perfect[p, 3]-perfect[k, 3])
        b2 = abs(perfect[p, 3]-perfect[k, 7])
        b3 = abs(perfect[p, 7]-perfect[k, 3])
        # theorically b4 is not needed
        return min(b1, b2, b3)


# mainpart
if __name__ == "__main__":
    argvs = sys.argv
    make_perfect()
    arrange_lat(3, 1)  # (R ,type) type= vacancy:1 , Cu: 2
    # input_empty()
    # and RATE or num (=N_X,N_Y)are all parameters
    # R unit is angstrome


# moving test
# check num of arrage
# if you do not need, delete 'dump'
#
# Use OVITO to check molphorogy
# FFT configuration check
# arraving time should be defined consistnatly
# change x_detec, A2, Time to get beta
