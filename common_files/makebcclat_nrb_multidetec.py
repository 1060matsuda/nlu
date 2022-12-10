# %%
# -*- coding:utf-8 -*-
import sys
import math
import numpy as np
import random
import yaml

# lattice constant a0 and strain settings
LAT_CONST = 2.855324  # A
POISSON = 0.293
STRAIN_X = 0.000
STRAIN_Y = -POISSON * STRAIN_X
STRAIN_Z = -POISSON * STRAIN_X

DIST_X = LAT_CONST
DIST_Y = LAT_CONST
DIST_Z = LAT_CONST

with open("config.yaml", "r") as yml:
    config = yaml.safe_load(yml)

# the size of the whole simulation cell
# N_Y and N_Z represents the length of the cell (unit: lattice constant)
num = 20  # determins the size of the plane 100
N_Y = num
N_Z = num

frequency_high = float(config["f1"])
frequency_low = float(config["f2"])  # GHz
"""[GHz]"""
cycles = 6
"""Specify how many cycles of low-frequency-wave you want to input."""

# detector location
N_X_D = int(1001)
"""Detector's position is determined based on the HIGHER frequency. 
This is because f+ and f- (which we want to detect) is around the higher frequency.
EX) 500GHz +- 10GHZ -> 490GHz and 510GHz"""

print("N_X_D")
print(N_X_D)
N_BUF = int(N_X_D-1)*3 +1

# The cell length including the buffer layer
#N_BUF = 700

# The cell length before(excluding) the buffer layer
N_X = int(N_X_D-1)*2 +1


# N_BIG is set when you place atomic more layers right beside the buffer kayer.
# Currently not used; N_BIG-N_BUF=0 is required.
N_BIG = N_BUF

detecsNum = int(20)
XDArray = np.zeros(detecsNum)
XD_interval = math.ceil(N_X_D/(detecsNum+1))

#N_X_D0 =int(N_X*(1-2**(-0)))
for i in range(detecsNum):
    XDArray[i] = int(XD_interval*(i+1))


# The whole number of atoms in the simulation setup
#n_atoms = N_X * N_Y *N_Z * 2
n_atoms = N_BIG * N_Y * N_Z * 2
print("n_atoms")
print(n_atoms)
# the number of atoms in the region where defects are introduced
n_arrange = (N_X_D-1) * N_Y * N_Z * 2
n_edge = N_X * N_Y * N_Z * 2
# Set here is the density of vacancies/precipitates
if "defects_rate" in config:
    # vacancy density rate. If RATE = 1, material will be nothing, or whole Cu
    RATE = float(config["defects_rate"])
else:
    RATE = 0.0
print("defects_rate is:")
print(RATE)

# %%
# Floor function. If n_arrange*RATE == 52.6, then it retuns 52
n_vacancy = math.floor(n_arrange * RATE)
cut_off = 5.6

# in order to get the whole simulation box to be non-cyclic, we must put vacuum zone at the both end in x-direction.
# the thickness of the vaccum is set as the variable "indent", seen below:
indent = 500  # indent should not be zero

print(n_vacancy)

# Set here are the fundamental arrangement of atoms; as for bcc crystal, put atoms at (0,0,0), (1/2,1/2,1/2), and the corresponding location of each unit cell
vert = [0., 0., 0.]
surf1 = [0.5*LAT_CONST, 0.5*LAT_CONST, 0.5*LAT_CONST]
#surf2 = [0.5*LAT_CONST,0.5*LAT_CONST,0.]
#surf3 = [0.5*LAT_CONST,0.,0.5*LAT_CONST]

#size_x = (DIST_X * (N_X+indent*2)) * (1.0 + STRAIN_X)
size_x = (DIST_X * (N_BIG+indent)) * (1.0 + STRAIN_X)
size_y = (DIST_Y * N_Y) * (1.0 + STRAIN_Y)
size_z = (DIST_Z * N_Z) * (1.0 + STRAIN_Z)


# the function below creates the list of atoms in perfect crystal
# Each atom has the list of:
# atom ID, coordinates in order of x, y, z, groupID, absent(0)/present(1) binary, distance over the cyclic boundaries
# group ID
# 1:Most of the atoms located at the left of the detector. These can be modified into porecipitates or vacancies
# 2:atoms at the source
# 3:atoms at the detector
# 4:Most of the atoms located at the right of the detector. These can't be modified
# 5:atoms at the right end of the simulation box (called "edge"). It depends on the crystalline structure whether to denote 5 or 6. In this program only 6 is used.
# 7,8,9:atoms in the region of the right side of the edge (=buffer layer).It depends on the crystalline structure whether to denote 7, 8, or 9. In this program, three numbers are treated equally.
# detailed information (I didn't understand well though) about 9 and 10 is described in the original makebcclat.py.

# atom_number x y z grouping_num TF inv_y inv_z
perfect = np.zeros((n_atoms, 8))


def make_perfect():
    """the function below creates the list of atoms in perfect crystal.
    Each atom has the list of:
    atom ID, coordinates in order of x, y, z, groupID, absent(0)/present(1) binary, distance over the cyclic boundaries.
    Group IDs reference:
    1:Most of the atoms located at the left of the detector. These can be modified into porecipitates or vacancies
    2:atoms at the source
    3:atoms at the detector
    4:Most of the atoms located at the right of the detector. These can't be modified
    5:atoms at the right end of the simulation box (called "edge"). It depends on the crystalline structure whether to denote 5 or 6. In this program only 6 is used.
    7,8,9:atoms in the region of the right side of the edge (=buffer layer).It depends on the crystalline structure whether to denote 7, 8, or 9. In this program, three numbers are treated equally.
    detailed information (I didn't understand well though) about 9 and 10 is described in the original makebcclat.py."""
    count = 0
    # make perfect lattice
    for ix in range(N_BIG):
        for iy in range(N_Y):
            for iz in range(N_Z):
                if ix == 0:
                    perfect[count] = [count+1, vert[0]+ix*LAT_CONST, vert[1]+iy*LAT_CONST, vert[2] +
                                      iz*LAT_CONST, 2, 1, vert[1]+(iy-N_Y)*LAT_CONST, vert[2]+(iz-N_Z)*LAT_CONST]
                    count += 1
                    perfect[count] = [count+1, surf1[0]+ix*LAT_CONST, surf1[1]+iy*LAT_CONST, surf1[2] +
                                      iz*LAT_CONST, 2, 1, surf1[1]+(iy-N_Y)*LAT_CONST, surf1[2]+(iz-N_Z)*LAT_CONST]
                    count += 1
                elif ix < N_X_D:
                    perfect[count] = [count+1, vert[0]+ix*LAT_CONST, vert[1]+iy*LAT_CONST, vert[2] +
                                      iz*LAT_CONST, 1, 1, vert[1]+(iy-N_Y)*LAT_CONST, vert[2]+(iz-N_Z)*LAT_CONST]
                    count += 1
                    perfect[count] = [count+1, surf1[0]+ix*LAT_CONST, surf1[1]+iy*LAT_CONST, surf1[2] +
                                      iz*LAT_CONST, 1, 1, surf1[1]+(iy-N_Y)*LAT_CONST, surf1[2]+(iz-N_Z)*LAT_CONST]
                    count += 1
                elif ix == N_X_D:
                    perfect[count] = [count+1, vert[0]+ix*LAT_CONST, vert[1]+iy*LAT_CONST, vert[2] +
                                      iz*LAT_CONST, 3, 1, vert[1]+(iy-N_Y)*LAT_CONST, vert[2]+(iz-N_Z)*LAT_CONST]
                    count += 1
                    perfect[count] = [count+1, surf1[0]+ix*LAT_CONST, surf1[1]+iy*LAT_CONST, surf1[2] +
                                      iz*LAT_CONST, 3, 1, surf1[1]+(iy-N_Y)*LAT_CONST, surf1[2]+(iz-N_Z)*LAT_CONST]
                    count += 1
                elif ix < N_X:
                    perfect[count] = [count+1, vert[0]+ix*LAT_CONST, vert[1]+iy*LAT_CONST, vert[2] +
                                      iz*LAT_CONST, 4, 1, vert[1]+(iy-N_Y)*LAT_CONST, vert[2]+(iz-N_Z)*LAT_CONST]
                    count += 1
                    perfect[count] = [count+1, surf1[0]+ix*LAT_CONST, surf1[1]+iy*LAT_CONST, surf1[2] +
                                      iz*LAT_CONST, 4, 1, surf1[1]+(iy-N_Y)*LAT_CONST, surf1[2]+(iz-N_Z)*LAT_CONST]
                    count += 1
                elif ix == N_X:
                    perfect[count] = [count+1, vert[0]+ix*LAT_CONST, vert[1]+iy*LAT_CONST, vert[2] +
                                      iz*LAT_CONST, 5, 1, vert[1]+(iy-N_Y)*LAT_CONST, vert[2]+(iz-N_Z)*LAT_CONST]
                    count += 1
                    perfect[count] = [count+1, surf1[0]+ix*LAT_CONST, surf1[1]+iy*LAT_CONST, surf1[2] +
                                      iz*LAT_CONST, 6, 1, surf1[1]+(iy-N_Y)*LAT_CONST, surf1[2]+(iz-N_Z)*LAT_CONST]
                    count += 1
                elif ix == N_X+1:
                    perfect[count] = [count+1, vert[0]+ix*LAT_CONST, vert[1]+iy*LAT_CONST, vert[2] +
                                      iz*LAT_CONST, 7, 1, vert[1]+(iy-N_Y)*LAT_CONST, vert[2]+(iz-N_Z)*LAT_CONST]
                    count += 1
                    perfect[count] = [count+1, surf1[0]+ix*LAT_CONST, surf1[1]+iy*LAT_CONST, surf1[2] +
                                      iz*LAT_CONST, 8, 1, surf1[1]+(iy-N_Y)*LAT_CONST, surf1[2]+(iz-N_Z)*LAT_CONST]
                    count += 1
                elif ix <= N_BUF:
                    perfect[count] = [count+1, vert[0]+ix*LAT_CONST, vert[1]+iy*LAT_CONST, vert[2] +
                                      iz*LAT_CONST, 9, 1, vert[1]+(iy-N_Y)*LAT_CONST, vert[2]+(iz-N_Z)*LAT_CONST]
                    count += 1
                    perfect[count] = [count+1, surf1[0]+ix*LAT_CONST, surf1[1]+iy*LAT_CONST, surf1[2] +
                                      iz*LAT_CONST, 9, 1, surf1[1]+(iy-N_Y)*LAT_CONST, surf1[2]+(iz-N_Z)*LAT_CONST]
                    count += 1
                else:
                    perfect[count] = [count+1, vert[0]+ix*LAT_CONST, vert[1]+iy*LAT_CONST, vert[2] +
                                      iz*LAT_CONST, 10, 1, vert[1]+(iy-N_Y)*LAT_CONST, vert[2]+(iz-N_Z)*LAT_CONST]
                    count += 1
                    perfect[count] = [count+1, surf1[0]+ix*LAT_CONST, surf1[1]+iy*LAT_CONST, surf1[2] +
                                      iz*LAT_CONST, 10, 1, surf1[1]+(iy-N_Y)*LAT_CONST, surf1[2]+(iz-N_Z)*LAT_CONST]
                    count += 1
                for arrayNum in range(len(XDArray)):
                    if int(ix) == XDArray[arrayNum]:
                        perfect[count-2][4]=arrayNum+10
                        perfect[count-1][4]=arrayNum+10

def center(x):
    if x % 2 == 0:
        center = x/2
    else:
        center = (x-1)/2
    return center


# the function below introduces vacancies or precipitates
# in case of vacancies, this will delete the atoms
# in case of precipitates, this will change the group ID of some atoms into 10
# If Vac_Cu == 1, then it'll introduce vacancies or voids. If Vac_Cu == 2, then it'll introduce precipitates
def arrange_lat(R, Vac_Cu):
    del_num = 0
    klist = [1]  # absolutely not be void in this position
    # When seed is specified as constant, output lattice will always be the same.
    # np.random.seed(seed=32)
    plist = np.random.randint(1, n_arrange+1, n_atoms)
    for j in range(n_atoms):  # any big value is ok
        if del_num >= n_vacancy:
            break
        #p = random.randint(1,n_arrange+1)
        p = plist[j]
        if R == 0:
            # Note that, currently, the latter condition means nothing.
            if perfect[p, 4] == 1 and \
               ((perfect[p, 1] >= (0.5)*LAT_CONST+cut_off and perfect[p, 1] <= (N_X_D)*LAT_CONST-cut_off)
               or
               (perfect[p, 1] >= (N_X_D+0.5)*LAT_CONST+cut_off and perfect[p, 1] <= (N_X)*LAT_CONST-cut_off)):
                for k in range(len(klist)):
                    # if atom k is Vac, perfect[k,5] should be 0
                    if perfect[klist[k], 5] == 0 and Vac_Cu == 1:
                        if neighbor(p, klist[k], R, cut_off) == 1:
                            break
                    # if atom k is Cu, perfect[k,4] should be 10
                    if perfect[klist[k], 4] == 10 and Vac_Cu == 2:
                        if neighbor(p, klist[k], R, cut_off) == 1:
                            break
                    if k == len(klist)-1:
                        if Vac_Cu == 1:
                            perfect[p, 5] = 0           # if Vac[p,5] ==0
                            del_num += 1
                            klist.append(p)
                            print("Atom "+str(p)+"has been deleted") #This video has been deleted #Pro-wanker
                            # print(len(klist)-1)
                            # print(del_num)
                        if Vac_Cu == 2:
                            perfect[p, 4] = 10           # if Cu [p,4] ==10
                            del_num += 1
                            # print(del_num)
                            klist.append(p)
                            # print(len(klist)-1)
        else:
            # Supposing the radius is more than 2.8A. If the cutoff is set to 2R(2*radius), then particles are scattered somewhat well. Currently virtually unused.
            #delta = 0.05
            #delta = 0.5
            # if perfect[p,4] ==1 and perfect[p,1] >= (N_X_D*(0.5-delta)+indent)*LAT_CONST+R+cut_off and perfect[p,1] <= (N_X_D*(0.5+delta)+indent)*LAT_CONST-(R+cut_off):]
            if perfect[p, 4] == 1 and \
               ((perfect[p, 1] >= (0.5)*LAT_CONST+cut_off+R and perfect[p, 1] <= (N_X_D)*LAT_CONST-cut_off-R)
               or
               (perfect[p, 1] >= (N_X_D+0.5)*LAT_CONST+cut_off+R and perfect[p, 1] <= (N_X)*LAT_CONST-cut_off-R)):
                for k in range(len(klist)):
                    if perfect[klist[k], 5] == 0 and Vac_Cu == 1:          # if Vac, then [k,5] ==0
                        if neighbor(p, klist[k], R, 2*R) == 1:
                            break
                    if perfect[klist[k], 4] == 10 and Vac_Cu == 2:          # if Cu, then [k,4] ==10
                        if neighbor(p, klist[k], R, 2*R) == 1:
                            break
                    if k == len(klist)-1:
                        for v in range(n_arrange):
                            if neighbor(p, v, R, 0) == 1:
                                if Vac_Cu == 1:
                                    perfect[v, 5] = 0  # if Vac, [v,5] ==0
                                    del_num += 1
                                    klist.append(v)
                                    # print(klist[del_num-1])
                                    # print(del_num)
                                    # print(len(klist)-1)
                                if Vac_Cu == 2:
                                    perfect[v, 4] = 10  # if Cu, [v,4] ==10
                                    del_num += 1
                                    # print(del_num)
                                    klist.append(v)
                                    # print(len(klist)-1)
    print("Number of deleted/modified atoms is")
    print(del_num)
    print("point-defect-equivalent density:")
    print(del_num/n_arrange)
    name(Vac_Cu+8+detecsNum)

# create Fe.lat


def name(type_num):
    f = open("Fe.lat", "w")
    f.write("# Fe\n")
    f.write("\n")
    #f.write("%d atoms\n" % perfect_num)
    f.write("\n")
    f.write("%d atom types\n" % type_num)
    f.write("%12.6f %12.6f xlo xhi\n" % (-DIST_X*indent*(1+STRAIN_X), size_x))
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
                i+1, perfect[i, 4], perfect[i, 1], perfect[i, 2], perfect[i, 3]))
    f.close()

    with open("Fe.lat") as f:
        l = f.readlines()

    l.insert(2, "%d atoms" % lat)
    with open("Fe.lat", mode='w') as f:
        f.writelines(l)

    add(type_num)

# create add.lat, which tells you the initial position of buffer layer atoms.


def add(type_num):
    f = open("add.lat", "w")
    f.write("# Fe\n")
    f.write("\n")
    #f.write("%d atoms\n" % perfect_num)
    f.write("\n")
    f.write("%d atom types\n" % type_num)
    f.write("%12.6f %12.6f xlo xhi\n" % (-DIST_X*indent*(1+STRAIN_X), size_x))
    f.write("0 %12.6f ylo yhi\n" % size_y)
    f.write("0 %12.6f zlo zhi\n" % size_z)
    f.write("\n")
    f.write("Atoms\n")
    f.write("\n")
    lat1 =0
    lat2 =n_edge
    #lat = n_atoms
    for i in range(n_atoms):
        if perfect[i, 5] == 1 and (perfect[i, 4] == 9 or perfect[i, 4] == 7 or perfect[i, 4] == 8):
            lat1 +=1
            lat2 +=1
            f.write("%12d %6d %12.6f %12.6f %12.6f\n" % (
                i+1, perfect[i, 4], perfect[i, 1], perfect[i, 2], perfect[i, 3]))
    f.close()

    with open("add.lat") as f:
        l = f.readlines()

    l.insert(2, "%d atoms" % lat1)
    with open("add.lat", mode='w') as f:
        f.writelines(l)


# currently not used. Detail is described in the original file.
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


# currently not used. Detail is described in the original file.
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
# if two atoms are neighbors, it means they don't belong to different clusters; in other words they are very nearby)
def neighbor(p, k, R, cut):
    if distance(p, k, 1)**2 + distance(p, k, 2)**2 + distance(p, k, 3)**2 <= (cut+R)**2:
        return 1
    else:
        return 0

# This function culculates the distance in x / y / z direction


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
defect_radius = float(config["defects_radius"])

if __name__ == "__main__":
    argvs = sys.argv
    make_perfect()
    arrange_lat(defect_radius, 1)  # (R ,type) type= vacancy:1 , Cu: 2

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
