# -*- coding:utf-8 -*-
import sys
import math
import numpy as np
import random

# lattice constant a0 and strain settings
LAT_CONST = 2.855324 
POISSON = 0.293
STRAIN_X = 0.000
STRAIN_Y = -POISSON * STRAIN_X
STRAIN_Z = -POISSON * STRAIN_X

DIST_X = LAT_CONST 
DIST_Y = LAT_CONST 
DIST_Z = LAT_CONST 

# the size of the whole simulation cell
#N_Y and N_Z represents the length of the cell (unit: lattice constant)
num =  20 # determins the size of the plane 100
N_Y = num
N_Z = num

#The cell length before(excluding) the buffer layer
N_X =  301

#The cell length including the buffer layer
N_BUF = 451

#detector location
N_X_D = 151

#N_BIG is set when you place atomic more layers right beside the buffer kayer.
#Currently not used; N_BIG-N_BUF=0 is required. 
N_BIG = 451

#The whole number of atoms in the simulation setup
#n_atoms = N_X * N_Y *N_Z * 2
n_atoms = N_BIG * N_Y *N_Z * 2

# the number of atoms in the region where defects are introduced
n_arrange = (N_X_D-1) * N_Y * N_Z * 2 

# Set here is the density of vacancies/precipitates
RATE = 0.000 # vacancy density rate. If RATE = 1, material will be nothing, or whole Cu
n_vacancy = math.floor(n_arrange * RATE) #Floor function. If n_arrange*RATE == 52.6, then it retuns 52
cut_off = 5.6

# in order to get the whole simulation box to be non-cyclic, we must put vacuum zone at the both end in x-direction.
# the thickness of the vaccum is set as the variable "indent", seen below:
indent = 5  #indent should not be zero

print(n_vacancy)

#Set here are the fundamental arrangement of atoms; as for bcc crystal, put atoms at (0,0,0), (1/2,1/2,1/2), and the corresponding location of each unit cell
vert = [0.,0.,0.]
surf1 = [0.5*LAT_CONST,0.5*LAT_CONST,0.5*LAT_CONST]
#surf2 = [0.5*LAT_CONST,0.5*LAT_CONST,0.]
#surf3 = [0.5*LAT_CONST,0.,0.5*LAT_CONST]

#size_x = (DIST_X * (N_X+indent*2)) * (1.0 + STRAIN_X)
size_x = (DIST_X * (N_X_D+indent*2)) * (1.0 + STRAIN_X)
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

perfect =np.zeros((n_atoms, 8)) # atom_number x y z grouping_num TF inv_y inv_z

def make_perfect():

    count = 0
    # make perfect lattice
    for ix in range(indent+1,N_BIG+indent+1):
        for iy in range(N_Y):
            for iz in range(N_Z): 
                if ix==indent+1:                   
                    perfect[count] = [count+1 ,vert[0]+ix*LAT_CONST ,vert[1]+iy*LAT_CONST ,vert[2]+iz*LAT_CONST ,2 , 1,vert[1]+(iy-N_Y)*LAT_CONST,vert[2]+(iz-N_Z)*LAT_CONST]
                    count +=1      
                    perfect[count] = [count+1 ,surf1[0]+ix*LAT_CONST ,surf1[1]+iy*LAT_CONST ,surf1[2]+iz*LAT_CONST ,2 , 1,surf1[1]+(iy-N_Y)*LAT_CONST,surf1[2]+(iz-N_Z)*LAT_CONST]
                    count +=1
                elif ix < N_X_D+indent:
                    perfect[count] = [count+1 ,vert[0]+ix*LAT_CONST ,vert[1]+iy*LAT_CONST ,vert[2]+iz*LAT_CONST ,1 , 1,vert[1]+(iy-N_Y)*LAT_CONST,vert[2]+(iz-N_Z)*LAT_CONST]
                    count +=  1
                    perfect[count] = [count+1 ,surf1[0]+ix*LAT_CONST ,surf1[1]+iy*LAT_CONST ,surf1[2]+iz*LAT_CONST ,1 , 1,surf1[1]+(iy-N_Y)*LAT_CONST,surf1[2]+(iz-N_Z)*LAT_CONST]
                    count +=  1
                elif ix==N_X_D+indent:
                    perfect[count] = [count+1 ,vert[0]+ix*LAT_CONST ,vert[1]+iy*LAT_CONST ,vert[2]+iz*LAT_CONST ,3 , 1,vert[1]+(iy-N_Y)*LAT_CONST,vert[2]+(iz-N_Z)*LAT_CONST]
                    count += 1 
                    perfect[count] = [count+1 ,surf1[0]+ix*LAT_CONST ,surf1[1]+iy*LAT_CONST ,surf1[2]+iz*LAT_CONST ,3 , 1,surf1[1]+(iy-N_Y)*LAT_CONST,surf1[2]+(iz-N_Z)*LAT_CONST]
                    count +=  1
                elif ix < N_X+indent:
                    perfect[count] = [count+1 ,vert[0]+ix*LAT_CONST ,vert[1]+iy*LAT_CONST ,vert[2]+iz*LAT_CONST ,4 , 0,vert[1]+(iy-N_Y)*LAT_CONST,vert[2]+(iz-N_Z)*LAT_CONST]
                    count +=  1
                    perfect[count] = [count+1 ,surf1[0]+ix*LAT_CONST ,surf1[1]+iy*LAT_CONST ,surf1[2]+iz*LAT_CONST ,4 , 0,surf1[1]+(iy-N_Y)*LAT_CONST,surf1[2]+(iz-N_Z)*LAT_CONST]
                    count +=  1
                elif ix == N_X+indent:
                    perfect[count] = [count+1 ,vert[0]+ix*LAT_CONST ,vert[1]+iy*LAT_CONST ,vert[2]+iz*LAT_CONST ,5 , 0,vert[1]+(iy-N_Y)*LAT_CONST,vert[2]+(iz-N_Z)*LAT_CONST]
                    count +=  1
                    perfect[count] = [count+1 ,surf1[0]+ix*LAT_CONST ,surf1[1]+iy*LAT_CONST ,surf1[2]+iz*LAT_CONST ,6 , 0,surf1[1]+(iy-N_Y)*LAT_CONST,surf1[2]+(iz-N_Z)*LAT_CONST]
                    count +=  1
                elif ix == N_X+indent+1:
                    perfect[count] = [count+1 ,vert[0]+ix*LAT_CONST ,vert[1]+iy*LAT_CONST ,vert[2]+iz*LAT_CONST ,7 , 0,vert[1]+(iy-N_Y)*LAT_CONST,vert[2]+(iz-N_Z)*LAT_CONST]
                    count +=  1
                    perfect[count] = [count+1 ,surf1[0]+ix*LAT_CONST ,surf1[1]+iy*LAT_CONST ,surf1[2]+iz*LAT_CONST ,8 , 0,surf1[1]+(iy-N_Y)*LAT_CONST,surf1[2]+(iz-N_Z)*LAT_CONST]
                    count +=  1
                elif ix <= N_BUF+indent:
                    perfect[count] = [count+1 ,vert[0]+ix*LAT_CONST ,vert[1]+iy*LAT_CONST ,vert[2]+iz*LAT_CONST ,9 , 0,vert[1]+(iy-N_Y)*LAT_CONST,vert[2]+(iz-N_Z)*LAT_CONST]
                    count +=  1
                    perfect[count] = [count+1 ,surf1[0]+ix*LAT_CONST ,surf1[1]+iy*LAT_CONST ,surf1[2]+iz*LAT_CONST ,9 , 0,surf1[1]+(iy-N_Y)*LAT_CONST,surf1[2]+(iz-N_Z)*LAT_CONST]
                    count +=  1
                else:
                    perfect[count] = [count+1 ,vert[0]+ix*LAT_CONST ,vert[1]+iy*LAT_CONST ,vert[2]+iz*LAT_CONST ,10 , 0,vert[1]+(iy-N_Y)*LAT_CONST,vert[2]+(iz-N_Z)*LAT_CONST]
                    count +=  1
                    perfect[count] = [count+1 ,surf1[0]+ix*LAT_CONST ,surf1[1]+iy*LAT_CONST ,surf1[2]+iz*LAT_CONST ,10 , 0,surf1[1]+(iy-N_Y)*LAT_CONST,surf1[2]+(iz-N_Z)*LAT_CONST]
                    count +=  1

def center(x):
    if x%2==0:
        center = x/2
    else:
        center = (x-1)/2
    return center


# the function below introduces vacancies or precipitates
# in case of vacancies, this will delete the atoms
# in case of precipitates, this will change the group ID of some atoms into 10
# If Vac_Cu == 1, then it'll introduce vacancies or voids. If Vac_Cu == 2, then it'll introduce precipitates
def arrange_lat(R,Vac_Cu):
    del_num = 0   
    klist = [1] # absolutely not be void in this position
    np.random.seed(seed=32)
    plist = np.random.randint(1,n_arrange+1,n_atoms)
    for j in range(n_atoms): # any big value is ok
        if del_num >= n_vacancy:
            break
        #p = random.randint(1,n_arrange+1)  
        p = plist[j]
        if R == 0:
            #Note that, currently, the latter condition means nothing.
            if perfect[p,4] ==1 and \
               (( perfect[p,1] >= (indent+0.5)*LAT_CONST+cut_off       and perfect[p,1] <= (N_X_D+indent)*LAT_CONST-cut_off) \
               or \
               ( perfect[p,1]  >= (N_X_D+0.5+indent)*LAT_CONST+cut_off and perfect[p,1] <= (N_X+indent)*LAT_CONST-cut_off)):
                for k in range(len(klist)):
                    if perfect[klist[k],5] ==0 and Vac_Cu ==1 :          # if Vac [k,5] ==0
                        if neighbor(p,klist[k],R,cut_off) == 1:
                            break
                    if perfect[klist[k],4] ==10 and Vac_Cu ==2 :          # if Cu [k,4] ==10
                        if neighbor(p,klist[k],R,cut_off) == 1:
                            break
                    if k == len(klist)-1:
                        if Vac_Cu ==1:    
                            perfect[p,5] =0           # if Vac[p,5] ==0
                            del_num +=1
                            klist.append(p)
                            print(len(klist)-1) 
                            #print(del_num)
                        if Vac_Cu ==2:
                            perfect[p,4] =10           # if Cu [p,4] ==10
                            del_num +=1
                            #print(del_num)
                            klist.append(p)
                            print(len(klist)-1)      
        else:
            #SUpposing the radius is more than 2.8A. If the cutoff is set to 2R(2*radius), then particles are scattered somewhat well. Currently virtually unused.
            #delta = 0.05
            #delta = 0.5
            #if perfect[p,4] ==1 and perfect[p,1] >= (N_X_D*(0.5-delta)+indent)*LAT_CONST+R+cut_off and perfect[p,1] <= (N_X_D*(0.5+delta)+indent)*LAT_CONST-(R+cut_off):]
            if perfect[p,4] ==1 and \
               (( perfect[p,1]>= (indent+0.5)*LAT_CONST+cut_off+R       and perfect[p,1] <= (N_X_D+indent)*LAT_CONST-cut_off-R) \
               or \
               ( perfect[p,1] >= (N_X_D+0.5+indent)*LAT_CONST+cut_off+R and perfect[p,1] <= (N_X+indent)*LAT_CONST-cut_off-R) ):
                for k in range(len(klist)):
                    if perfect[klist[k],5] ==0 and Vac_Cu ==1 :          # if Vac [k,5] ==0
                        if neighbor(p,klist[k],R,2*R) == 1:
                            break
                    if perfect[klist[k],4] ==10 and Vac_Cu ==2 :          # if Cu [k,4] ==10
                        if neighbor(p,klist[k],R,2*R) == 1:
                            break
                    if k == len(klist)-1:
                        for v in range(n_arrange):
                            if neighbor(p,v,R,0) ==1:
                                if Vac_Cu ==1:
                                    perfect[v,5] = 0  # if Vac [v,5] ==0
                                    del_num +=1 
                                    klist.append(v)
                                    #print(klist[del_num-1])
                                    #print(del_num)
                                    print(len(klist)-1) 
                                if Vac_Cu ==2:        
                                    perfect[v,4] = 10  # if Cu [v,4] ==10
                                    del_num +=1 
                                    #print(del_num)
                                    klist.append(v)
                                    print(len(klist)-1)
    
    name(Vac_Cu+8)

#create Fe.lat
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
    lat =0
    for i in range(n_atoms):
        if perfect[i,5] ==1 :
            lat +=1
            f.write("%12d %6d %12.6f %12.6f %12.6f\n" % (lat, perfect[i,4],perfect[i,1],perfect[i,2],perfect[i,3]))    
    f.close()

    with open("Fe.lat") as f:
        l = f.readlines()

    l.insert(2, "%d atoms" % lat)
    with open("Fe.lat", mode='w') as f:
        f.writelines(l)
    
    add(type_num)

#create add.lat, which tells you the initial position of buffer layer atoms.
def add(type_num):
    f = open("add.lat", "w")
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
    lat =0
    for i in range(n_atoms):
        if perfect[i,5] ==1 and (perfect[i,4] == 9 or perfect[i,4] == 7 or perfect[i,4] == 8):
            lat +=1
            f.write("%12d %6d %12.6f %12.6f %12.6f\n" % (lat, perfect[i,4],perfect[i,1],perfect[i,2],perfect[i,3]))
    f.close()

    with open("add.lat") as f:
        l = f.readlines()

    l.insert(2, "%d atoms" % lat)
    with open("add.lat", mode='w') as f:
        f.writelines(l)



#currently not used. Detail is described in the original file.
def input_empty():
    f = open("in.lammps", "a")
    for jx in range(indent):
        for jy in range(num):
            for jz in range(num):
                f.write("create_atoms %d single %12.6f %12.6f %12.6f\n" % (4, vert[0]+jx*LAT_CONST,vert[1]+jy*LAT_CONST,vert[2]+jz*LAT_CONST))
                f.write("create_atoms %d single %12.6f %12.6f %12.6f\n" % (4, surf1[0]+jx*LAT_CONST,surf1[1]+jy*LAT_CONST,surf1[2]+jz*LAT_CONST))

    for jx in range(N_X+indent, N_X+2*indent):
        for jy in range(num):
            for jz in range(num):
                f.write("create_atoms %d single %12.6f %12.6f %12.6f\n" % (4, vert[0]+jx*LAT_CONST,vert[1]+jy*LAT_CONST,vert[2]+jz*LAT_CONST))
                f.write("create_atoms %d single %12.6f %12.6f %12.6f\n" % (4, surf1[0]+jx*LAT_CONST,surf1[1]+jy*LAT_CONST,surf1[2]+jz*LAT_CONST))
    f.write("write_restart restart.equil\n")
    f.close()


#currently not used. Detail is described in the original file.
def inverse_lat():
    f = open("inverse.lat", "w")
    lat =0
    for i in range(perfect_num): 
        if perfect[i,5] ==0 : # if Cu [i,4] ==5
            lat +=1
            f.write("%12d %6d %12.6f %12.6f %12.6f\n" % (lat, perfect[i,4],perfect[i,1],perfect[i,2],perfect[i,3]))    
    f.close()
    

# This function determines whether the two atoms are closer than the cutoff length
# if two atoms are neighbors, it means they don't belong to different clusters; in other words they are very nearby)
def neighbor(p,k,R,cut):
     if distance(p,k,1)**2 + distance(p,k,2)**2 + distance(p,k,3)**2 <= (cut+R)**2:
         return 1
     else: 
         return 0

# This function culculates the distance in x / y / z direction
def distance(p,k,a):
    if a == 1:
        return abs(perfect[p,1]-perfect[k,1])
    elif a == 2:
        b1 = abs(perfect[p,2]-perfect[k,2])
        b2 = abs(perfect[p,2]-perfect[k,6])
        b3 = abs(perfect[p,6]-perfect[k,2])
        #theorically b4 is not needed
        return min(b1,b2,b3)
        
    else:
        b1 = abs(perfect[p,3]-perfect[k,3])
        b2 = abs(perfect[p,3]-perfect[k,7])
        b3 = abs(perfect[p,7]-perfect[k,3])
        #theorically b4 is not needed
        return min(b1,b2,b3)
            
            
# mainpart
if __name__ == "__main__":
    argvs = sys.argv
    make_perfect()
    arrange_lat(0,1)  #(R ,type) type= vacancy:1 , Cu: 2
    
    #input_empty()
                       #and RATE or num (=N_X,N_Y)are all parameters
                       #R unit is angstrome
                      

################moving test
#check num of arrage
#if you do not need, delete 'dump'
#
#Use OVITO to check molphorogy
#FFT configuration check
#arraving time should be defined consistnatly
#change x_detec, A2, Time to get beta
