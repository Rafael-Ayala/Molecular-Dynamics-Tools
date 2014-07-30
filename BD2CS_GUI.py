
import numpy as np
import copy as copy
import sys
import os

from Tkinter import *
from tkFileDialog import *

input1 = None
input2 = None
input3 = None
input4 = None
input5 = None
input6 = None


def button1Click(event):
    global input1
    file = askopenfilename()
    input1 = file

def button2Click(event):
    global input2
    file = askopenfilename()
    input2 = file

def button3Click(event):
    global input3
    file = askopenfilename()
    input3 = file

def button4Click(event):
    global input4
    file = askopenfilename()
    input4 = file

def button6Click(event):
	global input6
	dir = askdirectory()
	input6 = dir

def buttonOKClick(event):
    global input5
    userCutOff = input5_entry.get()
    input5 = userCutOff

def buttonRUNClick(event):
    root.destroy()

root = Tk()
root.title("BD2CS.py")

Container1=Frame(root)
Container1.pack()

button1 = Button(Container1)
button1.configure(text="Select mobile protein PDB")
button1.pack(side=TOP)
button1.bind("<Button-1>", button1Click)

button2 = Button(Container1)
button2.configure(text="Select fixed protein PDB")
button2.pack(side=TOP)
button2.bind("<Button-1>", button2Click)

button3 = Button(Container1)
button3.configure(text="Select complexes file")
button3.pack(side=TOP)
button3.bind("<Button-1>", button3Click)

button4 = Button(Container1)
button4.configure(text="Select PDB with centers of masses")
button4.pack(side=TOP)
button4.bind("<Button-1>", button4Click)

button6 = Button(Container1)
button6.configure(text="Select output folder")
button6.pack(side=TOP)
button6.bind("<Button-1>", button6Click)

input5_label=Label(Container1)
input5_label["text"]="Enter distance cut off and click OK"
input5_label.pack()

input5_entry = Entry(Container1)
input5_entry["text"]="Select"
input5_entry.pack()

buttonOK = Button(Container1)
buttonOK.configure(text="OK")
buttonOK.pack()
buttonOK.bind("<Button-1>", buttonOKClick)

buttonRUN = Button(Container1)
buttonRUN.configure(text="Run")
buttonRUN.pack()
buttonRUN.bind("<Button-1>", buttonRUNClick)

root.mainloop()

print 1

os.chdir(input6) # Set output folder

# 1st argument: mobile protein PDB file
# 2nd argument: fixed protein PDB file
# 3rd argument: file with the centers of masses and X and Y orientation vectors of the mobile protein in all complexes
# 4th argument: PDB files with the coordinates of all the centers of masses referred to the same (0, 0, 0)
#               point than the PDB file of the fixed protein
# 5th argument: maximum distance between the amide N atoms of 2 residues to consider they make a contact
# 6th argument: route of the desired output folder

# read_pdb:
# Input parameters:
#   - file: string with the route of a PDB file
# Reads a PDB file and returns a list with one element for each atom. Each element contains:
#   - 1st element: atom type
#   - 2nd element: residue type
#   - 3rd element: residue number
#   - 4th element: coordinates of the atom, stored as a list of 3 elements
# All of them are stored as strings

def read_pdb(file):
 	global_list = []
 	for line in open(str(file)):
		list = line.split()
		id = list[0]
		if ((id == 'ATOM') or (id == 'HETATOM')):
            # check if there is a column to store the chain, to use the correct indices to extract each element
			if list[4].isalpha():
				atom_type = list[2]
				residue_type = list[3]
				residue_number = list[5]
				position = list[6:9]
				global_list.append([atom_type, residue_type, residue_number, position])
			else:
				atom_type = list[2]
				residue_type = list[3]
				residue_number = int(list[4])
				position = list[5:8]
				global_list.append([atom_type, residue_type, residue_number, position])
	return global_list

# center_of_masses:
# Input parameters:
#   - PDB_input: list returned by read_pdb
#   - gravimetric: if True, considers the masses of the atoms. If False, assumes all the atoms have the same mass
# Returns a list of 3 floats with the coordinates of the center of masses of the atoms stored in PDB_input

def center_of_masses(PDB_input, gravimetric=True):
    if gravimetric==True:
        masses=[]
        x_w=[]
        y_w=[]
        z_w=[]
        for i in range(len(PDB_input)):
            if PDB_input[i][0][0]=='H':
                masses.append(1)
                x_w.append(1*float(PDB_input[i][3][0]))
                y_w.append(1*float(PDB_input[i][3][1]))
                z_w.append(1*float(PDB_input[i][3][2]))
            elif PDB_input[i][0][0]=='C':
                masses.append(12)
                x_w.append(12*float(PDB_input[i][3][0]))
                y_w.append(12*float(PDB_input[i][3][1]))
                z_w.append(12*float(PDB_input[i][3][2]))
            elif PDB_input[i][0][0]=='N':
                masses.append(14)
                x_w.append(14*float(PDB_input[i][3][0]))
                y_w.append(14*float(PDB_input[i][3][1]))
                z_w.append(14*float(PDB_input[i][3][2]))
            elif PDB_input[i][0][0]=='O':
                masses.append(16)
                x_w.append(16*float(PDB_input[i][3][0]))
                y_w.append(16*float(PDB_input[i][3][1]))
                z_w.append(16*float(PDB_input[i][3][2]))
            elif PDB_input[i][0][0]=='S':
                masses.append(32)
                x_w.append(32*float(PDB_input[i][3][0]))
                y_w.append(32*float(PDB_input[i][3][1]))
                z_w.append(32*float(PDB_input[i][3][2]))
            elif PDB_input[i][0][0]=='F':
                masses.append(56)
                x_w.append(56*float(PDB_input[i][3][0]))
                y_w.append(56*float(PDB_input[i][3][1]))
                z_w.append(56*float(PDB_input[i][3][2]))
        return [sum(x_w)/sum(masses), sum(y_w)/sum(masses), sum(z_w)/sum(masses)]
    else:
        x_center_of_masses=sum([float(set_of_coords[0]) for set_of_coords in atomic_positions])/len(atomic_positions)
        y_center_of_masses=sum([float(set_of_coords[1]) for set_of_coords in atomic_positions])/len(atomic_positions)
        z_center_of_masses=sum([float(set_of_coords[2]) for set_of_coords in atomic_positions])/len(atomic_positions)
        return [x_center_of_masses, y_center_of_masses, z_center_of_masses]

# transform_coords:
# Input parameters:
#   - original_coords: list of the same format as the list returned by read_PDB
#   - complex_center_of_mass: center of mass of the mobile protein in one of the complexes
#   - cc_atom_transformation_matrix: numpy array storing a 3x3 matrix with the X, Y and Z rotation vectors of the mobile protein in one
#                                    of the complexes. Each row must be a vector.
# Returns a list of the same format as the list returned by read_PDB that represents the coordinates of each atom of the
# mobile protein in one of the complexes. A translation and a transformation are applied to the original coordinates of the mobile
# protein in order to calculate the new coordinates.

def transform_coords(original_coords, complex_center_of_mass, cc_atom_transformation_matrix):
	new_coords = copy.copy(original_coords)
	original_center_of_masses=center_of_masses(original_coords)
    
    # loop to go over all atoms
    
	for i in range(len(original_coords)):
        
        # extract the initial coordinates of an atom
        
		cc_atom_x_initial_coord=float(original_coords[i][3][0])
		cc_atom_y_initial_coord=float(original_coords[i][3][1])
		cc_atom_z_initial_coord=float(original_coords[i][3][2])
        
        # calculate the translation needed to move the center of mass of the protein to its position in the complex and apply it to
        # the current atom
        
		cc_atom_translation_coords=[cc_atom_x_initial_coord+(complex_center_of_mass[0]-original_center_of_masses[0]),cc_atom_y_initial_coord+(complex_center_of_mass[1]-original_center_of_masses[1]),cc_atom_z_initial_coord+(complex_center_of_mass[2]-original_center_of_masses[2])]
        
        # apply the transformation needed to rotate the protein to its orientation in the complex to the current atom
        
		cc_atom_transformation_coords=np.inner(cc_atom_translation_coords, np.transpose(cc_atom_transformation_matrix))
        
        # create an element to be added to the returned list corresponding to the current atom
        
		new_element = [new_coords[i][0:3], cc_atom_transformation_coords]
        
        # add the element to the returned list
        
		new_coords[i]=new_element
	return new_coords

# get_all_complexes_atomic_coords:
# Input parameters:
#   - original_coords: list of the same format as the list returned by read_PDB with the coordinates of the mobile protein
#   - complexes_center_of_mass: numpy array with 3 columens storing the coordinates of the center of mass of the mobile protein
#                               in all complexes
#   - complexes_x_vector: numpy array with 3 columns storing the coordinates of the X rotation vector of the mobile protein in
#                         all complexes
#   - complexes_y_vector: numpy array with 3 columns storing the coordinates of the Y rotation vector of the mobile protein in
#                         all complexes
#   - complexes_z_vector: numpy array with 3 columns storing the coordinates of the Z rotation vector of the mobile protein in
#                         all complexes
# Calculates the transformation matrix to be passed to transform_coords and returns a list containing a number of elements equal to the
# number of complexes. Each element is a list of the same format as the list returned by read_PDB, and stores the coordinates of all the
# atoms of the mobile protein in one of the complexes.

def get_all_complexes_atomic_coords(original_coords, complexes_center_of_mass, complexes_x_vector, complexes_y_vector, complexes_z_vector):
	lists_coords_every_complex=[]
	for j in range(len(complexes_x_vector)):
        
        # defines the transformation matrix
        
		cc_atom_transformation_matrix=np.array([complexes_x_vector[j],complexes_y_vector[j], complexes_z_vector[j]])
        
        # extracts the center of masses of the jth complex
        
		complex_center_of_mass=complexes_center_of_mass[j]
        
        # applies the translation and transformation needed to calculate the coordinates of all the atoms of the mobile protein
        # in the jth complex and appends the result to the returned list
        
		lists_coords_every_complex.append(transform_coords(original_coords, complex_center_of_mass, cc_atom_transformation_matrix))
	return lists_coords_every_complex

# distance_atoms:
# Input parameters:
#   - coordinates_1_input: list with the 3 coordinates of one point
#   - coordinates_2_input: list with the 3 coordinates of another point
# Converts all the coordinates to floats and calculates the distance between both points

def distance_atoms(coordinates_1_input, coordinates_2_input):
	coordinates_1=[float(x) for x in coordinates_1_input]
	coordinates_2=[float(x) for x in coordinates_2_input]
	distance_value=np.sqrt(np.power((coordinates_2[0]-coordinates_1[0]),2)+np.power((coordinates_2[1]-coordinates_1[1]),2)+np.power((coordinates_2[2]-coordinates_1[2]),2))
	return distance_value

# contact_residues:
# Input parameters:
#   - fixed_protein_coords: list with the coordinates of all the atoms of the fixed protein, with the same format as the list
#                           returned by read_PDB
#   - original_coords: list of the same format as the list returned by read_PDB with the coordinates of the mobile protein
#   - all_complexes_mobile_coords: list returned by get_all_complexes_atomic_coords
#   - cut_off: maximum distance between the amide N atoms of 2 residues to consider they make a contact
# Calculates the frequency with which each residue of the mobile protein appears in the interface of the complex and returns a list where
# each element is a list of 2 elements (the 1st element is the residue number and the 2nd element is the number of complexes where it
# is in the interface

def contact_residues(fixed_protein_coords, original_coords, all_complexes_mobile_coords, cut_off):
	residues_contact_count=[]
    
    # loop to go over each residue and calculate the number of complexes where it is in the interface
    
	for i in range(int(original_coords[len(original_coords)-1][2])):
        
        # calculates the indices of the elements of original_coords that store data of N atoms of the ith residue
        
		print "Residue"+str(i)+ "analyzed"
        
		N_atom_in_residue_indeces=[idx for idx in range(len(original_coords)) if (original_coords[idx][2]==(i+1) and original_coords[idx][0][0]=='N')]
        
        # initializes the count of complexes where the ith residue is in the interface
        
		contact_count=0
        
        # loop to go over each complex for the ith residue
        
		for j in range(len(all_complexes_mobile_coords)): # inicializa un bucle para contar el porcentaje de aparicion de un residuo
			contact_logical=None
            
            # extracts the coordinates of the amide N atom of the ith residue in the jth complex
            # [N_atom_in_residue_indeces[0]] is used instead of simply [N_atom_in_residue_indeces] because some residues contain more
            # than 1 N atom. In this case, the first element of [N_atom_in_residue_indeces] always corresponds to the amide N atom
            
			atom_mobile_protein_coords=all_complexes_mobile_coords[j][N_atom_in_residue_indeces[0]][1]
            
            # loop to calculate the distance between the amide N atom of the ith residue in the jth complex and all the N atoms
            # in the fixed protein
            
			for l in range(len(fixed_protein_coords)):
                
                # calculate the distance between amide N atoms only
                
				if fixed_protein_coords[l][0]=='N':
                    
                    # extracts the coordinates of the amide N atom of the fixed protein whose distance to the amide N atom of the
                    # ith residue in the jth complex is being calculated
                    
					atom_fixed_protein_coords=fixed_protein_coords[l][3]
                    
                    # calculate the distance
                    
					distance_between_atoms=distance_atoms(atom_mobile_protein_coords, atom_fixed_protein_coords)
                    
                    # if the distance is smaller than the set cut_off, add 1 to the contact count of the ith residue and stop the
                    # comparisons for the jth complex, since it has already been proved the ith residue is in the interface in the
                    # jth comple. Then start the calculations for the ith residue in the (j+1)th complex, or if the jth complex was
                    # the last one, for the (i+1)th residue in the 1st complex
                    
					if distance_between_atoms<cut_off:
						print distance_between_atoms
						contact_count+=1
						contact_logical=True
						break
		residues_contact_count.append([i+1, contact_count])
	return residues_contact_count

# load the PDB file of the mobile protein

initial_coords=read_pdb(input1)

# calculates the center of masses of the initial coordinates of the mobile protein

original_center_of_masses_mobile_protein=center_of_masses(initial_coords)

centers_of_masses_PDB_full_array=read_pdb(input4)

# extract only the spatial coordinates of the centers of masses of all the complexes. The coordinates are stored as strings

centers_of_masses_PDB_as_strings=[mass_center_line[3] for mass_center_line in centers_of_masses_PDB_full_array]

# convert the previous coordinates to floats

centers_of_masses_PDB=[]

for mass_center_line in centers_of_masses_PDB_as_strings:
    centers_of_masses_PDB.append([float(mass_center_coordinate) for mass_center_coordinate in mass_center_line])

# load the file with the centers of masses (with a different reference than the fixed protein PDB!!) and the rotation vectors of
# all the complexes

complexes_full_array=np.loadtxt(input3, comments="#", unpack=False, dtype='float')

# defines an array with 3 columns to store the X rotation vector coordinates of all the complexes

complexes_x_vector = complexes_full_array[:,5:8] # matriz con las coordenadas de los vectores de orientacion X

# defines an array with 3 columns to store the Y rotation vector coordinates of all the complexes

complexes_y_vector = complexes_full_array[:,8:11] # matriz con las coordenadas de los vectores de orientacion Y

# calculates the coordinates of the Z rotation vector for all the complexes and stores them in an array with 3 columns

complexes_z_vector = [0 for i in range(len(complexes_x_vector))]
for i in range(len(complexes_x_vector)):
	a=np.cross(complexes_x_vector[i], complexes_y_vector[i])
	complexes_z_vector[i]=a

# load the PDB file of the fixed protein

fixed_protein_coords=read_pdb(input2)

# calculates the atomic coordinates of every atom in all the complexes

allcomplexes=get_all_complexes_atomic_coords(initial_coords, centers_of_masses_PDB, complexes_x_vector, complexes_y_vector, complexes_z_vector)

input3=float(input5)

# determine the frequency with which each residue appears in the interface and stores it in a list

contactresidues=contact_residues(fixed_protein_coords, initial_coords, allcomplexes, input3)

# writen the frequency of each residue in a format suitable to define an attribute in Chimera

output_file=open('output.txt','a')

output_file.write('attribute: frequenceInComplex\n')
output_file.write('match mode: 1-to-1\n')
output_file.write('recipient: residues\n')

for x in range(len(contactresidues)):
    output_file.write('\t:%s\t%s\n' % (x+1, contactresidues[x][1]))


