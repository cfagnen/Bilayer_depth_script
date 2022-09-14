#!/usr/bin/python3

#Available at https://github.com/cfagnen/Bilayer_depth_script
#Author: Charline Fagnen

# indicate the orientation of your system with z_orientation
# lipid_dict : add the missing lipids you use in your simulations

#VARIABLES
########################
z_orientation = 1 # 1 : if the system is oriented in the same direction as z axis, -1 if it is upside-down
lipid_dict = {
        "POPC ": 1,
        "POPE ": 2,
        "POPS ": 3,
        "POP2 ": 4,
        "DPSM ": 5,
        "CHOL ": 6,
        "PIPC ": 7,
        "PUPE ": 8,
        "PAPS ": 9,
        "PAPC ": 10
        }
########################

if z_orientation == 1:
    z_orientation = 1
    print("OK, your system is oriented in the same direction than z axis.")
elif z_orientation == -1:
    z_orientation = -1
    print("OK, your system is upside down.")
else:
    print("Error: z_orientation parameter must be 1 (same direction than z axis) or -1 (upside down).")
    sys.exit()

import numpy as np
import sys, os
import matplotlib
matplotlib.use('agg') #allows plotting without use of X-server, to allow use over remote connection
import matplotlib.pyplot as plt
np.set_printoptions(suppress=True)
import MDAnalysis as mda
import re
import warnings

#Allow to ignore the warning about the missing mass of atoms
warnings.filterwarnings("ignore", message="Failed to guess the mass for the following atom types:")


#Generate list of all the .gro files in the directory
gro_list = np.array([f for f in os.listdir('frames') if f.endswith('.gro')])
number_of_files = len(gro_list)

bins = 75 #number of bins in heatmap

depths = np.zeros(3) #initialise table of depths
distri_leaflets = np.zeros(4) #initialise table of distri

os.mkdir('figures')

def initialize_heatmap():
    heatmap = np.zeros((bins, bins)) #initialise heatmap
    heatmap[:] = np.nan #initialise heatmap
    hm_residues_per_bin = np.zeros((bins, bins)) #initialise count of residues per bin for heatmap averaging
    return heatmap, hm_residues_per_bin #package heatmap and residues data into a tuple for feeding into function

upper_heatmap_data = initialize_heatmap()
lower_heatmap_data = initialize_heatmap()
frame_n = 1 #initialise frame counter


def import_fatslim_results():
    atoms = u.atoms
    leaflet1_name="membrane_1_leaflet_1"
    leaflet2_name="membrane_1_leaflet_2"
    nbline = 0
    totnbline = 0
    with open("fatslim/" + leafletsepfilename,"r") as file:
        for line in file:
            nbline = nbline + 1
            totnbline = totnbline +1
            if re.search(leaflet2_name, line):
                nblim = nbline

    #Select all lines except the last one (because dimensions could be different)
    line_begin_leaflet2 = nblim + 1
    line_end_leaflet2 = totnbline - 1 
    line_begin_leaflet1 = 2
    line_end_leaflet1 = nblim - 2 

    #Last line
    line_true_end_leaflet1 = nblim - 2
    line_true_end_leaflet2 = totnbline - 1 
    
    #Extraction of the data
    leaflet1temp = []
    leaflet1term = []
    leaflet2temp = []
    leaflet2term = []
    i = 0
    with open("fatslim/" + leafletsepfilename, "r+") as f:
        for line in f:
            i = i + 1
            if i in range(line_begin_leaflet1, line_end_leaflet1):
                leaflet1temp.append(line.split())
            if i in range(line_begin_leaflet2, line_end_leaflet2):
                leaflet2temp.append(line.split())
            if i == int(line_true_end_leaflet1):
                leaflet1term.append(line.split())
            if i == int(line_true_end_leaflet2):
                leaflet2term.append(line.split())
    
    #Leaflet1: Change the shape of the data (because here one line is one tuple but it is false)
    leaflet1_dim1_old = np.shape(leaflet1temp)[0]
    leaflet1_dim2_old = np.shape(leaflet1temp)[1]
    leaflet1_nbtot_elements = leaflet1_dim1_old * leaflet1_dim2_old
    leaflet1_atom_list_temp = np.reshape(leaflet1temp, (leaflet1_nbtot_elements,))

    #Leaflet2: Change the shape of the data (because here one line is one tuple but it is false)
    leaflet2_dim1_old = np.shape(leaflet2temp)[0]
    leaflet2_dim2_old = np.shape(leaflet2temp)[1]
    leaflet2_nbtot_elements = leaflet2_dim1_old * leaflet2_dim2_old
    leaflet2_atom_list_temp = np.reshape(leaflet2temp, (leaflet2_nbtot_elements,))

    #Leaflet1 Last line: Change the shape of the data (because here one line is one tuple but it is false)
    leaflet1term_dim1_old = np.shape(leaflet1term)[0]
    leaflet1term_dim2_old = np.shape(leaflet1term)[1]
    leaflet1term_nbtot_elements = leaflet1term_dim1_old * leaflet1term_dim2_old
    leaflet1_atom_list_term = np.reshape(leaflet1term, (leaflet1term_nbtot_elements,))

    #Leaflet2 Last line: Change the shape of the data (because here one line is one tuple but it is false)
    leaflet2term_dim1_old = np.shape(leaflet2term)[0]
    leaflet2term_dim2_old = np.shape(leaflet2term)[1]
    leaflet2term_nbtot_elements = leaflet2term_dim1_old * leaflet2term_dim2_old
    leaflet2_atom_list_term = np.reshape(leaflet2term, (leaflet2term_nbtot_elements,))

    #Concatenation 
    leaflet1_atom_list = np.concatenate((leaflet1_atom_list_temp, leaflet1_atom_list_term), axis=0)
    leaflet2_atom_list = np.concatenate((leaflet2_atom_list_temp, leaflet2_atom_list_term), axis=0)
    leaflet1_atom_list_int = leaflet1_atom_list.astype(int)
    leaflet2_atom_list_int = leaflet2_atom_list.astype(int)

    #Add of -1 because the count begins to 0, so the first atom is the atom with index 0
    leaflet1_good_indices = np.add(leaflet1_atom_list_int, -1)
    leaflet2_good_indices = np.add(leaflet2_atom_list_int, -1)
    
    #Selection of the atoms ROH and PO4 atoms with MDAnalysis 
    leaflet1fatslim = atoms[leaflet1_good_indices].select_atoms("name PO4 or name ROH")
    leaflet2fatslim = atoms[leaflet2_good_indices].select_atoms("name PO4 or name ROH")

    #Atoms not distributed by Fatslim
    all_po4_roh_atoms = u.select_atoms(f'name PO4 or name ROH')
    atoms_po4_roh_excluded = all_po4_roh_atoms - leaflet1fatslim - leaflet2fatslim
    coords_leaflet1fatslim = leaflet1fatslim.positions
    coords_leaflet2fatslim = leaflet2fatslim.positions
    coords_atoms_po4_roh_excluded  = atoms_po4_roh_excluded.positions

    #distribute the excluded atom in the same group that the closest atom
    line_to_check = 0
    leaflet1 = coords_leaflet1fatslim.astype(int)
    leaflet2 = coords_leaflet2fatslim.astype(int)
    while line_to_check < len(coords_atoms_po4_roh_excluded):
        reference = coords_atoms_po4_roh_excluded[line_to_check][[0, 1, 2]]
        reference_new_shape = np.reshape(reference, (1,3))
        reference_to_add = reference_new_shape.astype(int)

        difference_from_leaflet1 = coords_leaflet1fatslim - reference
        distance_from_leaflet1 = (difference_from_leaflet1[:, 0]**2 + difference_from_leaflet1[:, 1]**2 + difference_from_leaflet1[:, 2]**2)**0.5

        difference_from_leaflet2 = coords_leaflet2fatslim - reference
        distance_from_leaflet2 = (difference_from_leaflet2[:, 0]**2 + difference_from_leaflet2[:, 1]**2 + difference_from_leaflet2[:, 2]**2)**0.5
 
        mindist_from_leaflet1 = distance_from_leaflet1.min()
        mindist_from_leaflet2 = distance_from_leaflet2.min()

        if mindist_from_leaflet2 >= mindist_from_leaflet1: # so the atom belongs to leaflet1
                leaflet1 = np.concatenate((leaflet1, reference_to_add))
        else:   # so the atom belongs to leaflet2
                leaflet2 = np.concatenate((leaflet2, reference_to_add))
        line_to_check += 1

    nb_atoms_leaflet1 = len(leaflet1)
    nb_atoms_leaflet2 = len(leaflet2)
    nb_atoms_considered = nb_atoms_leaflet1 + nb_atoms_leaflet2

    return leaflet1, leaflet2, nb_atoms_leaflet1, nb_atoms_leaflet2, nb_atoms_considered


def figure_distribution(coordleaflet1, coordleaflet2):
    if  z_orientation == 1:
        plt.plot(coordleaflet1[:, 0], coordleaflet1[:, 2], "o", color='red', markeredgecolor="k", markersize=4)
        plt.plot(coordleaflet2[:, 0], coordleaflet2[:, 2], "o", color='blue', markeredgecolor="k", markersize=4)
        plt.title("Frame " + numberframe)
        plt.xlabel('x')
        plt.ylabel('z')
        plt.savefig(f"figures/graph{numberframe}.png")
        plt.clf()
    else:
        plt.plot(coordleaflet1[:, 0], coordleaflet1[:, 2], "o", color='blue', markeredgecolor="k", markersize=4)
        plt.plot(coordleaflet2[:, 0], coordleaflet2[:, 2], "o", color='red', markeredgecolor="k", markersize=4)
        plt.title("Frame " + numberframe)
        plt.xlabel('x')
        plt.ylabel('z')
        plt.savefig(f"figures/graph{numberframe}.png")
        plt.clf()


def get_depth(leaflet): #identifies the top and bottom 10% of residues by height, and returns the difference of their means
        top_cutoff = np.percentile(leaflet[:, 2], 90)
        top_ten_percent = leaflet[:, 2][leaflet[:, 2] > top_cutoff]
        average_top = np.mean(top_ten_percent)
        bottom = np.min(leaflet[:, 2])
        depth = average_top - bottom
        return depth


for path in gro_list:

    #First we parse the resid and coordinates from .gro file
    frame = open("frames/" + path)
    numberframe = str(''.join(list(filter(str.isdigit, path))))
    print("Analysis in progress of the frame:", numberframe)
    leafletsepfilebeg = "bilayer_leaflet_"
    leafletsepfileend = "_0000.ndx"
    leafletsepfilename = leafletsepfilebeg + numberframe + leafletsepfileend


    data = frame.read().splitlines() #read the frame to a list, one line per list
    numlen = 0
    coords = np.zeros((0,5))
    firstrow = data[0] # get timestamp
    time_starts_at = str.find(firstrow, ' t=') + 3 #find the position of the string 't=', which precedes the timestamp
    if time_starts_at > 10:     #this will be true if there is a timestamp
        timestamp = float(firstrow[time_starts_at:])
    else:   #if there is no timestamp, then set timestamp to 0
        timestamp = 0

    frame.close()


    u = mda.Universe("frames/" + path)
    
    leaflet1coord, leaflet2coord, nb_atoms_leaflet1, nb_atoms_leaflet2, nb_atoms_considered = import_fatslim_results()

    leaflet1_mean_height = np.mean(leaflet1coord[:, 2])
    leaflet2_mean_height = np.mean(leaflet2coord[:, 2])
   

    #Define which leaflet is the upper and which one is the lower one
    if (leaflet1_mean_height >  leaflet2_mean_height) and ( z_orientation == 1 ): #Assign leaflets to upper or lower
        upper_leaflet = leaflet1coord
        lower_leaflet = leaflet2coord
        percent_atoms_upper = nb_atoms_leaflet1 * 100 / nb_atoms_considered
        percent_atoms_lower = nb_atoms_leaflet2 * 100 / nb_atoms_considered
    elif (leaflet2_mean_height > leaflet1_mean_height) and ( z_orientation == 1 ):
        upper_leaflet = leaflet2coord
        lower_leaflet = leaflet1coord
        percent_atoms_upper = nb_atoms_leaflet2 * 100 / nb_atoms_considered
        percent_atoms_lower = nb_atoms_leaflet1 * 100 / nb_atoms_considered
    elif (leaflet1_mean_height >  leaflet2_mean_height) and ( z_orientation == -1 ):
        upper_leaflet = leaflet2coord
        lower_leaflet = leaflet1coord
        #*-1 again to have the upper and lower membrane in the good direction in the png files
        upper_leaflet[:,2] *= -1
        lower_leaflet[:,2] *= -1
        percent_atoms_upper = nb_atoms_leaflet2 * 100 / nb_atoms_considered
        percent_atoms_lower = nb_atoms_leaflet1 * 100 / nb_atoms_considered
    elif (leaflet2_mean_height >  leaflet1_mean_height) and ( z_orientation == -1 ):
        upper_leaflet = leaflet1coord
        lower_leaflet = leaflet2coord
        upper_leaflet[:,2] *= -1
        lower_leaflet[:,2] *= -1
        percent_atoms_upper = nb_atoms_leaflet1 * 100 / nb_atoms_considered
        percent_atoms_lower = nb_atoms_leaflet2 * 100 / nb_atoms_considered
    elif leaflet1_mean_height == leaflet2_mean_height:
        print('Error: Leaflet heights are the same')
        sys.exit()


    numberframefl = float(numberframe)
    distri_leaflets = np.vstack((distri_leaflets, [numberframefl ,percent_atoms_upper, percent_atoms_lower, nb_atoms_considered]))
    figure_distribution(leaflet1coord,leaflet2coord)

    #*-1 again to have absolute value for heatmap
    if z_orientation == -1:
        upper_leaflet[:,2] *= -1
        lower_leaflet[:,2] *= -1
 

    figures_list = np.array([f for f in os.listdir('figures') if f.endswith('.png')])
    number_of_figures = len(figures_list)
    percentage_advanced = number_of_figures * 100 / number_of_files
    print(number_of_figures, " Frames / ",number_of_files, " Frames : ", percentage_advanced, "%")
    #Now obtain the depths
    upper_depth = get_depth(upper_leaflet)
    lower_depth = get_depth(lower_leaflet)

    depths = np.vstack((depths, [timestamp, upper_depth, lower_depth]))
    #Generate heatmap of z-coordinates over time
    lastline_split = np.array(data[-1].split(' ')) #extract box size from last line
    lastline_numbers = lastline_split[lastline_split != '']
    range_x = float(lastline_numbers[0]) * 10 # convert nm in angstroms
    range_y = float(lastline_numbers[1]) * 10

    def make_heatmap(coords, heatmap_data, bins, range_x, range_y): #adds to a pre-initialised heatmap
        heatmap = heatmap_data[0]
        hm_residues_per_bin = heatmap_data[1]

        x_binwidth = range_x * 1.5 / bins
        y_binwidth = range_y * 1.5 / bins
        x_start = -0.25 * range_x #heatmap origin, which is outside the actual simulation in order to capture all information
        y_start = -0.25 * range_y 

        for ybin in range(1, bins+1):                                               #in each row
            ybin_start = coords[:, 1] > (y_start + y_binwidth * (ybin - 1))   #boolean for y > start of bin
            ybin_end = coords[:, 1] <= (y_start + y_binwidth * ybin)              #boolean for y < end of bin
            ybin_selector = ybin_start * ybin_end
            this_ybin = coords[ybin_selector, :]                        #select residues in this ybin


            for xbin in range(1, bins+1):                                           #go through the bins
                xbin_start = this_ybin[:, 0] > (x_start + x_binwidth * (xbin - 1))
                xbin_end = this_ybin[:, 0] <= (x_start + x_binwidth * xbin)
                xbin_selector = xbin_start * xbin_end
                thisbin = this_ybin[xbin_selector, :]

                if np.isnan(heatmap[xbin - 1, ybin -1]): #if the heatmap bin is currently empty

                    total_z_this_bin = sum(thisbin[:, 2])
                    hm_residues_per_bin[xbin - 1, ybin -1] = len(thisbin[:, 2])
                else: #if the heatmap bin already contains data
                    total_z_this_bin = (heatmap[xbin - 1, ybin -1] * hm_residues_per_bin[xbin - 1, ybin -1]) + sum(thisbin[:, 2]) #generate a new total from the previous contents and the current frame
                    hm_residues_per_bin[xbin - 1, ybin -1] = hm_residues_per_bin[xbin - 1, ybin -1] + len(thisbin[:, 2]) #increase the residue count accordingly

                if len(thisbin[:, 2]) > 0: #if there was anything in this bin in this frame
                    heatmap[xbin - 1, ybin -1] = total_z_this_bin / hm_residues_per_bin[xbin - 1, ybin -1] #update the average value in the heatmap

        return (heatmap, hm_residues_per_bin)

    upper_heatmap_data = make_heatmap(upper_leaflet, upper_heatmap_data, bins, range_x, range_y)
    lower_heatmap_data = make_heatmap(lower_leaflet, lower_heatmap_data, bins, range_x, range_y)

    frame_n = frame_n + 1

depths = np.delete(depths, 0, axis=0) #remove initialisation column
depths = depths[depths[:, 0].argsort()] #sort depths by timestamp

distri_leaflets = np.delete(distri_leaflets, 0, axis=0) #remove initialisation column
distri_leaflets = distri_leaflets[distri_leaflets[:, 0].argsort()] #sort depths by timestamp

#Export depths and distribution of the leaflets to files
np.savetxt("depths.dat", depths)
np.savetxt("distri_leaflets.dat", distri_leaflets) #, fmt="%s")

#Write mean depths to file
f = open('depths.txt', 'w+')
f.write("Depth of upper leaflet dome = " + str(np.mean(depths[:, 1])) + " A\n" +
        "Depth of lower leaflet dome = " + str(np.mean(depths[:, 2])) + " A")
f.close()

#Output graph of change in depth over course of simulation
u_DoT = plt.plot(depths[:, 0], depths[:, 1], linewidth=0.5)
plt.savefig("upper_leaflet_depth_over_time.svg")
plt.close()

l_DoT = plt.plot(depths[:, 0], depths[:, 2], linewidth=0.5)
plt.savefig("lower_leaflet_depth_over_time.svg")
plt.close()

#Control of the number of atoms take into account
u_DoT = plt.plot(distri_leaflets[:, 0], distri_leaflets[:, 1], linewidth=0.5, color='blue', label='Upper Leaflet')
u_DoT = plt.plot(distri_leaflets[:, 0], distri_leaflets[:, 2], linewidth=0.5, color='red', label ='Lower Leaflet')
plt.title("Distribution of lipids in the leaflets")
plt.xlabel('Frames')
plt.ylabel('Percentage(%)')
plt.legend()
plt.savefig("upper_vs_lower.svg")
plt.close()

u_DoT = plt.plot(distri_leaflets[:, 0], distri_leaflets[:, 3], linewidth=0.5, color='purple')
plt.title("Number of lipids taken into account")
plt.xlabel('Frames')
plt.ylabel('Number of lipids')
plt.savefig("control_number_lipids.svg")
plt.close()

#Output heatmap data
np.savetxt("upper_heatmap.dat", upper_heatmap_data[0])
np.savetxt("rescounts_upper_heatmap.dat", upper_heatmap_data[1])

np.savetxt("lower_heatmap.dat", lower_heatmap_data[0])
np.savetxt("rescounts_lower_heatmap.dat", lower_heatmap_data[1])

hmplot = plt.imshow(upper_heatmap_data[0])
plt.colorbar(hmplot)
plt.savefig("upper_heatmap.pdf")
plt.savefig("upper_heatmap.svg")
plt.close()

hmplot = plt.imshow(lower_heatmap_data[0])
plt.colorbar(hmplot)
plt.savefig("lower_heatmap.pdf")
plt.savefig("lower_heatmap.svg")
plt.close()

ctplot = plt.contourf(upper_heatmap_data[0])
plt.colorbar(ctplot)
plt.savefig("upper_contour.pdf")
plt.savefig("upper_contour.svg")
plt.close()

ctplot = plt.contourf(lower_heatmap_data[0])
plt.colorbar(ctplot)
plt.savefig("lower_contour.pdf")
plt.savefig("lower_contour.svg")
plt.close()

