README

Available at https://github.com/cfagnen/Bilayer_depth_script
Author: Charline Fagnen
These scripts compute the depth of the two leaflets of a bilayer.
The script command_fatslim.sh is based on the Fatslim tool developed by S. Buchoux. (https://github.com/FATSLiM/fatslim)
The script depth_script.py is inspired by the one written by J. Chong (https://github.com/jiehanchong/membrane-depth-analysis) 


The process for determining the depth of the membrane is divided into two scripts.
-command_fatslim.sh : a shell script detecting which leaflet the lipids belong to.
-depth_script.py : a python3 script computing the depths of the upper and lower leaflets.

Required:
---------

- FATSLiM (Fast Analysis Toolbox for Simulations of Lipid Membranes)
  (https://pythonhosted.org/fatslim/index.html)
  S. Buchoux. FATSLiM: a fast and robust software to analyze MD simulations of membranes. Bioinformatics. (2017), doi:10.1093/bioinformatics/btw563 . 
- MDAnalysis
  (https://www.mdanalysis.org/)
   R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L. Seyler, D. L. Dotson, J. Domanski, S. Buchoux, I. M. Kenney, and O. Beckstein. MDAnalysis: A Python package for the rapid analysis of molecular dynamics simulations. In S. Benthall and S. Rostrup, editors, Proceedings of the 15th Python in Science Conference, pages 98-105, Austin, TX, 2016. SciPy, doi:10.25080/majora-629e541a-00e.
   N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein. MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations. J. Comput. Chem. 32 (2011), 2319-2327, doi:10.1002/jcc.21787. PMCID:PMC3144279
- Gromacs
  H. Bekker, H.J.C. Berendsen, E.J. Dijkstra, S. Achterop, R. van Drunen, D. van der Spoel, A. Sijbers, and H. Keegstra et al., “Gromacs: A parallel computer for molecular dynamics simulations”; pp. 252–256 in Physics computing 92. Edited by R.A. de Groot and J. Nadrchal. World Scientific, Singapore, 1993.


Before to begin:
----------------

The scripts need the frames composed of all the membrane atoms (only membrane) located in a directory called frames.
The scripts also need the definition of the membrane in a ndx file called all_membrane.ndx.

It can be done with gromacs:
#Create the directory with the frames containing only the membrane atoms.

gmx make_ndx -f trajectory.tpr -o all_membrane.ndx
x|y|z # where x,y,z, etc... are the indices of the differents lipids
name XX Membrane
q
mkdir frames
cd frames
gmx trjconv -fit progressive -f ../trajectory.xtc -s ../trajectory.tpr -n ../all_membrane.ndx -sep -o frame.gro
cd ..

#Note: sometimes, the files obtained by this command failed
In this case:
gmx trjconv -fit progressive -f ../trajectory.xtc -s ../trajectory.tpr -n ../all_membrane.ndx -sep -o frame.pdb
And the convert them in .gro:


#Create the ndx files describing only the membrane atoms. 
Example:

gmx make_ndx -f frames/frame0.gro -o membrane.ndx
  0 System              : 71965 atoms
  1 Other               : 71965 atoms
  2 POPC                : 40236 atoms
  3 POPE                : 15324 atoms
  4 POPS                :  1896 atoms
  5 CHOL                : 10216 atoms
  6 POP2                :  2544 atoms
  7 DPSM                :  1749 atoms

#Save the complete system and rename it headgroups, removing other groups.
> 0
Copied index group 0 'System'
  8 System              : 71965 atoms
> name 8 headgroups
> del 0-7
Removed group 0 'System'
Removed group 1 'Other'
Removed group 2 'POPC'
Removed group 3 'POPE'
Removed group 4 'POPS'
Removed group 5 'CHOL'
Removed group 6 'POP2'
Removed group 7 'DPSM'
> q



command_fatslim.sh
------------------


This script attributes a leaflet to each lipid.
If the original cutoff value does not help to distinguish the leaflets, the cutoff will be incremented by 0.1 until it reaches the cutoff_max.
If fatslim does not succeed in detecting the leaflets with cutoff_max, the frame will be moved to a new directory.

Parameters to set (lines 5-9 in the script):
nbfirstframe=0            # Number of the first frame
nblastframe=7500          # Number of the last frame
nbthreads_for_fatslim=6   # Number of the threads to use
cutoff_original=2         # Fatslim parameter to determine
cutoff_max=7              # Highest cutoff value allowed


Determination of cutoff_original
The cutoff is defined like it by fatslim's developers:
"For each lipid, the point cloud considered is consituted by all the beads within a user-tweakable cutoff distance (see --cutoff) from the reference lipid."

The cutoff_original should be 2 or 3. Sometimes 2 is too restrictive and 3 can save us time (by avoiding unsuccessful attempts).

It is possible to test fatslim with your system for one frame:
#Try with cutoff (cutoff_original) = 3
fatslim membranes -c frames/frame0.gro -n membrane.ndx -o fatslim/frame0.xvg --output-index fatslim/bilayer_leaflet_0.ndx --cutoff 3  --nthreads 6
#Creation of the .gro files for each leaflet to visualize the seperation in VMD
gmx  editconf -f frames/frame0.gro -n fatslim/bilayer_leaflet_0_0000.ndx -o leaflet1.gro (Select 0)
gmx  editconf -f frames/frame0.gro -n fatslim/bilayer_leaflet_0_0000.ndx -o leaflet2.gro (Select 1)


Launch (Approximatively: 2h05min for 5 replicas launched in parallel - 5 replicas * 7500 frames - 6 threads - box size: 44*44*24 nm^3):
./command_fatslim.sh

This script creates 2 directory: frames_unresolved_by_fatslim and fatslim.
The first contains the frame where fatslim fails. If fatslim is effective with all the frames, the directory is empty.
The directory fatslim contains the results of fatslim. The depth_script.py will use the bilayer_leaflet_{frame_number}_0000.ndx files that indicate the indices of the atoms belonging to each leaflet.

Three files are created to check that fatslim works well:
- nb_of_assigned_atoms_per_frame.dat (First column = frame number, Second column = number of atoms clustered in one leaflet by fatslim)
- percent_assigned_atoms_per_frame.dat (First column = frame number, Second column = percentage of atoms clustered in one leaflet by fatslim)
- nb_of_clusters_detected_per_frame.dat (First column = frame number, Second column = number of clusters detected by fatslim)



depth_script.py (Python3)
-------------------------

This script:
- assigns a leaflet to the unassigned atoms (the same than the closest atom assigned to a leaflet with fatslim),
- computes the depth of the upper and the lower leaflets.
 
The user needs to give one variable (line 7):
	z_orientation = 1 # put 1 if the system is oriented in the same direction as z axis
	z_orientation = -1 # put -1 if the system is upside-down
and check if all the lipids of the membrane studied are indicated in the lipid dictionary (line 8), if it is not the case, add it/them:
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

Launch: (Approximatively: 3h15min for 5 replicas launched in parallel - 5 replicas * 7500 frames - 6 threads - box size: 44*44*24 nm^3):

python depth_script.py

The script creates one directory called figures, where for each frame, one figure is computed showing the definition of the lower (red) and upper (blue) leaflets according to x and z coordinates. 

You can create a movie from these pictures with the following command (installation of ffmpeg required):
cd figures
ffmpeg -f image2 -i graph%d.png -r 24 -vcodec mpeg4 -b 15000k ../leaflets_definition.mp4
cd ..

The depth_script.py creates 19 files: 
depths.dat : timestamp in the first column, upper leaflet depth in the second column, and lower leaflet depth in the third column.
distri_leaflets.dat : Frame number in the first column, the second and third columns: percentage of lipids in upper and lower leaflets, respectively.
depths.txt : average depth of upper and lower leaflets.
upper_leaflet_depth_over_time.svg : Depth of the upper leaflet versus time.
lower_leaflet_depth_over_time.svg : Depth of the lower leaflet versus time.
upper_vs_lower.svg : Figure describing the percentage of lipids in each leaflet.
control_number_lipids.svg : Number of lipids used in the depth computation by frame.
upper_heatmap.dat, lower_heatmap.dat : matrix of the depths at each point for the upper and lower leaflets, respectively.
rescounts_upper_heatmap.dat, rescounts_lower_heatmap.dat = Value of bins for upper and lower maps, respectively.


