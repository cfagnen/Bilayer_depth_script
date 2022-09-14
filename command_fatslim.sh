#!/bin/bash

#Available at https://github.com/cfagnen/Bilayer_depth_script
#Author: Charline Fagnen

#VARIABLES
########################
nbfirstframe=0
nblastframe=1
nbthreads_for_fatslim=6
cutoff_original=3
cutoff_max=7
########################

mkdir fatslim frames_unresolved_by_fatslim

#number total of atoms in the membrane
nbatoms_total=$(cat frames/frame${nbfirstframe}.gro | head -2 | tail -1)
#Fatslim results accepted only if more than 95% of the atoms are considered (/1 is to truncate)
limit_nbatoms=$(echo "${nbatoms_total} * 0.95 / 1" | bc )

#FATSLIM
for i in `seq ${nbfirstframe} ${nblastframe}`
do
#Launch fatslim
fatslim membranes -c frames/frame${i}.gro -n membrane.ndx -o fatslim/frame${i}.xvg --output-index fatslim/bilayer_leaflet_${i}.ndx --cutoff ${cutoff_original}  --nthreads ${nbthreads_for_fatslim}
#Check the number of atoms clustered
nbatoms_clustered=$(grep -v "mem" fatslim/bilayer_leaflet_${i}_0000.ndx | tr ' ' '\n' | sed '/^$/d' | wc | awk '{print $1}')
nb_clusters=$(grep "mem" fatslim/bilayer_leaflet_${i}_0000.ndx | wc | awk '{print $1}')


#If the number of clustered atoms < 95% and if the number of clusters detected is not 2, relaunch fatslim with a cutoff incremented by 0.1
cutoff_update=$cutoff_original
while [ "${nbatoms_clustered}" -le "${limit_nbatoms}" ] || [ "${nb_clusters}" != 2 ]
do
   if [ "$(bc <<< "$cutoff_max > $cutoff_update")" == "1" ]
   then
       rm fatslim/bilayer_leaflet_${i}_0000.ndx
       rm fatslim/frame${i}.xvg
       cutoff_update=$(echo "${cutoff_update} + 0.1 " | bc )
       fatslim membranes -c frames/frame${i}.gro -n membrane.ndx -o fatslim/frame${i}.xvg --output-index fatslim/bilayer_leaflet_${i}.ndx --cutoff ${cutoff_update}  --nthreads ${nbthreads_for_fatslim}
       nbatoms_clustered=$(grep -v "mem" fatslim/bilayer_leaflet_${i}_0000.ndx | tr ' ' '\n' | sed '/^$/d' | wc | awk '{print $1}')
       nb_clusters=$(grep "mem" fatslim/bilayer_leaflet_${i}_0000.ndx | wc | awk '{print $1}')
   fi
done

#If cutoff max reached, move the .gro file of the frame in the directory frames_unresolved_by_fatslim with its results
if (( $(echo "$cutoff_max == $cutoff_update" | bc -l) ))
then
   mv fatslim/bilayer_leaflet_${i}_0000.ndx frames_unresolved_by_fatslim/
   mv fatslim/frame${i}.xvg frames_unresolved_by_fatslim/
   mv frames/frame${i}.gro frames_unresolved_by_fatslim/
fi

done

rm fatslim/frame*.xvg

#Final : number of atoms taken into account by frame
for i in `seq ${nbfirstframe} ${nblastframe}`
do
if [ -f fatslim/bilayer_leaflet_${i}_0000.ndx ]
then
nbatoms_clustered_final=$(grep -v "mem" fatslim/bilayer_leaflet_${i}_0000.ndx | tr ' ' '\n' | sed '/^$/d' | wc | awk '{print $1}')
echo "${i} ${nbatoms_clustered_final}"
fi
done > nb_of_assigned_atoms_per_frame.dat

#Final: percentage of atoms taken into account by frame
for i in `seq ${nbfirstframe} ${nblastframe}`
do
if [ -f fatslim/bilayer_leaflet_${i}_0000.ndx ]
then
nbatoms_clustered_final=$(grep -v "mem" fatslim/bilayer_leaflet_${i}_0000.ndx | tr ' ' '\n' | sed '/^$/d' | wc | awk '{print $1}')
nbatoms_clustered_final_percent=$(echo "$nbatoms_clustered_final $nbatoms_total" | awk '{printf "%.6f \n", $1*100/$2}')
echo "${i} ${nbatoms_clustered_final_percent}"
fi
done > percent_assigned_atoms_per_frame.dat

#Final: number of clusters (=leaflets)
for i in `seq ${nbfirstframe} ${nblastframe}`
do
if [ -f fatslim/bilayer_leaflet_${i}_0000.ndx ]
then
nb_clusters_final=$(grep "mem" fatslim/bilayer_leaflet_${i}_0000.ndx | wc | awk '{print $1}')
echo "${i} ${nb_clusters_final}"
fi
done > nb_of_clusters_detected_per_frame.dat

