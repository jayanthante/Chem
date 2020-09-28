#!/bin/sh
# This program XYZifier converts NWChem-Venus output into a coordinate file which can
# be viewed as movie in Molden.
# Done by Dr. Jayanth K Ajay
# ./xyzifier.sh input_file.out framerate output_file.xyz
# $1 = input_file.out
# $2 = framerate
# $3 = output_file.xyz
frate=$2
if [ -z $1 ]; then
	echo "Filename not given"
	echo "Please run xyzifier as ./xyzifier.sh <input> <n=every nth coordinates to be taken> <output>"
	exit 1
else 
	if [ ! -f $1 ];	then 
		echo "$1 does not exist!!"
		exit 2
	fi
fi
if [ -z $2 ]; then
	echo "Framerate not given"
	echo "Please run xyzifier as ./xyzifier.sh <input> <n=every nth coordinates to be taken> <output>"
	exit 3
fi
if [ -z $3 ]; then
	echo "Output filename not given"
	newname=$(basename $1 .out).xyz
	read -r -p "Want to print the output to $newname ? (y/n)" yesno
	case "$yesno" in
	[Yy]) 
		outname=$newname;;
	[Nn])  
		echo "Quitting ...";
		exit 4;;
 	*) 
 		echo "Quitting ..";
 		exit 5;;
    esac
else
	newname=$3
fi
# reading the number of atoms from the .out file
natoms=$(grep "NUMBER OF ATOMS=" $1 |awk '{print $4}') 

# reading the atomic symbols as an array from the .out file
symbols=($(grep "labels" $1 |awk '{print $8}'))

# getting the number of steps from the .out file
nsteps=$(grep "XXXX" $1| wc -l)

# getting the space coordinates from the outfile and printing them after printing the number of atoms
# stores the formatted output to lala.tmp
# grep '       Q' -A $natoms pro.out | sed "s/ \+Q \+P/$natoms\n/" | sed '/--/d'| awk '{print "    "$1"    "$2"    "$3}' | sed 's/ -/-/g' > lala.tmp
grep '       Q' -A $natoms $1 | sed "s/ \+Q \+P/$natoms\n/" | awk '{print "    "$1"    "$2"    "$3}' | sed 's/ -/-/g' > lala.tmp

# the number of lines needed for one timestep in lala.tmp
# so a frame at a time contains N_atoms + 3 lines
framesize=$(($natoms + 3))

# Atoms are listed from 3rd line
# this is for inserting the atomic symbols in front of corresponding coordinates
atomframe=3
# for ii in `seq 0 $(($natoms-1))`; do
#	sed -i "$atomframe~$framesize s/^/${symbols[$ii]}/" lala.tmp
#	atomframe=$(($atomframe+1))
# 	echo "symbols  $ii ${symbols[$ii]}"
# done

# This inserts the atomic symbol for each of the given coordinates for all the timesteps
for atom in ${symbols[@]}; do
#	echo $atom;
 	sed -i "$atomframe~$framesize s/^/$atom/" lala.tmp #inserts atomic symbols
	atomframe=$(($atomframe+1)) # chooses the lines corresponding to a particular atom in the symbols array
done

# Removes the existing lili.tmp
rm -rf lili.tmp

# This is for printing every nth frame only;eg. only every 5th frame is printed
for j in `seq 0 $nsteps`; do
	if [ $(($j*$frate*$framesize)) -le $(( ($nsteps-1)*$framesize )) ]; then
		startline=$(($j*$frate*$framesize + 1)) # the position where the nth frame starts
		endline=$(($j*$frate*$framesize + $framesize)) # the position where the nth frame ends
		sed -n "$startline, $endline p" lala.tmp >> lili.tmp #copies the nt
	fi
done

# final touches
sed -i "s/ \+$natoms/$natoms/g" lili.tmp
sed -i '/--/d' lili.tmp 

mv lili.tmp $newname
rm -rf lala.tmp
