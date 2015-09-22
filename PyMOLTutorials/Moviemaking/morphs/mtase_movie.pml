#
# generic_morph_to_movie.pml -- a skeleton script that uses a morph and holo structure
# 			     to make a more biologically interesting movie
#

# README
#
# This script loads the morph from generic_morph.pml.  It then loads
# the holo structure to obtain the ligands and any solvent nearby the
# binding pocket.  Lastly, this movie inserts and shows the ligands
# in the molecular morph movie.
#
# This script will output a multi-state PDB file and an MPEG file
# which contains the movie called,
# 	morph_movie.mpeg
#

# CODE
#
# Load the multi-state pdb file
load morph_out.pdb

#
# remove all solvent molecules (waters, etc)
#
remove solvent

#
# load the holo structure
#
load ../pdbs/2zvj.pdb, tmp

#
# create two objects, one called ligand and the other called waters.  These objects are 
# created FROM state 1 of the tmp object just loaded and will be inserted into state 30
# of the molecular morph movie.
#
create ligand, tmp and organic, 1, 30
create waters, tmp and (solvent within 5 of organic), 1, 30

#
# we are done with the tmp object, so we remove it
#
delete tmp

#
# we use the old-style mset command to make a movie that does a state-sweep.  This
# will start in state 1 for 30 frames, then transition up to state 30 from 1.  Next, it will
# stop at state 30 for 30 frames.  Then, it goes back to state one, taking another 30 frames.
#
# Here's a tip: the mset command takes two types of arguments eg, 
# 	 (1) 1x30, which means repeat state 1 for 30 frames, and
#	 (2) 1-30, which means transition from state 1 to 30
# Just write out your movie sequence as,
#      mset 1x30 1-30 30x30 30-1
# and then insert a space before any 'x' or a '-' to
# convert this into a valid PyMOL command.
#
mset 1 x30 1 -30 30 x30 30 -1

mplay
reset
orient ligand
show_as sticks, ligand
zoom center, 12
turn x, 90

cmd.movie.produce("morph_movie.mpeg", quiet=0)
