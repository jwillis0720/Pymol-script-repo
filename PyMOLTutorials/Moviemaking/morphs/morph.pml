#
# generic_morph.pml -- a skeleton PyMOL script for the generic morphing of two objects
#

# README
#
# This script loads two files from the ../pdbs directory and 
# creates a morph from the first protein to the second.  It 
# then writes the morph to disk as a multi-state PDB file and
# quits PyMOL.
#
# Please run generic_morph_to_movie.pml after this script to 
# create the biologically relevant movie.
#
# This script writes the file
#      morph_out.pdb
# a multi-state PDB file to the working directory.
#

# CODE
#
# Load the 1st object
#
load file1.pdb, startConf

#
# Load the 2nd object
#
load file2.pdb, finalConf

#
# remove the solvent
#
remove solvent

#
# Align the two objects, saving the alignment object, aln.
#
align startConf, finalConf, object=aln

#
# Create the input object for RigiMOL from the startConf
# and finalConf objects.  startConf goes into state 1 of
# the new object and finalConf goes into state 2.
#
create morph_in, startConf, 1, 1
create morph_in, finalConf, 1, 2

#
# now we import the morph command and run it
#
from epymol import rigimol

#
# execute RigiMOL on the morph_in structure and put the
# output into a new object called "morph_out"
#
rigimol.morph( "morph_in", "morph_out", refinement=10, async=1 )

#
# uncomment the following two lines to stop the movie and save the file
#
# mstop
# save morph_out.pdb, morph_out, state=0

