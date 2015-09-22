# methyltransferase morph between 2zlb and 2zvj
# load start conformation and end conformation
load ../pdbs/2zlb.pdb, startConf
load ../pdbs/2zvj.pdb, finalConf

# the waters don't map 1-1 from structure to
# structure, so we remove them
remove solvent

# align the two structures
align startConf, finalConf, object=aln, cycles=0

# create the 2-state object morph_in,
# putting state 1 of startConf into
# state 1 of moph_in, and state 1 of
# finalConf into state 2 of morph_in, eg
# morph_in[state1] = startConf[state1]
# morph_in[state2] = finalConf[state1]
create morph_in, startConf in aln, 1, 1
create morph_in, finalConf in aln, 1, 2

# import the rigimol module
from epymol import rigimol

# execute RigiMOL on the morph_in structure and put the
# output into a new object called "morph_out"
rigimol.morph( "morph_in", "morph_out", refinement=2, async=1 )
