# load start and final conformations
load ../pdbs/2zlb.pdb, startConf
load ../pdbs/2zvj.pdb, finalConf

# remove solvent atoms, they don’t map 1-1 across structures
remove solvent

# align the two structures
align startConf, finalConf, object=aln

# make 100-state morph_in from state 1 of startConf and state 1 of finalConf
for x in range(51): cmd.create("m_out", "startConf", 1, x)
for x in range(51,101): cmd.create("m_out", "finalConf", 1, x)

# perform the linear smoothing across all states
smooth m_out, 1, 100, 1, 100

# import the rigimol module
from epymol import rigimol

# compute the morph and store the resulting
# trajectory in a new object called "morph_out"
rigimol.refine(2, "m_out")
