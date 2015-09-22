#
# custom_view.py -- A PyMOL script to create a custom view for
#     any protein with a ligand.
#

# create a named selection for the ligand
select ligand, org

# show the ligand as sticks
as sticks, ligand

# orient the scene on the complete first small organic molecule
orient bm. first ligand

# create a named selection for all atoms with 5 Ang. of the ligand selection
select pocket, poly within 5 of ligand

# show the pocket selection only as surface
as surface, pocket

# find polar contacts from the pocket to the ligand
distance polar_contacts, poly, ligand, mode=2

# set all surfaces to draw as light grey
set surface_color, grey70

# make surfaces transparent for polymer atoms
set transparency, 0.5, poly

# turn on valences for the ligand
set_bond valence, 1, ligand

# position the labels 3.5 Ang. closer to the camera to un-obscure them
set label_position, [0, 0, 3.5]

# label the ligand with its name
label first ligand, "Molecule: %s" % (resn)
