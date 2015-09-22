#
# residue_properties.py -- a PyMOL script to map residue-level properties
#    onto the alpha carbons of a protein for a color-based representation
#    of the data.
#

# read in surface area values from disk

f = open("1eaz_b.txt", 'rb').readlines()

# set each atom's b-factor value to the one read from disk

alter *, b = f.pop(0)

# color by b-factor

preset.b_factor_putty("all")



python
f = open("1eaz_b.txt", 'wb')
for x in d:
    f.write("%3.2f\n" % x)
f.close()
python end

