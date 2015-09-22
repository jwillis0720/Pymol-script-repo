'''
See more here: http://www.pymolwiki.org/index.php/com

DESCRIPTION
        get the center of mass of selection or move selection to the origin.
 
ARGUMENTS
        selection = string: a valid PyMOL selection {default: all}
        center = 0 or 1: if center=1 center the selection {default: 0}
        returns: center of mass: [ xCOG, yCOG, zCOG ]
 
SEE ALSO
        get_extent, get_position, http://pymolwiki.org/index.php/Center_Of_Mass

'''
from pymol import cmd
from pymol import stored
from chempy import cpv

def COG(selection='all', center=0, quiet=1):

        model = cmd.get_model(selection)
        nAtom = len(model.atom)
 
        COG = cpv.get_null()
 
        for a in model.atom:
                COG = cpv.add(COG, a.coord)
        COG = cpv.scale(COG, 1./nAtom)
 
        if not int(quiet):
                print ' COG: [%8.3f,%8.3f,%8.3f]' % tuple(COG)
 
        if int(center):
                cmd.alter_state(1, selection, "(x,y,z)=sub((x,y,z), COG)",
                        space={'COG': COG, 'sub': cpv.sub})
 
        return COG
 
cmd.extend("COG", COG)
