
from collections import OrderedDict
import cmd
def get_sane_pairing(pairing):
    l = {}
    for position in pairing:
        try:
            position_n = int(position[0])
        except ValueError:
            position_n = position[0]
        try:
            position_l = longer_names[position[1]]
        except KeyError:
            print 'Warning: {} not found'.format(position[1])
            continue
        l[position_n] = position_l
    return l

cmd.extend('get_sane_pairing',get_sane_pairing)

holder = []
def printer(resi,resn,name):
    holder.append((resn ,resi))

def get_resi_in_selection(objsel="(all)"):
    print objsel
    myspace = {'myfunc':printer}
    cmd.iterate('{}'.format(objsel),'myfunc(resi,resn,name)',space=myspace)
    for res_pair in list(OrderedDict.fromkeys(holder)):
        print "{},{}".format(res_pair[1],res_pair[0]),
        print "\t",


cmd.extend('get_resi_in_selection',get_resi_in_selection)