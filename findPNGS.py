from pymol import cmd as _cmd
from collections import OrderedDict

def find_pngs(selection):
	#print _cmd.iterate()
	chains = cmd.get_chains(selection)
	for chain in chains:
		#rint chain
		#seq_dict = OrderedDict()
		myspace = {'resis': OrderedDict()}
		cmd.iterate(selection + " and chain {}".format(chain), 'resis[resi] = resn', space=myspace)
		resn_resi =  myspace['resis']
		residues = resn_resi.keys()
		residue_names = resn_resi.values()
		pseudo_enumerate = 0
		for i,j in zip(residues, residue_names):
			if j == 'ASN' or j == "S56" or j == "TES":
				next_aa = residue_names[pseudo_enumerate + 1]
				third_aa = residue_names[pseudo_enumerate + 2]
				if next_aa != "PRO":
					if third_aa in ['THR','SER']:
						cmd.select("Glycan_{}{}".format(chain,i),'chain {} and resi {}'.format(chain,i))
			pseudo_enumerate += 1
	cmd.group("PNGS",'Glycan_*')
cmd.extend('find_PNGS',find_pngs)
