from pymol import cmd

def merge(selection='*', output="merged_output"):
	for i in cmd.get_object_list(selection=selection):
		cmd.create(output,i,1,-1)

cmd.extend('merge',merge)
