from pymol import cmd

def useRosettaRadii(): 
	cmd.alter("element C", "vdw=2.00") 
	cmd.alter("element N", "vdw=1.75") 
	cmd.alter("element O", "vdw=1.55") 
	cmd.alter("element H", "vdw=1.00") 
	cmd.alter("element P", "vdw=2.15") 
	cmd.alter("element S", "vdw=1.90") 
	cmd.alter("element RE", "vdw=1.40") 
	cmd.alter("element CU", "vdw=1.40") 
	cmd.set("sphere_scale", 1.0)

def expandRadii(delta=1.0, sel='all'): 
	for a in cmd.get_model(sel).atom:
		r = float(a.vdw) + float(delta) 
		cmd.alter("index "+a.index,'vdw='+r) 
		cmd.rebuild(sel,"spheres")

def contractRadii(delta=1.0, sel='all'): 
	for a in cmd.get_model(sel).atom:
		r = float(a.vdw) - float(delta) 
		cmd.alter("index "+a.index,'vdw='+r) 
		cmd.rebuild(sel,"spheres")

def useTempRadii(sel="all"): 
	for ii in range(30): 
		radius = "%0.1f"%(float(ii+1)/10) 
		cmd.alter(sel+" and b="+radius,"vdw="+radius) 
		cmd.rebuild()

def showpacking(sel="all"): 
	useRosettaRadii() 
	useTempRadii(sel+" and resn CAV") 
	cmd.hide('everything', sel+" and resn CAV" ) 
	cmd.show('spheres', sel+" and resn CAV" )

cmd.extend('useRosettaRadii', useRosettaRadii) 
cmd.extend('expandRadii',expandRadii) 
cmd.extend('contractRadii',contractRadii) 
cmd.extend("useTempRadii",useTempRadii) 
cmd.extend("showpacking",showpacking)
