# -*- coding: UTF-8 -*-
import os
from dbstep.constants import BOHR_TO_ANG

"""
writer

Contains classes and methods to write data to files

"""


class Logger:
	"""Enables output to terminal and to text file"""
	# Designated initializer 
	def __init__(self,filein,suffix,append):
		# Create the log file at the input path
		self.log = open(filein+"_"+append+"."+suffix, 'w' )
	# Write a message only to the log and not to the terminal
	def Writeonlyfile(self, message):
		self.log.write(message+"\n")


class WriteCubeData:
	""" Write new cube file of translated, rotated molecule for PyMOL """
	def __init__(self, file, cube):
		self.FORMAT = 'cube'
		oldfile = open(file+"."+self.FORMAT,"r")
		oldlines = oldfile.readlines()
		molfile = open(file+"_radius."+self.FORMAT,"w")
		# write new coordinates to file
		for line in oldlines[0:2]:
			molfile.write(line)

		dims = [cube.xdim,cube.ydim,cube.zdim]
		molfile.write("{:5} {:11.6f} {:11.6f} {:11.6f} {:4}".format(len(cube.ATOMNUM),cube.ORIGIN[0] / BOHR_TO_ANG, cube.ORIGIN[1] / BOHR_TO_ANG, cube.ORIGIN[2] / BOHR_TO_ANG, 1)+'\n')
		for i in range(len(cube.INCREMENTS)):
			molfile.write("{:5} {:11.6f} {:11.6f} {:11.6f}".format(dims[i],cube.INCREMENTS[i][0] / BOHR_TO_ANG,cube.INCREMENTS[i][1] / BOHR_TO_ANG,cube.INCREMENTS[i][2] / BOHR_TO_ANG)+'\n')
		for i in range(len(cube.CARTESIANS)):
			x = cube.CARTESIANS[i][0] / BOHR_TO_ANG
			y = cube.CARTESIANS[i][1] / BOHR_TO_ANG
			z = cube.CARTESIANS[i][2] / BOHR_TO_ANG
			molfile.write("{:5} {:11.6f} {:11.6f} {:11.6f} {:11.6f}".format(cube.ATOMNUM[i],float(cube.ATOMNUM[i]),x,y,z)+'\n')


		for line in cube.DENSITY_LINE:
			molfile.write(line)

		# there may well be a fast way to do this directly from the 3D array of X,Y,Z points
		# see http://paulbourke.net/dataformats/cube/
		# for x in xvals:
		# 	for y in yvals:
		# 		width = []
		# 		for z in zvals:
		# 			width.append((x**2 + y**2)**0.5)
		# 		# for cube print formatting
		# 		list_of_widths = itertools.zip_longest(*(iter(width),) * 6)
		# 		for widths in list_of_widths:
		# 			outline = ''
		# 			for val in widths:
		# 				if val != None:
		# 					outline += str('  {:10.5E}'.format(val))
		# 			molfile.write('\n'+outline)
		molfile.close()


def pymol_export(file, mol, spheres, cylinders, isoval, visv, viss):
	"""Outputs a python script that can be imported into PyMol (with 'run script.py')"""
	base, ext = os.path.splitext(file)

	log = Logger(base, "py", "steric")
	log.Writeonlyfile('from pymol.cgo import *')
	log.Writeonlyfile('from pymol import cmd\n')
	
	if visv == 'circle' or viss: 
		log.Writeonlyfile('def cgoCircle(r=3.5, cr=1.0, cg=0.4, cb=0.8, w=8.0, z=0.0, name=None):'
		'\n  """'
		'\n  Create a CGO circle'
		'\n  from https://pymolwiki.org/index.php/CgoCircle'
		'\n  PARAMS'
		'\n        x, y, z'
		'\n          X, Y and Z coordinates of the origin'
		'\n        r'
		'\n          Radius of the circle'
		'\n        cr, cg, cb'
		'\n          Color triplet, [r,g,b] where r,g,b are all [0.0,1.0].'
		'\n        w'
		'\n          Line width of the circle'
		'\n  RETURNS'
		'\n        the CGO object (it also loads it into PyMOL, too).'
		'\n'
		'\n  """'
		'\n  x, y = 0.0, 0.0'
		'\n' 
		'\n  r = abs(float(r))'
		'\n  cr = abs(float(cr))'
		'\n  cg = abs(float(cg))'
		'\n  cb = abs(float(cb))'
		'\n  w = float(w)'
		'\n'
		'\n  obj = [ BEGIN, LINES, COLOR, cr, cg, cb ]'
		'\n  for i in range(180):'
		'\n        obj.append( VERTEX )'
		'\n        obj.append(r*math.cos(i) + x )'
		'\n        obj.append(r*math.sin(i) + y )'
		'\n        obj.append(z)'
		'\n        obj.append( VERTEX )'
		'\n        obj.append(r*math.cos(i+0.1) + x )'
		'\n        obj.append(r*math.sin(i+0.1) + y )'
		'\n        obj.append(z)'
		'\n  obj.append(END)'
		'\n  if cb == 0.0:' 
		'\n        if cr == 0.0:'
		'\n              cName = "BminC_"+str(z)'
		'\n        else:'
		'\n              cName = "BmaxC_"+str(z)'
		'\n  else:'
		'\n        if name is not None:'
		'\n            cName = "circle_"+str(r)+name'
		'\n        else:'
		'\n            cName = "circle_"+str(r)'
		'\n  cmd.load_cgo( obj, cName )'
		'\n  cmd.set("cgo_line_width", w, cName )'
		'\n  return obj')
		log.Writeonlyfile('\ncmd.extend( "cgoCircle", cgoCircle )')
	
	#visualize volumes 
	if len(spheres) == 1:
		sphere = spheres[0]
		sphere_r = abs(float(sphere.split(',')[-2].strip()))
		if visv == 'circle':
			log.Writeonlyfile('cgoCircle(r='+str(sphere_r)+',name="a")')
			log.Writeonlyfile('cgoCircle(r='+str(sphere_r)+',name="b")')
			log.Writeonlyfile('cgoCircle(r='+str(sphere_r)+',name="c")')
			log.Writeonlyfile('cmd.rotate([1,0,0],90,object="circle_'+str(sphere_r)+'b",origin=[0,0,0])')
			log.Writeonlyfile('cmd.rotate([1,0,0],90,object="circle_'+str(sphere_r)+'c",origin=[0,0,0])')
			log.Writeonlyfile('cmd.rotate([0,0,1],90,object="circle_'+str(sphere_r)+'c",origin=[0,0,0])')
		else:
			sphere_id = "sphere_"+str(sphere_r)
			log.Writeonlyfile('\nsphere = [')
			log.Writeonlyfile(sphere)
			log.Writeonlyfile(']')
			log.Writeonlyfile('cmd.load_cgo(sphere,"'+sphere_id+'")')
	else:
		for n,sphere in enumerate(spheres):
			sphere_r = abs(float(sphere.split(',')[-2].strip()))
			if visv == 'circle':
				log.Writeonlyfile('cgoCircle(r='+str(sphere_r)+')')
			else:
				sphere_id = "sphere_"+str(sphere_r)
				log.Writeonlyfile('\nsphere = [')
				log.Writeonlyfile(sphere)
				log.Writeonlyfile(']')
				log.Writeonlyfile('cmd.load_cgo(sphere,"'+sphere_id+'")')
		
	#visualize sterimol param cylinders (cylinders ordered: Bmax, Bmin then L is always last)
	log.Writeonlyfile('\ncylinder = [')
	for n,cyl in enumerate(cylinders): log.Writeonlyfile(cyl)
	log.Writeonlyfile(']\ncmd.load_cgo(cylinder, '+'"axes"'+')\n')
	
	#visualize circles outlining Bmin, Bmax - need z, r and rgb colors
	if viss:
		for i in range(len(cylinders)-1): 
			cyl_vals = cylinders[i].split(',')
			z = cyl_vals[3].strip()
			x = float(cyl_vals[4].strip())
			y = float(cyl_vals[5].strip())
			r = str((x**2+y**2)**0.5)
			cr = cyl_vals[11].strip()
			cg = cyl_vals[12].strip()
			log.Writeonlyfile('cgoCircle(r='+r+',z='+z+',cr='+cr+',cg='+cg+',cb=0.0)')

	full_path = os.path.abspath(file)
	name, ext = os.path.splitext(full_path)

	if ext == '.cube':
		log.Writeonlyfile('\ncmd.load("'+full_path+'", '+'"dens"'+')')
		log.Writeonlyfile('\ncmd.load("'+name+'_radius.cube", '+'"distances"'+')')
		log.Writeonlyfile('\ncmd.isosurface("isodens", "dens", '+str(isoval)+')')

	log.Writeonlyfile('\ncmd.load("'+name+'_transform.xyz")')
	log.Writeonlyfile('cmd.show_as("spheres", "'+base.split('/')[-1]+'_transform")')
	log.Writeonlyfile('cmd.set("sphere_transparency", 0.5)')
	log.Writeonlyfile('cmd.set("orthoscopic", "on")')


def xyz_export(file,mol):
	"""Write xyz coordinates of molecule to file"""
	name, ext = os.path.splitext(file)
	log = Logger(name, "xyz", "transform")
	log.Writeonlyfile(str(len(mol.ATOMTYPES)))
	log.Writeonlyfile(name.split('/')[-1].split('\\')[-1])
	coords = ''
	for i in range(len(mol.ATOMTYPES)):
		coords += mol.ATOMTYPES[i]+'\t'
		for j in range(3):
			coords += "{0:.8f}".format(mol.CARTESIANS[i][j])+'\t'
		coords +='\n'
	log.Writeonlyfile(coords)