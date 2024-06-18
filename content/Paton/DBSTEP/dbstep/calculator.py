# -*- coding: UTF-8 -*-
import math
import numpy as np
from scipy.spatial.transform import Rotation as R
import sys
from dbstep.constants import metals


"""
calculator

Performs calculations for finding angles, translation and rotation of molecules
"""


def unit_vector(vector):
	""" Returns the unit vector of the vector """
	return vector / np.linalg.norm(vector)
	

def point_vec(coords, spec_atom_2):
	"""returns coordinate vector between any number of atoms """
	point = np.array([0.0,0.0,0.0])
	
	for atom in spec_atom_2:
		point += coords[atom-1]
		
	return point


def rotate_mol(
		coords, spec_atom_1, end_point,
		verbose=False, atom3=False, cube_origin=False):
	"""Aligns spec_atom_1-end_point bond to the z-axis via rotation about the
	x and y axes.

	Args:
		coords (np.ndarray): 3D coordinates of the all the molecule's atoms
		spec_atom_1 (int): non-zero based index of the atom we're looking down
		end_point (np.ndarray): 3D coordinate of the atom we're looking to
		verbose (bool): should information about rotation be printed or not
		atom3 (int, optional): atom to specify rotation around z axis to align atom3 to the positive x direction & y=0
		cube_origin (np.ndarray, optional): origin a cube file

	Returns:
		Rotated version of coords and cube_origin (if provided).
	"""

	atom_1_index = spec_atom_1 - 1
	intersecting_vector = end_point - coords[atom_1_index]

	new_coords = np.copy(coords)
	new_cube_origin = np.copy(cube_origin)

	yaw, pitch, roll = 0, 0, 0

	yaw = angle_between_axis(intersecting_vector, 1, 2)
	if yaw != 0:
		print_rotation_info('x', yaw, verbose)
		end_point = apply_rotation(end_point, np.array([yaw, 0, 0]))
		intersecting_vector = end_point - coords[atom_1_index]

	pitch = angle_between_axis(intersecting_vector, 0, 2)
	if pitch != 0:
		print_rotation_info('y', pitch, verbose)
		end_point = apply_rotation(end_point, np.array([0, pitch, 0]))

	if atom3 is not False:
		atom3_coords = coords[int(atom3) - 1]
		intersecting_vector = atom3_coords - end_point
		roll = angle_between_axis(intersecting_vector, 1, 0)
		if roll != 0:
			roll = check_rotated_atom3_x_direction(atom3_coords, roll)
			print_rotation_info('z', roll, verbose)

	if yaw != 0 or pitch != 0 or roll != 0:
		# rotation of all of coords done here
		three_rotations = np.array([yaw, pitch, roll])
		new_coords = apply_rotation(new_coords, three_rotations)
		if cube_origin is not False:
			# effectively rotates whole grid that will be generated later
			new_cube_origin = apply_rotation(cube_origin, three_rotations)
	else:
		if verbose:
			print("   No rotation necessary :)")

	if cube_origin is False:
		return new_coords
	else:
		return new_coords, new_cube_origin


def angle_between_axis(vector, axis_from_index, axis_to_index):
	""" Returns the angle in radians needed to rotate axis_from to be parallel to axis_to. """
	v1_u = unit_vector(np.array(vector))

	angle = np.arctan2(v1_u[axis_from_index], v1_u[axis_to_index])

	# I've found through testing that these from-to combinations need their
	# rotations sign to be reversed.
	reverse_angle_combinations = [(0, 2), (1, 0), (2, 1)]
	if (axis_from_index, axis_to_index) in reverse_angle_combinations:
		angle = -angle

	return angle


def print_rotation_info(axis, radians, verbose):
	"""Prints rotation information if verbose and radians is non-zero."""
	if verbose:
		print(
			f'   Rotating molecule about {axis}-axis '
			f'{np.degrees(radians):.2f} degrees.')


def apply_rotation(item_to_rotate, radians):
	"""Rotates a vector or matrix about x, y and z axes specified by radians array."""
	rot = R.from_euler('xyz', radians)
	return np.round(rot.apply(item_to_rotate), 8)


def check_rotated_atom3_x_direction(atom3_coords, roll):
	"""If x is in negative direction, subtract or add pi to roll and return the result."""
	new_atom3_coords = apply_rotation(atom3_coords, np.array([0, 0, roll]))

	if new_atom3_coords[0] < 0:
		plus_pi, minus_pi = roll + np.pi, roll - np.pi
		roll = plus_pi if abs(plus_pi) < abs(minus_pi) else minus_pi
	return roll


def translate_mol(MOL, options, origin):
	"""# Translates molecule to place center atom at cartesian origin [0,0,0]"""
	coords, atoms, spec_atom = MOL.CARTESIANS, MOL.ATOMTYPES, options.spec_atom_1
	base_id = spec_atom - 1
	base_atom = atoms[base_id]
	try:
		displacement = coords[base_id] - origin
		if np.linalg.norm(displacement) == 0:
			if options.verbose: print("\n   Molecule is defined with {}{} at the origin".format(base_atom,(base_id+1)))
		else:
			if options.verbose == True: print("\n   Translating molecule by {} to set {}{} at the origin".format(-displacement, base_atom, (base_id+1)))
		for n, coord in enumerate(coords):
			coords[n] = coords[n] - displacement
	except:
		   sys.exit("   WARNING! Unable to find an atom to set at the origin")
	return coords


def translate_dens(mol, options, xmin, xmax, ymin, ymax, zmin, zmax, xyz_max, origin):
	""" Translates molecule so that a specified atom (spec_atom) is at the origin. Defaults to a metal if no atom is specified."""
	coords, atoms, cube_origin = mol.CARTESIANS, mol.ATOMTYPES,mol.ORIGIN
	spec_atom = options.spec_atom_1
	for n, atom in enumerate(atoms):
		if not spec_atom:
			if atom in metals:
				base_id, base_atom = n, atom
		else:
			if n+1 == spec_atom:
				base_id, base_atom = n, atom
	try:
		displacement = coords[base_id] - origin
		if np.linalg.norm(displacement) == 0:
			if options.verbose: print("\n   Molecule is already defined with {}{} at the origin".format(base_atom,(base_id+1)))
		else:
			if options.verbose: print("\n   Translating molecule by {} to set {}{} at the origin".format(-displacement, base_atom, (base_id+1)))
		for n, coord in enumerate(coords):
			coords[n] = coords[n] - displacement
		cube_origin = cube_origin + displacement
		xmin -= displacement[0]
		xmax -= displacement[0]
		ymin -= displacement[1]
		ymax -= displacement[1]
		zmin -= displacement[2]
		zmax -= displacement[2]
		xyz_max = max(xmax, ymax, zmax, abs(xmin), abs(ymin), abs(zmin))
	except:
		   sys.exit("   WARNING! Unable to find an atom (e.g. metal) to set at the origin")
	return [coords, cube_origin, xmin, xmax, ymin, ymax, zmin, zmax, xyz_max]