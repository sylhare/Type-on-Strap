import pytest

import numpy as np
from numpy.random import default_rng
from scipy.spatial.transform import Rotation as R

from dbstep import calculator


@pytest.fixture
def matrix_to_rotate(rng):
	"""A matrix with shape (2,3) filled with random ints ranging from 0-100"""
	return rng.integers(0, 100, size=(2, 3))


@pytest.fixture(params=[1] * 6)
def rotation_radians(request, rng):
	"""An np.npdarry of rotations about the x, y then z axis in radians. Parametrized to give 6 different outputs."""
	# 1/3 chance that rotation about a certain axis is 0
	return rng.uniform(0.1, 2*np.pi, size=3) * rng.integers(0, 2, size=3) * request.param


@pytest.fixture
def rng():
	"""numpy's random number generator."""
	return default_rng()


def test_apply_rotation(matrix_to_rotate, rotation_radians):
	"""Cross validates overall euler rotation with multiple matrix rotations gotten from SciPy's from_rotvec."""
	rot_matrix_x = np.array(R.from_rotvec([rotation_radians[0], 0, 0]).as_matrix())
	rot_matrix_y = np.array(R.from_rotvec([0, rotation_radians[1], 0]).as_matrix())
	rot_matrix_z = np.array(R.from_rotvec([0, 0, rotation_radians[2]]).as_matrix())

	expected = rot_matrix_z @ rot_matrix_y @ rot_matrix_x @ matrix_to_rotate.T
	actual = calculator.apply_rotation(matrix_to_rotate, rotation_radians)

	assert np.isclose(expected.T, actual).all()


quadrant_points_2d = [[0.5, 1.], [0.5, -1.], [-0.5, 1.], [-0.5, -1]]

@pytest.fixture(
	params=quadrant_points_2d,
	ids=[str(point) for point in quadrant_points_2d])
def quadrant_point(request):
	"""parametrized with four quadrants of a 2D circle."""
	return request.param


axis_combinations = [
	('z', 'x', 'y'), ('x', 'z', 'y'),
	('y', 'x', 'z'), ('x', 'y', 'z'),
	('z', 'y', 'x'), ('y', 'z', 'x')]

@pytest.fixture(
	params=axis_combinations,
	ids=[str(comb) for comb in axis_combinations])
def axis_combination(request):
	"""parametrized with (axis_from, axis_to, zero_axis).

	NOTE:
		Characters are used in axis combinations to improve human readability.
	"""
	axis_indices = {'x': 0, 'y': 1, 'z': 2}
	comb = request.param
	return axis_indices[comb[0]], axis_indices[comb[1]], axis_indices[comb[2]]


@pytest.fixture
def vector_to_rotate(quadrant_point, axis_combination):
	"""Gives full 3D vector using quadrant_point and axis_combination."""

	axis_from_index, axis_to_index, zero_axis_index = axis_combination

	full_vector = quadrant_point.copy()
	full_vector.insert(zero_axis_index, 0)
	return full_vector


def test_angle_between_axis(axis_combination, vector_to_rotate):
	axis_from_index, axis_to_index, zero_axis_index = axis_combination

	angle = calculator.angle_between_axis(vector_to_rotate, axis_from_index, axis_to_index)
	rotation_vector = [0, 0, 0]
	rotation_vector[zero_axis_index] = angle
	rotated_vector = calculator.apply_rotation(vector_to_rotate, rotation_vector)

	error_message = (
		f'failed to rotate {vector_to_rotate} using from: {axis_from_index}, to: {axis_to_index}\n')
	assert rotated_vector[axis_from_index] == 0, error_message


def test_rotate_mol(vector_to_rotate):
	coords = np.array([[0,0,0], vector_to_rotate])
	spec_atom_1 = 1
	end_point = vector_to_rotate
	atom3 = False

	new_coords = calculator.rotate_mol(coords, spec_atom_1, end_point)

	assert new_coords[1][0] == 0
	assert new_coords[1][1] == 0
