import pytest

from dbstep import Dbstep
import numpy as np

class DbstepShell:
	"""A reference to the dbstep class so its helper methods can be accessed without running its __init__"""
	def __init__(self):
		self.__class__ = Dbstep.dbstep


class TestGetSpecAtoms:
	"""Tests _get_spec_atoms"""

	@pytest.mark.parametrize("input_atom_1, input_atom_2, expected_atom_1, expected_atom_2", [
		(None, None, 1, [2]),
		(1, "4,6,1,4", 1, [4, 6, 1, 4]),
		(3, [3, 4, "5"], 3, [3, 4, 5])
	])
	def test_get_spec_atoms(self, input_atom_1, input_atom_2, expected_atom_1, expected_atom_2):
		options = Dbstep.set_options({'atom1': input_atom_1, 'atom2': input_atom_2})
		db = DbstepShell()
		db._get_spec_atoms(options)
		assert options.spec_atom_1 == expected_atom_1
		assert options.spec_atom_2 == expected_atom_2

	@pytest.mark.parametrize("input_atom_1, input_atom_2, expected_exception", [
		(1, "aa", ValueError),
		("aa", [1], ValueError),
		(3, [1, 3, "4a"], ValueError),
		(3, [3, 4, [1, 2, 3]], TypeError),
	])
	def test_get_spec_atoms_exception(self, input_atom_1, input_atom_2, expected_exception):
		options = Dbstep.set_options({'atom1': input_atom_1, 'atom2': input_atom_2})
		db = DbstepShell()
		with pytest.raises(expected_exception):
			db._get_spec_atoms(options)


class SetupAgainstVerloop:
	"""Gets data ready to be put into parametrized methods in TestAgainstVerloop."""

	file_names = [
		'H.xyz', 'Me.xyz', 'Et.xyz', 'iPr.xyz', 'nBu.xyz', 'CH2iPr.xyz', 'cHex.xyz', 'nPr.xyz', 'Ad.xyz',
		'tBu.xyz', 'CH2tBu.xyz', 'CHEt2.xyz', 'CHiPr2.xyz', 'CHPr2.xyz', 'CEt3.xyz', 'Ph.xyz', 'Bn.xyz',
		'4ClPh.xyz', '4MePh.xyz', '4MeOPh.xyz', '35diMePh.xyz', '1Nap.xyz']

	verloop_ls = [
		2.15, 3.2, 4.2, 4.2, 6.26, 5.25, 6.26, 5.25, 6.26, 4.2, 5.25, 5.25, 5.25, 6.26, 5.25, 6.37, 4.62,
		7.69, 7.42, 8.29, 6.37, 6.37]
	verloop_b1s = [
		1.09, 1.7, 1.7, 1.99, 1.7, 1.7, 2., 1.7, 3.25, 2.86, 1.7, 1.99, 2.14, 1.99, 2.86, 1.71, 1.7,
		1.75, 1.71, 1.81, 1.71, 1.71]
	verloop_b5s = [
		1.09, 2.13, 3.26, 3.26, 4.63, 4.54, 3.58, 3.58, 3.58, 3.26, 4.54, 4.54, 4.54, 5.76, 4.54, 3.2, 6.11,
		3.2, 3.2, 3.2, 4.38, 5.59]

	def __init__(self):
		self.dbstep_ls, self.dbstep_b1s, self.dbstep_b5s = self._dbstep_outputs()

		self.file_names_and_ls = []
		self.file_names_and_b1s = []
		self.file_names_and_b5s = []

		self._fill_parameter_lists()

	def _dbstep_outputs(self):
		"""Reads and calculates dbstep parameters for all of the files in self.file_names"""
		dbstep_ls, dbstep_b1s, dbstep_b5s = [], [], []
		for file in self.file_names:
			db_obj = Dbstep.dbstep('dbstep/examples/' + file, sterimol=True, measure='classic',  commandline=True)
			dbstep_ls.append(np.round(db_obj.L + 0.4, 2))
			dbstep_b1s.append(np.round(db_obj.Bmin, 2))
			dbstep_b5s.append(np.round(db_obj.Bmax, 2))
		return dbstep_ls, dbstep_b1s, dbstep_b5s

	def _fill_parameter_lists(self):
		"""Fills the file name and Sterimol parameters lists so they can be passed into pytest.mark.parametrize."""
		values = zip(
			self.file_names,
			self.verloop_ls, self.verloop_b1s, self.verloop_b5s,
			self.dbstep_ls, self.dbstep_b1s, self.dbstep_b5s)
		for file_name, verloop_l, verloop_b1, verloop_b5, dbstep_l, dbstep_b1, dbstep_b5 in values:
			self.file_names_and_ls.append((file_name, verloop_l, dbstep_l))
			self.file_names_and_b1s.append((file_name, verloop_b1, dbstep_b1))
			self.file_names_and_b5s.append((file_name, verloop_b5, dbstep_b5))


class TestAgainstVerloop:
	"""Compares dbstep sterimol parameters to parameters calculated by Verloop's fortran code.

		NOTE: I decided to give each parameter their own tolerance variable so they can
		be changed in isolation.
	"""
	sav = SetupAgainstVerloop()

	@pytest.mark.parametrize("name, verloop_l, dbstep_l", sav.file_names_and_ls)
	def test_l_against_verloop(self, name, verloop_l, dbstep_l):
		tolerance = 0.01
		self.compare_with_tolerance(verloop_l, dbstep_l, tolerance)

	@pytest.mark.parametrize("name, verloop_b1, dbstep_b1", sav.file_names_and_b1s)
	def test_b1_against_verloop(self, name, verloop_b1, dbstep_b1):
		tolerance = 0.01
		self.compare_with_tolerance(verloop_b1, dbstep_b1, tolerance)

	@pytest.mark.parametrize("name, verloop_b5, dbstep_b5", sav.file_names_and_b5s)
	def test_b5_against_verloop(self, name, verloop_b5, dbstep_b5):
		tolerance = 0.01
		self.compare_with_tolerance(verloop_b5, dbstep_b5, tolerance)

	@staticmethod
	def compare_with_tolerance(value1, value2, tolerance):
		"""Compares to values and asserts whether they are +/- tolerance away from each other.

			NOTE: I didn't use np/math.isclose because I like it to print how far
			values are from each other on failure.
		"""
		assert np.round(abs(value1 - value2), 2) <= tolerance
