#!/usr/bin/env python
# n.b. chmod 755

# TODO: Give ekstra kwd through gaussian - use to define solvent.

import sys
import os
import fortranformat as ff
import numpy as np

from rdkit import Chem

__XTB_PATH__ = os.environ.get("XTB_PATH", None)
if __XTB_PATH__ is None:
    raise RuntimeError("Path to xtb executable not set via environment variable 'XTB_PATH'.")


def parse_ifile(ifile):
    """ parse the ifile from xTB"""

    pt = Chem.GetPeriodicTable()

    with open(ifile, "r") as ifile:
        gau_input = ifile.readlines()

    tokens = gau_input[0].split()

    natoms = int(tokens[0])
    nderiv = int(tokens[1])
    chrg = int(tokens[2])
    spin = int(tokens[3])

    coords = np.empty((natoms, 3))
    atomtypes = []
    for i, line in enumerate(gau_input[1 : 1 + natoms]):
        line = line.split()
        atomtypes.append(pt.GetElementSymbol(int(line[0])))
        coords[i] = np.array(list(map(float, line[1 : 1 + 3]))) * 0.529177249
    return natoms, nderiv, chrg, spin, atomtypes, coords


def parse_ofile(ofile, energy, natoms, dipole, gradient=None, hessian=None):
    """
    parse that outfule for Gaussian to read.
    """
    headformat = ff.FortranRecordWriter("4D20.12")
    bodyformat = ff.FortranRecordWriter("3D20.12")

    f = open(ofile, "w")
    head = [energy, dipole[0], dipole[1], dipole[2]]
    headstring = headformat.write(head)
    f.write(headstring + "\n")

    if gradient is None:
        gradient = np.zeros((natoms, 3))

    for i in range(natoms):
        output = bodyformat.write(gradient[i])
        f.write(output + "\n")

    # polarizability and dipole derivatives are set to zero
    polarizability = np.zeros((2, 3))
    dipole_derivative = np.zeros((3 * natoms, 3))

    for i in range(2):
        output = bodyformat.write(polarizability[i])
        f.write(output + "\n")

    for i in range(3 * natoms):
        output = bodyformat.write(dipole_derivative[i])
        f.write(output + "\n")

    if hessian is not None:  # Only needed if ndreiv = 1
        #print(hessian.shape)
        tril = np.tril_indices(hessian.shape[0])
        tril_hessian = hessian[tril]
        for window in tril_hessian.reshape(int(tril_hessian.shape[0] / 3), 3):
            output = bodyformat.write(window)
            f.write(output + "\n")

    f.close()


def write_xyz(natoms, atomtypes, coords):
    """ Write .xyz file """
    xyz = f"{natoms} \n \n"
    for atomtype, coord in zip(atomtypes, coords):
        xyz += f"{atomtype}  {' '.join(list(map(str, coord)))} \n"

    with open("mol.xyz", "w") as inp:
        inp.write(xyz)


def get_energy(output):
    """ """
    for line in output.split("\n"):
        if "TOTAL ENERGY" in line:
            return float(line.split()[3])


def get_dipole(output):
    """ """
    tmp_output = output.split("molecular dipole:")
    del tmp_output[0]

    dipole_data = tmp_output[0].split("\n")
    dipole_line = dipole_data[3].split()

    dipole = np.array([dipole_line[1], dipole_line[2], dipole_line[3]], dtype=float)
    return dipole


def get_gradient(natoms):
    """ """
    with open("gradient") as grad_file:
        gradient_file = grad_file.read()

    gradient_data = gradient_file.split("\n")
    del gradient_data[: 2 + natoms]

    gradient = np.empty((natoms, 3))
    for i, line in enumerate(gradient_data[:natoms]):
        line = line.split()
        gradient[i, :] = line
    return gradient


def get_hessian(natoms):
    """ """
    hess_file = open("hessian", "r")

    i = 0
    hessian = np.empty(3 * natoms * 3 * natoms) 
    for line in hess_file:
        if "$hessian" in line:
            continue

        for elm in line.strip().split():
            hessian[i] = float(elm)
            i += 1
    
    hess_file.close()

    hessian = hessian.reshape((3 * natoms, 3 * natoms))
    return hessian


def run_xtb(natoms, nderiv, chrg, spin, atomtypes, coords, solvent=None):
    """ """

    write_xyz(natoms, atomtypes, coords)

    os.environ["OMP_NUM_THREADS"] = str(2)
    os.environ["MKL_NUM_THREADS"] = str(2)
    cmd = f"xtb mol.xyz --chrg {chrg} --uhf {spin - 1} --gfn 2 "
    if nderiv == 1:
        cmd += "--grad "
    elif nderiv == 2:
        cmd += "--hess --grad "
    
    if solvent is not None:
        method, solvent = solvent.split('=')
        cmd += f"--{method} {solvent} "

    print(cmd)
    with os.popen(cmd) as xtb:
        output = xtb.read()

    energy = get_energy(output)
    dipole = get_dipole(output)

    return energy, dipole


def clean_dir():
    """ delete all files """
    files = [
        "energy",
        "charges",
        "mol.xyz",
        "xtbrestart",
        "gradient",
        "hessian",
        "vibspectrum",
        "wbo",
        "mol.engrad",
        "xtbtopo.mol",
        "g98.out",
        "xtbhess.xyz",
    ]
    for _file in files:
        if os.path.exists(_file):
            os.remove(_file)


if __name__ == "__main__":
    
    solvent = None

    if len(sys.argv[1]) > 7: # given ekstra kwd
        for i, kwd in enumerate(sys.argv):
            if kwd == "R":
                break
             
            if "gbsa" or "alpb" in kwd:
                solvent = kwd
        
        ifile=sys.argv[i+1]    
        ofile=sys.argv[i+2]
    
    else:
        ifile = sys.argv[2]
        ofile = sys.argv[3]

    (natoms, nderiv, chrg, spin, atomtypes, coords) = parse_ifile(ifile)
    energy, dipole = run_xtb(natoms, nderiv, chrg, spin, atomtypes, coords, solvent=solvent)

    if nderiv == 0:
        parse_ofile(ofile, energy, natoms, dipole)
    elif nderiv == 1:
        gradient = get_gradient(natoms)
        parse_ofile(ofile, energy, natoms, dipole, gradient=gradient)
    elif nderiv == 2:
        gradient = get_gradient(natoms)
        hessian = get_hessian(natoms)
        parse_ofile(ofile, energy, natoms, dipole, gradient=gradient, hessian=hessian)

    clean_dir()
