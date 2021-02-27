from mmelemental.models import forcefield as ff
import numpy
import pytest
import os

def rewrite(filename):
    import json

    with open(filename, "r") as fp:
        data = json.load(fp)

    with open(filename, "w") as fp:
        json.dump(data, fp, indent=4)

def rand_angles(npts):
    return numpy.random.rand(npts) * numpy.pi

natoms = 10

###########################################
####### NON-BONDED ########################
def test_nonbonded():
    lj = ff.nonbonded.potentials.LennardJones(
        epsilon=numpy.random.rand(natoms), sigma=numpy.random.rand(natoms)
    )
    nonbonded = ff.nonbonded.NonBonded(params=lj, charges=numpy.random.rand(natoms))
    assert nonbonded.form == "LennardJones"
    return nonbonded

############################################
######## BONDS #############################
def test_bonds():
    linear = ff.bonded.bonds.potentials.Harmonic(spring=numpy.random.rand(natoms))
    bonds = ff.bonded.Bonds(params=linear, lengths=numpy.random.rand(natoms))
    assert bonds.form == "Harmonic"
    return bonds

###########################################
####### ANGLES ############################
def test_angles():
    linear = ff.bonded.angles.potentials.Harmonic(spring=numpy.random.rand(natoms))
    angles = ff.bonded.Angles(params=linear, angles=rand_angles(natoms), angles_units="radians")
    assert angles.form == "Harmonic"
    return angles

###########################################
####### DIHEDRALS #########################
def test_dihedrals():
    linear = ff.bonded.dihedrals.potentials.Harmonic(spring=numpy.random.rand(natoms))
    dihedrals = ff.bonded.Dihedrals(params=linear, angles=rand_angles(natoms), angles_units="radians")
    assert dihedrals.form == "Harmonic"
    return dihedrals

def test_forcefield():
    nonbonded = test_nonbonded()
    bonds = test_bonds()
    angles = test_angles()
    dihedrals = test_dihedrals()
    mm_ff = ff.ForceField(nonbonded=nonbonded, bonds=bonds, angles=angles, dihedrals=dihedrals)
    mm_ff.to_file("forcefield.json")
    rewrite("forcefield.json")
    os.remove("forcefield.json")
