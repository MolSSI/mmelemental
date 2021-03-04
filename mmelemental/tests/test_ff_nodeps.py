from mmelemental.models import forcefield as ff
import numpy
import pytest
import os

natoms = 10

def rewrite(filename):
    import json

    with open(filename, "r") as fp:
        data = json.load(fp)

    with open(filename, "w") as fp:
        json.dump(data, fp, indent=4)


def rand_angles(npts):
    return numpy.random.rand(npts) * numpy.pi

###########################################
####### NON-BONDED ########################
def test_nonbonded():
    lj = ff.nonbonded.potentials.LennardJones(
        epsilon=numpy.random.rand(natoms), sigma=numpy.random.rand(natoms)
    )
    nonbonded = ff.nonbonded.NonBonded(params=lj)
    assert nonbonded.form == "LennardJones"
    return nonbonded


############################################
######## BONDS #############################
def test_bonds_linear():
    linear = ff.bonded.bonds.potentials.Harmonic(
        spring=numpy.random.rand(natoms), lengths=numpy.random.rand(natoms)
    )
    bonds = ff.bonded.Bonds(params=linear)
    assert bonds.form == "Harmonic"
    return bonds


def test_bonds_hybrid():
    natoms_by_2 = int(natoms / 2)
    lengths = numpy.random.rand(natoms)
    linear = ff.bonded.bonds.potentials.Harmonic(
        spring=numpy.random.rand(natoms_by_2), lengths=lengths[:natoms_by_2]
    )
    gromos = ff.bonded.bonds.potentials.Gromos96(
        spring=numpy.random.rand(natoms_by_2), lengths=lengths[natoms_by_2:]
    )

    bonds = ff.bonded.Bonds(params=[linear, gromos])
    assert bonds.form == ["Harmonic", "Gromos96"]
    return bonds


###########################################
####### ANGLES ############################
def test_angles():
    linear = ff.bonded.angles.potentials.Harmonic(
        spring=numpy.random.rand(natoms),
        angles=rand_angles(natoms),
        angles_units="radians",
    )
    angles = ff.bonded.Angles(params=linear)
    assert angles.form == "Harmonic"
    return angles


###########################################
####### DIHEDRALS #########################
def test_dihedrals():
    linear = ff.bonded.dihedrals.potentials.Harmonic(
        spring=numpy.random.rand(natoms),
        angles=rand_angles(natoms),
        angles_units="radians",
    )
    dihedrals = ff.bonded.Dihedrals(params=linear)
    assert dihedrals.form == "Harmonic"
    return dihedrals


def test_forcefield():
    nonbonded = test_nonbonded()
    bonds = test_bonds_linear()
    angles = test_angles()
    dihedrals = test_dihedrals()
    mm_ff = ff.ForceField(
        nonbonded=nonbonded,
        bonds=bonds,
        angles=angles,
        dihedrals=dihedrals,
        charges=numpy.random.rand(natoms),
    )
    mm_ff.to_file("forcefield.json")
    rewrite("forcefield.json")
    os.remove("forcefield.json")


def test_forcefield_hybrid():
    nonbonded = test_nonbonded()
    bonds = test_bonds_hybrid()
    angles = test_angles()
    dihedrals = test_dihedrals()
    mm_ff = ff.ForceField(
        nonbonded=nonbonded,
        bonds=bonds,
        angles=angles,
        dihedrals=dihedrals,
        charges=numpy.random.rand(natoms),
    )
    mm_ff.to_file("forcefield.json")
    rewrite("forcefield.json")
