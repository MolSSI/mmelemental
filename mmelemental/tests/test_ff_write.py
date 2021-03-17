from mmelemental.models import forcefield as ff
import numpy
import pytest
import os

natoms = 10
nbonds = natoms - 1
nangles = natoms - 2
ndihedrals = natoms - 3


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
    nonbonded = ff.nonbonded.NonBonded(params=lj, form="LennardJones")
    assert nonbonded.form == "LennardJones"
    return nonbonded


############################################
######## BONDS #############################
def test_bonds_linear():
    print(ff.bonded.Bonds.supported_potentials())
    linear = ff.bonded.bonds.potentials.Harmonic(spring=numpy.random.rand(nbonds))
    bonds = ff.bonded.Bonds(
        params=linear, lengths=numpy.random.rand(nbonds), form="Harmonic"
    )
    assert bonds.form == "Harmonic"
    return bonds


def test_bonds_gromos96():
    lengths = numpy.random.rand(nbonds)
    gromos = ff.bonded.bonds.potentials.Gromos96(
        spring=numpy.random.rand(nbonds),
    )

    bonds = ff.bonded.Bonds(params=gromos, lengths=lengths, form="Gromos96")
    assert bonds.form == "Gromos96"
    return bonds


###########################################
####### ANGLES ############################
def test_angles():
    linear = ff.bonded.angles.potentials.Harmonic(
        spring=numpy.random.rand(nangles),
    )
    angles = ff.bonded.Angles(
        params=linear,
        angles=rand_angles(nangles),
        angles_units="radians",
        form="Harmonic",
    )
    assert angles.form == "Harmonic"
    return angles


###########################################
####### DIHEDRALS #########################
def test_dihedrals():
    linear = ff.bonded.dihedrals.potentials.Harmonic(
        spring=numpy.random.rand(ndihedrals),
    )
    dihedrals = ff.bonded.Dihedrals(
        params=linear,
        angles=rand_angles(ndihedrals),
        angles_units="radians",
        form="Harmonic",
    )
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
    bonds = [test_bonds_linear(), test_bonds_gromos96()]
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
