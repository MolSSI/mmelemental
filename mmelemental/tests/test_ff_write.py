from mmelemental.models import forcefield as ff
import numpy
import pytest
import os

natoms = 10
nbonds = natoms - 1
nangles = natoms - 2
ndihedrals = natoms - 3

# Data
rand_npts = numpy.random.rand(natoms)
rand_nbond = numpy.random.rand(nbonds)
rand_nangle = numpy.random.rand(nangles)
rand_ndih = numpy.random.rand(ndihedrals)
rand_int_ndih = numpy.random.randint(ndihedrals, size=ndihedrals)

# Pytest parametrization args
nb_pytest_params = "form,params", [
    ("LennardJones", {"epsilon": rand_npts, "sigma": rand_npts})
]
bonds_pytest_params = "form,lengths,params", [
    ("Harmonic", rand_nbond, {"spring": rand_nbond}),
    ("Gromos96", rand_nbond, {"spring": rand_nbond}),
]
angles_pytest_params = "form,lengths,params", [
    ("Harmonic", rand_nangle, {"spring": rand_nangle})
]
dihedrals_pytest_params = "form,lengths,params", [
    (
        "Charmm",
        rand_ndih,
        {"energy": rand_ndih, "periodicity": rand_int_ndih, "phase": rand_ndih},
    ),
    (
        "CharmmMulti",
        rand_ndih,
        {
            "energy": (rand_ndih, rand_ndih),
            "periodicity": (rand_int_ndih, rand_int_ndih),
            "phase": (rand_ndih, rand_ndih),
        },
    ),
]
dihedrals_improper_pytest_params = "form,lengths,params", [
    (
        "Charmm",
        rand_ndih,
        {"energy": rand_ndih, "periodicity": rand_int_ndih, "phase": rand_ndih},
    ),
    (
        "CharmmMulti",
        rand_ndih,
        {
            "energy": (rand_ndih, rand_ndih),
            "periodicity": (rand_int_ndih, rand_int_ndih),
            "phase": (rand_ndih, rand_ndih),
        },
    ),
]


def rewrite(filename):
    import json

    with open(filename, "r") as fp:
        data = json.load(fp)

    with open(filename, "w") as fp:
        json.dump(data, fp, indent=4)


def rand_angles(npts):
    return rand_data * numpy.pi


###########################################
####### NON-BONDED ########################
@pytest.mark.parametrize(*nb_pytest_params)
def test_nonbonded(form, params):
    potential = getattr(ff.nonbonded.potentials, form)(**params)
    nonbonded = ff.nonbonded.NonBonded(params=potential, form=form)
    assert nonbonded.form == form
    return nonbonded


############################################
######## BONDS #############################
@pytest.mark.parametrize(*bonds_pytest_params)
def test_bonds(form, lengths, params):
    indices = [(i, i + 1, 1.0) for i in range(nbonds)]
    potential = getattr(ff.bonded.bonds.potentials, form)(**params)
    bonds = ff.bonded.Bonds(
        params=potential,
        lengths=lengths,
        form=form,
        connectivity=indices,
    )
    assert bonds.form == form
    return bonds


###########################################
####### ANGLES ############################
@pytest.mark.parametrize(*angles_pytest_params)
def test_angles(form, lengths, params):
    indices = [(i, i + 1, i + 2) for i in range(nangles)]
    potential = getattr(ff.bonded.angles.potentials, form)(**params)
    angles = ff.bonded.Angles(
        params=potential,
        angles=lengths,
        angles_units="radians",
        connectivity=indices,
        form=form,
    )
    assert angles.form == form
    return angles


###########################################
####### DIHEDRALS #########################
@pytest.mark.parametrize(*dihedrals_pytest_params)
def test_dihedrals(form, lengths, params):
    indices = [(i, i + 1, i + 2, i + 3) for i in range(ndihedrals)]
    potential = getattr(ff.bonded.dihedrals.potentials, form)(**params)
    dihedrals = ff.bonded.Dihedrals(
        params=potential,
        angles=lengths,
        angles_units="radians",
        connectivity=indices,
        form=form,
    )
    assert dihedrals.form == form
    return dihedrals


###########################################
####### IMPROPER DIHEDRALS ################
@pytest.mark.parametrize(*dihedrals_pytest_params)
def test_improper_dihedrals(form, lengths, params):
    indices = [(i, i + 1, i + 2, i + 3) for i in range(ndihedrals)]
    potential = getattr(ff.bonded.dihedrals_improper.potentials, form)(**params)
    dihedrals_improper = ff.bonded.DihedralsImproper(
        params=potential,
        angles=lengths,
        angles_units="radians",
        connectivity=indices,
        form=form,
    )
    assert dihedrals_improper.form == form
    return dihedrals_improper


@pytest.mark.parametrize(*dihedrals_pytest_params)
def test_forcefield(form, lengths, params):
    nonbonded = test_nonbonded(
        form="LennardJones", params={"epsilon": rand_npts, "sigma": rand_npts}
    )
    bonds = test_bonds(
        form="Harmonic", lengths=rand_nbond, params={"spring": rand_nbond}
    )
    angles = test_angles(
        form="Harmonic", lengths=rand_nangle, params={"spring": rand_nangle}
    )
    dihedrals = test_dihedrals(form, lengths, params)
    dihedrals_improper = test_improper_dihedrals(form, lengths, params)
    mm_ff = ff.ForceField(
        nonbonded=nonbonded,
        bonds=bonds,
        angles=angles,
        dihedrals=dihedrals,
        dihedrals_improper=dihedrals_improper,
        charges=numpy.random.rand(natoms),
        symbols=["H" for _ in range(natoms)],
    )
    mm_ff.to_file("forcefield-single.json")
    rewrite("forcefield-single.json")
    os.remove("forcefield-single.json")


@pytest.mark.skip("Cannot run this now")
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
        symbols=["H" for _ in range(natoms)],
    )
    mm_ff.to_file("forcefield.json")
    rewrite("forcefield.json")
    os.remove("forcefield.json")
