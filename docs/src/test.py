import mmelemental

angles = mmelemental.models.forcefield.Angles(
    **{
        "form": "Harmonic",
        "params": {"spring": [0.09565291598722435]},
        "angles": [104.52],
        "connectivity": [
            [0, 1, 2],
        ],
    }
)

bonds = mmelemental.models.forcefield.Bonds(
    **{
        "form": "Harmonic",
        "params": {"spring": [2512.08, 2512.08]},
        "lengths": [0.9572, 0.9572],
        "connectivity": [[0, 1, 1.0], [0, 2, 1.0]],
    }
)

nb = mmelemental.models.forcefield.NonBonded(
    **{
        "form": "LennardJones",
        "params": {"epsilon": [0.636386, 0.0, 0.0], "sigma": [1.575305, 0.0, 0.0]},
    }
)

ff = mmelemental.models.ForceField(
    **{
        "symbols": ["H", "H", "O"],
        "charges": [-0.834, 0.417, 0.417],
        "masses": [16.0, 1.008, 1.008],
        "exclusions": "3",
        "defs": ["OW", "HW", "HW"],
        "nonbonded": nb,
        "bonds": bonds,
        "angles": angles,
    }
)

print(ff)
