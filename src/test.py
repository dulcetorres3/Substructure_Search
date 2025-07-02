from molecule import MoleculeScreen
import pytest


def test_constructor_types():
    """Test that nodes and bonds must be within a list container."""

    with pytest.raises(TypeError):
        molecule = MoleculeScreen((1, 2, 3), [(1, 2, 1), (2, 3, 1)])
    with pytest.raises(TypeError):
        molecule = MoleculeScreen((1, 2, 3), ((1, 2, 1), (2, 3, 1)))


def test_alternative_constructor():
    """Test that an error is raised when a file is not found"""

    with pytest.raises(FileNotFoundError):
        molecule = MoleculeScreen.from_sdf("epoxide.sdf")


@pytest.mark.parametrize(
    ("filepath"),
    [
        pytest.param("../test_compounds/benzofuran.sdf", id="benzofuran"),
        pytest.param("../test_compounds/ibuprofen.sdf", id="ibuprofen"),
    ],
)
def test_instantiation(filepath):
    """
    Tests that Molecule instances are properly constructed from filepaths to
    .sdf files.
    """
    mol = MoleculeScreen.from_sdf(filepath)
    assert isinstance(mol, MoleculeScreen)
    assert hasattr(mol, "nodes")
    assert hasattr(mol, "bonds")
    assert hasattr(mol, "graph")
    assert hasattr(mol, "color_dict")


@pytest.mark.parametrize(
    ("compound, expected"),
    [
        pytest.param("../test_compounds/benzene.sdf", True, id="sp2"),
        pytest.param("../test_compounds/cyclohexane.sdf", False, id="not_sp2"),
    ],
)
def test_sp2(compound, expected):
    """Tests that the hybdrization of carbon atoms are appropriately determined."""
    test_compound = MoleculeScreen.from_sdf(compound, True)
    assert test_compound.hetero_sp2_hybridized("C1") == expected


@pytest.mark.parametrize(
    ("compound, expected"),
    [
        pytest.param("../test_compounds/benzofuran.sdf", True, id="sp2"),
        pytest.param("../test_compounds/oxetane.sdf", False, id="not_sp2"),
    ],
)
def test_sp2_hetero(compound, expected):
    """Tests that the hybridization of hetero atoms are appropriately determined."""
    test_compound = MoleculeScreen.from_sdf(compound, True)
    assert test_compound.hetero_sp2_hybridized("O1") == expected


@pytest.mark.parametrize(
    ("compound, electrons, expected"),
    [
        pytest.param("../test_compounds/oxetane.sdf", 6, True, id="huckel"),
        pytest.param("../test_compounds/oxetane.sdf", 9, False, id="not_huckel"),
    ],
)
def test_huckel_electrons(compound, electrons, expected):
    """Tests huckel values are succesfully determined."""
    test_compound = MoleculeScreen.from_sdf(compound, True)
    assert test_compound.huckel_electrons(electrons) == expected


@pytest.mark.parametrize(
    ("compound, expected"),
    [
        pytest.param("../test_compounds/benzofuran.sdf", True, id="aromatic"),
        pytest.param("../test_compounds/cyclohexane.sdf", False, id="not_aromatic"),
    ],
)
def test_aromatic_ring_detections(compound, expected):
    """Tests that aromatic rings are appropriately determined."""
    test_compound = MoleculeScreen.from_sdf(compound, True)
    assert test_compound.aromatic_ring_detection() == expected


@pytest.mark.parametrize(
    ("compound, substructure, expected"),
    [
        pytest.param(
            "../test_compounds/benzofuran.sdf",
            "../test_compounds/benzene.sdf",
            True,
            id="hit",
        ),
        pytest.param(
            "../test_compounds/benzofuran.sdf",
            "../test_compounds/methyl.sdf",
            False,
            id="no_hit",
        ),
    ],
)
def test_screen(compound, substructure, expected):
    """Tests that substructures are correctly identified within a compound."""
    test_compound = MoleculeScreen.from_sdf(compound, True)
    test_substructure = MoleculeScreen.from_sdf(substructure, True)
    assert test_compound.screen(test_substructure) == expected


@pytest.mark.parametrize(
    ("compound_1, compound_2, expected"),
    [
        pytest.param(
            "../test_compounds/benzofuran.sdf",
            "../test_compounds/benzofuran.sdf",
            True,
            id="equal",
        ),
        pytest.param(
            "../test_compounds/benzofuran.sdf",
            "../test_compounds/benzene.sdf",
            False,
            id="not_equal",
        ),
    ],
)
def test_eq_operator(compound_1, compound_2, expected):
    """
    Tests the equivalence operator on different Molecule instances based on
    their molecular fingerprint.
    """
    mol1 = MoleculeScreen.from_sdf(compound_1)
    mol2 = MoleculeScreen.from_sdf(compound_2)
    assert (mol1 == mol2) == expected
