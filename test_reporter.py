import pytest
from reporter import get_properties, lipinski_check, passes_lipinski

def test_get_properties():
    # Ethanol
    ethanol = get_properties("CCO")
    assert ethanol["formula"] == "C2H6O"
    assert ethanol["canonical_smiles"] == "CCO"
    assert ethanol["molecular_weight"] == pytest.approx(46.07, abs=0.01)
    assert ethanol["heavy_atoms"] == 3

    # Aspirin
    aspirin = get_properties("CC(=O)Oc1ccccc1C(=O)O")
    assert aspirin["formula"] == "C9H8O4"
    assert aspirin["molecular_weight"] == pytest.approx(180.16, abs=0.01)
    assert aspirin["rings"] == 1
    assert aspirin["h_donors"] == 1
    assert aspirin["h_acceptors"] == 3

    # Canonical SMILES equivalence
    assert get_properties("OCC")["canonical_smiles"] == get_properties("CCO")["canonical_smiles"]

    # Invalid / empty inputs
    for bad_smiles in ["", "   ", "C(C)(C)(C)(C)C"]:
        with pytest.raises(ValueError):
            get_properties(bad_smiles)


def test_lipinski_check():
    # Aspirin passes all rules
    aspirin_props = get_properties("CC(=O)Oc1ccccc1C(=O)O")
    assert lipinski_check(aspirin_props) == []

    # C40 alkane violates both MW > 500 and logP > 5
    alkane_40 = get_properties("C" * 40)
    alkane_violations = lipinski_check(alkane_40)
    assert "MW > 500" in alkane_violations
    assert "logP > 5" in alkane_violations

    # Hand-built dict with donors = 6 and benign values elsewhere
    high_donors_props = {
        "molecular_weight": 250.0,
        "logp": 2.0,
        "h_donors": 6,
        "h_acceptors": 4,
    }
    assert lipinski_check(high_donors_props) == ["H-bond donors > 5"]

    # MW boundary: 500.0 is allowed, 500.01 violates MW > 500
    at_limit_props = {
        "molecular_weight": 500.0,
        "logp": 2.0,
        "h_donors": 1,
        "h_acceptors": 2,
    }
    assert lipinski_check(at_limit_props) == []

    just_over_props = {
        "molecular_weight": 500.01,
        "logp": 2.0,
        "h_donors": 1,
        "h_acceptors": 2,
    }
    assert lipinski_check(just_over_props) == ["MW > 500"]


def test_passes_lipinski():
    assert passes_lipinski([]) is True
    assert passes_lipinski(["MW > 500"]) is True
    assert passes_lipinski(["MW > 500", "logP > 5"]) is False