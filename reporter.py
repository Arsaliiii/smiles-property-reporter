import argparse
import sys
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors


def get_properties(smiles: str) -> dict:
    if not isinstance(smiles, str) or not smiles.strip():
        raise ValueError("SMILES string cannot be empty")

    mol = Chem.MolFromSmiles(smiles)
    if mol is None or mol.GetNumAtoms() == 0:
        raise ValueError(f"Invalid SMILES string: {smiles!r}")

    return {
        "canonical_smiles": Chem.MolToSmiles(mol),
        "formula": rdMolDescriptors.CalcMolFormula(mol),
        "molecular_weight": Descriptors.MolWt(mol),
        "logp": Descriptors.MolLogP(mol),
        "tpsa": Descriptors.TPSA(mol),
        "h_donors": Lipinski.NumHDonors(mol),
        "h_acceptors": Lipinski.NumHAcceptors(mol),
        "rotatable_bonds": Lipinski.NumRotatableBonds(mol),
        "heavy_atoms": mol.GetNumHeavyAtoms(),
        "rings": rdMolDescriptors.CalcNumRings(mol),
    }

def lipinski_check(props: dict) -> list[str]:
    violations: list[str] = []

    if props["molecular_weight"] > 500:
        violations.append("MW > 500")
    if props["logp"] > 5:
        violations.append("logP > 5")
    if props["h_donors"] > 5:
        violations.append("H-bond donors > 5")
    if props["h_acceptors"] > 10:
        violations.append("H-bond acceptors > 10")

    return violations

def fmt2(value: float) -> str:
    return f"{round(value, 2) + 0.0:.2f}"

def passes_lipinski(violations: list[str]) -> bool:
    return len(violations) <= 1

def main():
    parser = argparse.ArgumentParser(
        description="Compute molecular properties and Lipinski Rule of 5 summary from a SMILES string."
    )
    parser.add_argument("smiles", help="SMILES string of the molecule")
    args = parser.parse_args()

    try:
        props = get_properties(args.smiles)
        violations = lipinski_check(props)
        passed = passes_lipinski(violations)
    except ValueError as err:
        print(f"Error: {err}", file=sys.stderr)
        sys.exit(1)

    status = "PASS" if passed else "FAIL"
    violations_str = ", ".join(violations) if violations else "None"

    print(f"SMILES:           {args.smiles}")
    print(f"Canonical SMILES: {props['canonical_smiles']}")
    print(f"Formula:          {props['formula']}")
    print(f"Molecular weight: {fmt2(props['molecular_weight'])}")
    print(f"LogP:             {fmt2(props['logp'])}")
    print(f"TPSA:             {fmt2(props['tpsa'])}")
    print(f"H-bond donors:    {props['h_donors']}")
    print(f"H-bond acceptors: {props['h_acceptors']}")
    print(f"Rotatable bonds:  {props['rotatable_bonds']}")
    print(f"Heavy atoms:      {props['heavy_atoms']}")
    print(f"Rings:            {props['rings']}")
    print(f"Lipinski status:  {status} (violations: {violations_str})")

if __name__ == "__main__":
    main()