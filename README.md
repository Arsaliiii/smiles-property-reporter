# SMILES Property Reporter

A small command-line tool that takes a molecule written as a SMILES string and prints a report of its basic properties, followed by a Lipinski rule-of-five check.

## Description

Given a SMILES string such as `CC(=O)Oc1ccccc1C(=O)O` (aspirin), the program uses [RDKit](https://www.rdkit.org/) to report:

- the canonical SMILES
- the molecular formula
- the molecular weight
- logP (how oily or water-loving the molecule is)
- TPSA (topological polar surface area)
- the number of hydrogen-bond donors and acceptors
- the number of rotatable bonds, heavy atoms and rings
- a Lipinski rule-of-five result

Lipinski's rule of five is a rough screen for oral drug-likeness. A molecule violates a rule if its molecular weight is above 500, its logP is above 5, it has more than 5 H-bond donors, or it has more than 10 H-bond acceptors. Here, a molecule passes with 1 violation or fewer.

Invalid or empty input gives a single clear error message instead of a report.

## Requirements

- Python 3
- RDKit
- pytest (for the tests)

With conda:

```
conda install -c conda-forge rdkit pytest
```

## Usage

```
python reporter.py "SMILES"
```

Put the SMILES in quotes. Characters such as `(`, `)`, `=` and `#` are special in the terminal and can break the command without quotes.

### Example 1: aspirin (pass)

```
python reporter.py "CC(=O)Oc1ccccc1C(=O)O"
```
```
SMILES:           CC(=O)Oc1ccccc1C(=O)O
Canonical SMILES: CC(=O)Oc1ccccc1C(=O)O
Formula:          C9H8O4
Molecular weight: 180.16
LogP:             1.31
TPSA:             63.60
H-bond donors:    1
H-bond acceptors: 3
Rotatable bonds:  2
Heavy atoms:      13
Rings:            1
Lipinski status:  PASS (violations: None)
```

### Example 2: a 40-carbon alkane (fail)

```
python reporter.py "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
```
```
SMILES:           CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
Canonical SMILES: CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
Formula:          C40H82
Molecular weight: 563.10
LogP:             15.85
TPSA:             0.00
H-bond donors:    0
H-bond acceptors: 0
Rotatable bonds:  37
Heavy atoms:      40
Rings:            0
Lipinski status:  FAIL (violations: MW > 500, logP > 5)
```

### Example 3: invalid input

```
python reporter.py "C(C)(C)(C)(C)C"
```
```
Error: Invalid SMILES string: 'C(C)(C)(C)(C)C'
```

The carbon in that SMILES has five bonds, which is not allowed. RDKit may also print its own valence message before the error line. An empty string gives `Error: SMILES string cannot be empty`.

## How it works

### `get_properties(smiles)`
Parses the SMILES with `Chem.MolFromSmiles` and returns a dictionary of the properties listed above. It raises `ValueError` for an empty or whitespace-only string, and for a SMILES that RDKit cannot parse (RDKit returns `None` in that case). Values are returned unrounded, and rounding happens only when printing.

### `lipinski_check(props)`
Takes the dictionary from `get_properties` and returns a list of the rules that were violated, for example `["MW > 500", "logP > 5"]`. An empty list means no violations. The comparisons are strict, so a value exactly at the limit is not a violation.

### `passes_lipinski(violations)`
Returns `True` if the list has one violation or fewer, and `False` otherwise.

### `main()`
Reads the SMILES from the command line with `argparse`, computes everything first, and then prints the report. If the input is invalid, it prints one `Error:` line to stderr and exits, so a failed run never prints a half-finished report.

## Limitations

- The rule of five is a rough screen. It does not predict whether a molecule is actually a good drug.
- It handles one molecule at a time.
- The logP value is RDKit's Crippen estimate, not a measured value.

## Tests

`test_reporter.py` has one test for each of `get_properties`, `lipinski_check` and `passes_lipinski`, including the invalid inputs and the exact limits of the Lipinski rules. Run them with:

```
pytest
```
