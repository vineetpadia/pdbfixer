"""Offline coverage for matching observed chains to their complete sequences."""

from openmm import app
import pytest

from pdbfixer.pdbfixer import PDBFixer, Sequence


def make_fixer(chains, sequences):
    fixer = PDBFixer.__new__(PDBFixer)
    fixer.topology = app.Topology()
    for chain_id, residues in chains:
        chain = fixer.topology.addChain(chain_id)
        for name, residue_id, insertion_code in residues:
            fixer.topology.addResidue(name, chain, residue_id, insertionCode=insertion_code)
    fixer.sequences = [Sequence(chain_id, residues) for chain_id, residues in sequences]
    return fixer


@pytest.mark.parametrize('blank_code', ['', ' '])
def test_missing_residue_with_blank_insertion_code(blank_code):
    fixer = make_fixer(
        [('A', [('ALA', '1', blank_code), ('SER', '3', blank_code)])],
        [('A', ['ALA', 'GLY', 'SER'])],
    )
    fixer.findMissingResidues()
    assert fixer.missingResidues == {(0, 1): ['GLY']}


def test_multiple_insertion_codes_preserve_gaps_and_terminal_substitutions():
    fixer = make_fixer(
        [('A', [('ALA', '7', ''), ('GLY', '8', ' '), ('SER', '8', 'A'),
                ('THR', '8', 'B'), ('VAL', '10', '')])],
        [('A', ['ASN', 'ALA', 'GLY', 'SER', 'THR', 'LEU', 'VAL', 'MSE'])],
    )
    fixer.findMissingResidues()
    assert fixer.missingResidues == {(0, 0): ['ASN'], (0, 4): ['LEU'], (0, 5): ['MET']}


def test_internal_gaps_and_missing_termini():
    fixer = make_fixer(
        [('A', [('ALA', '4', ''), ('THR', '6', '')])],
        [('A', ['SER', 'ALA', 'GLY', 'THR', 'VAL'])],
    )
    fixer.findMissingResidues()
    assert fixer.missingResidues == {(0, 0): ['SER'], (0, 1): ['GLY'], (0, 2): ['VAL']}


def test_first_matching_offset_is_used():
    fixer = make_fixer(
        [('A', [('ALA', '1', ''), ('GLY', '2', '')])],
        [('A', ['ALA', 'GLY', 'ALA', 'GLY'])],
    )
    fixer.findMissingResidues()
    assert fixer.missingResidues == {(0, 2): ['ALA', 'GLY']}


@pytest.mark.parametrize('sequence_chain,sequence', [
    ('A', ['ALA', 'VAL']),
    ('A', ['ALA']),
    ('B', ['ALA', 'GLY']),
])
def test_unmatched_sequences_clear_previous_results(sequence_chain, sequence):
    fixer = make_fixer(
        [('A', [('ALA', '1', ''), ('GLY', '2', '')])],
        [(sequence_chain, sequence)],
    )
    fixer.missingResidues = {(0, 0): ['SER']}
    fixer.findMissingResidues()
    assert fixer.missingResidues == {}


def test_duplicate_chain_ids_use_first_unassigned_match_and_keep_chain_indices():
    fixer = make_fixer(
        [('A', []), ('A', [('ALA', '1', '')]), ('A', [('ALA', '1', '')])],
        [('A', ['GLY', 'ALA']), ('A', ['ALA', 'SER']), ('A', ['ALA', 'TRP'])],
    )
    fixer.findMissingResidues()
    assert fixer.missingResidues == {(1, 0): ['GLY'], (2, 1): ['SER']}


def test_duplicate_chain_ids_skip_nonmatching_chains():
    fixer = make_fixer(
        [('A', [('GLY', '1', '')]), ('A', [('ALA', '1', '')])],
        [('A', ['VAL', 'ALA']), ('A', ['GLY', 'SER'])],
    )
    fixer.findMissingResidues()
    assert fixer.missingResidues == {(0, 1): ['SER'], (1, 0): ['VAL']}
