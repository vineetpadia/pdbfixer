"""Regression tests for preparing hydrogen variants without minimization."""

import pytest
import openmm as mm
from openmm import app, unit

from pdbfixer import PDBFixer


def _make_fixer(source, include_oxygen=True, reverse_definition_bonds=False):
    fixer = PDBFixer.__new__(PDBFixer)
    fixer.topology = app.Topology()
    chain = fixer.topology.addChain('A')
    residue = fixer.topology.addResidue('XXX', chain)
    carbon = fixer.topology.addAtom('C', app.element.carbon, residue)
    if include_oxygen:
        oxygen = fixer.topology.addAtom('O', app.element.oxygen, residue)
        fixer.topology.addBond(carbon, oxygen)
    other = fixer.topology.addResidue('ALA', chain)
    other_carbon = fixer.topology.addAtom('C', app.element.carbon, other)
    fixer.positions = [mm.Vec3(0, 0, 0)]*fixer.topology.getNumAtoms()*unit.nanometer
    fixer.templates = {}
    fixer._standardTemplates = {'ALA'}
    fixer.platform = None

    atoms = [('C', 'C', False), ('O', 'O', True), ('HL', 'H', True),
             ('HO', 'H', True), ('HN', 'H', False)]
    bonds = [('C', 'O'), ('C', 'HL'), ('O', 'HO'), ('C', 'HN')]
    if reverse_definition_bonds:
        bonds = [(b, a) for a, b in bonds]
    definitions = {'XXX': (atoms, bonds)} if source == 'downloaded' else {}
    if source == 'registered':
        template = app.Topology()
        template_residue = template.addResidue('XXX', template.addChain())
        template_atoms = {name: template.addAtom(name, app.element.get_by_symbol(symbol), template_residue)
                          for name, symbol, leaving in atoms}
        for name1, name2 in bonds:
            template.addBond(template_atoms[name1], template_atoms[name2])
        positions = [mm.Vec3(0, 0, 0)]*template.getNumAtoms()*unit.nanometer
        fixer.registerTemplate(template, positions, [leaving for name, symbol, leaving in atoms])
    fixer._downloadNonstandardDefinitions = lambda: definitions
    return fixer, residue, carbon, other_carbon, definitions


def _capture_variants(monkeypatch, fixer):
    captured = []

    def add_hydrogens(self, **kwargs):
        captured.extend(kwargs['variants'])

    monkeypatch.setattr(app.Modeller, 'addHydrogens', add_hydrogens)
    fixer.addMissingHydrogens()
    return captured


@pytest.mark.parametrize('source', ['downloaded', 'registered'])
@pytest.mark.parametrize('include_oxygen', [False, True])
@pytest.mark.parametrize('reverse_external_bond', [False, True])
@pytest.mark.parametrize('reverse_definition_bonds', [False, True])
def test_leaving_hydrogen_variants(monkeypatch, source, include_oxygen,
                                  reverse_external_bond, reverse_definition_bonds):
    """External bonds omit leaving H, but internal bonds and ordinary H remain."""
    fixer, residue, carbon, other, definitions = _make_fixer(source, include_oxygen, reverse_definition_bonds)
    fixer.topology.addBond(other, carbon) if reverse_external_bond else fixer.topology.addBond(carbon, other)
    expected = ([('HO', 'O')] if include_oxygen else []) + [('HN', 'C')]
    assert fixer._describeVariant(residue, definitions) == expected
    assert _capture_variants(monkeypatch, fixer) == [expected, None]


@pytest.mark.parametrize('source', ['downloaded', 'registered'])
def test_variant_bonds_follow_topology_changes(monkeypatch, source):
    """Reusing a fixer after a topology edit must not reuse an old bond index."""
    fixer, residue, carbon, other, definitions = _make_fixer(source)
    all_hydrogens = [('HL', 'C'), ('HO', 'O'), ('HN', 'C')]
    assert _capture_variants(monkeypatch, fixer) == [all_hydrogens, None]
    fixer.topology.addBond(carbon, other)
    assert _capture_variants(monkeypatch, fixer) == [[('HO', 'O'), ('HN', 'C')], None]
    replacement, _, _, _, _ = _make_fixer(source)
    fixer.topology = replacement.topology
    fixer.positions = replacement.positions
    assert _capture_variants(monkeypatch, fixer) == [all_hydrogens, None]


def test_variant_preparation_scans_bonds_once(monkeypatch):
    """The amount of bond traversal stays linear as custom residues are added."""
    fixer, _, _, _, _ = _make_fixer('downloaded', include_oxygen=False)
    chain = next(fixer.topology.chains())
    for _ in range(20):
        residue = fixer.topology.addResidue('XXX', chain)
        fixer.topology.addAtom('C', app.element.carbon, residue)
        fixer.positions.append(mm.Vec3(0, 0, 0)*unit.nanometer)
    bonds = fixer.topology.bonds
    traversals = []

    def counted_bonds():
        traversals.append(True)
        return bonds()

    monkeypatch.setattr(fixer.topology, 'bonds', counted_bonds)
    variants = _capture_variants(monkeypatch, fixer)
    assert variants == [[('HL', 'C'), ('HN', 'C')], None] + [[('HL', 'C'), ('HN', 'C')]]*20
    assert len(traversals) == 1


def test_standard_variants_do_not_scan_bonds(monkeypatch):
    fixer, residue, _, _, _ = _make_fixer('downloaded')
    residue.name = 'ALA'
    fixer._downloadNonstandardDefinitions = lambda: {}

    def unexpected_bond_scan():
        pytest.fail('Standard variants do not require external bond information')

    monkeypatch.setattr(fixer.topology, 'bonds', unexpected_bond_scan)
    assert _capture_variants(monkeypatch, fixer) == [None, None]
