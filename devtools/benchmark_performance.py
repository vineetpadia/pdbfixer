"""Offline performance benchmarks. Run from the repository root with PYTHONPATH=.

Each case is warmed once before collecting samples. Synthetic cases isolate
scaling costs; they do not represent complete molecular repair timings.
"""

import argparse
import gc
import json
import os
from pathlib import Path
import platform
import statistics
from time import perf_counter
from unittest.mock import patch

os.environ.setdefault('OPENMM_CPU_THREADS', '1')

import numpy as np
import openmm as mm
from openmm import app, unit
import pdbfixer
from pdbfixer.pdbfixer import PDBFixer, Sequence


DATA = Path(pdbfixer.__file__).parent / 'tests' / 'data'


def sequence_case(size, insertion_codes=False):
    fixer = PDBFixer.__new__(PDBFixer)
    fixer.topology = app.Topology()
    chain = fixer.topology.addChain('A')
    if insertion_codes:
        for i in range(size):
            fixer.topology.addResidue('GLY', chain, str(i//2+1), insertionCode='A' if i%2 else '')
        fixer.sequences = [Sequence('A', ['GLY']*size)]
    else:
        fixer.topology.addResidue('GLY', chain, '1')
        fixer.sequences = [Sequence('A', ['ALA']*size)]

    def run():
        fixer.findMissingResidues()
        assert fixer.missingResidues == {}

    return run


def variant_case(size):
    fixer = PDBFixer.__new__(PDBFixer)
    fixer.topology = app.Topology()
    chain = fixer.topology.addChain('A')
    previous = None
    for _ in range(size):
        residue = fixer.topology.addResidue('XXX', chain)
        atom = fixer.topology.addAtom('C', app.element.carbon, residue)
        if previous is not None:
            fixer.topology.addBond(previous, atom)
        previous = atom
    fixer.positions = [mm.Vec3(i*0.15, 0, 0) for i in range(size)]*unit.nanometer
    fixer.templates = {}
    fixer._standardTemplates = set()
    fixer.platform = None
    definitions = {'XXX': ([('C', 'C', False), ('HL', 'H', True), ('HN', 'H', False)],
                           [('C', 'HL'), ('C', 'HN')])}
    fixer._downloadNonstandardDefinitions = lambda: definitions

    def check_variants(modeller, **kwargs):
        assert kwargs['variants'] == [[('HN', 'C')]]*size

    def run():
        # Isolate PDBFixer's variant preparation from OpenMM's minimization.
        with patch.object(app.Modeller, 'addHydrogens', check_variants):
            fixer.addMissingHydrogens()

    return run


def topology_case(missing_residues):
    fixer = PDBFixer(filename=str(DATA / '4JSV.pdb'))
    fixer.findMissingResidues()
    if not missing_residues:
        fixer.missingResidues = {}
    fixer.findMissingAtoms()

    def run():
        topology, positions, new_atoms, _ = fixer._addAtomsToTopology(True, True)
        assert topology.getNumAtoms() == len(positions)
        assert len(new_atoms) == (1794 if missing_residues else 4)

    return run


def complete_structure_case():
    fixer = PDBFixer(filename=str(DATA / 'test.pdb'))
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    return fixer.addMissingAtoms


def repair():
    fixer = PDBFixer(filename=str(DATA / '1BHL.pdb'), platform=mm.Platform.getPlatformByName('CPU'))
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms(seed=1)
    assert fixer.topology.getNumAtoms() == 1225
    assert np.isfinite(fixer.positions.value_in_unit(unit.nanometer)).all()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repeat', type=int, default=5)
    parser.add_argument('--include-repair', action='store_true', help='Also run the slower CPU loop-repair pipeline.')
    args = parser.parse_args()
    if args.repeat < 1:
        parser.error('--repeat must be positive')
    cases = {
        'unmatched_sequence_32000': sequence_case(32000),
        'insertion_codes_8000': sequence_case(8000, insertion_codes=True),
        'hydrogen_variants_2000': variant_case(2000),
        '4JSV_topology_no_missing_residues': topology_case(False),
        '4JSV_topology_with_missing_residues': topology_case(True),
        'complete_structure_add_missing_atoms': complete_structure_case(),
    }
    if args.include_repair:
        cases['1BHL_cpu_repair_pipeline'] = repair
    results = {}
    for name, run in cases.items():
        run()
        samples = []
        for _ in range(args.repeat):
            # Topology objects contain cycles.  Collect outside the timed region
            # so cleanup from a previous sample does not distort the next one.
            gc.collect()
            gc.disable()
            try:
                start = perf_counter()
                run()
                samples.append(perf_counter()-start)
            finally:
                gc.enable()
        results[name] = {'median_seconds': statistics.median(samples), 'samples_seconds': samples}
    print(json.dumps({
        'source': str(Path(pdbfixer.__file__).resolve()),
        'python': platform.python_version(),
        'openmm': mm.__version__,
        'numpy': np.__version__,
        'cpu_threads': os.environ['OPENMM_CPU_THREADS'],
        'repeat': args.repeat,
        'results': results,
    }, indent=2))


if __name__ == '__main__':
    main()
