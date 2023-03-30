import numpy as np
from rmgpy.molecule.molecule import Molecule
from rmgpy.molecule.group import GroupAtom, Group, GroupBond

def functional_group_analysis(smiles, max_n_radius_neighbor=1, max_num_heavy_atoms_in_functional_group=5, min_num_heavy_atoms_in_functional_group=2, max_num_heavy_atoms_in_ring=10):
    functional_group_smiles_set = set()

    molecule = make_rmg_mol(smiles)

    if molecule is None:
        return functional_group_smiles_set
    
    num_heavy_atoms_in_molecule = sum(not atom.is_hydrogen() for atom in molecule.atoms)
    if num_heavy_atoms_in_molecule <= min_num_heavy_atoms_in_functional_group:
        functional_group_smiles_set.add(smiles)
        return functional_group_smiles_set
    
    sssr = molecule.get_smallest_set_of_smallest_rings()
    # monorings, polyrings = molecule.get_disparate_cycles()
    all_rings = []
    for ring in sssr:
        if len(ring) <= max_num_heavy_atoms_in_ring:
            all_rings.append(ring)
    
    all_ring_atoms = set()
    for ring in all_rings:
        all_ring_atoms.update(ring)

    sampled_functional_group_smiles_set = get_ring_functional_groups(molecule, rings=all_rings)
    functional_group_smiles_set.update(sampled_functional_group_smiles_set)
        
    for atom in molecule.atoms:
        if not atom.is_hydrogen():
            sampled_functional_group_smiles = get_n_radius_functional_group(atom, molecule, all_ring_atoms, max_n_radius_neighbor=max_n_radius_neighbor, max_num_heavy_atoms_in_functional_group=max_num_heavy_atoms_in_functional_group, min_num_heavy_atoms_in_functional_group=min_num_heavy_atoms_in_functional_group)
            if sampled_functional_group_smiles is not None:
                functional_group_smiles_set.add(sampled_functional_group_smiles)

    return functional_group_smiles_set

def make_rmg_mol(smiles):
    try:
        molecule = Molecule().from_smiles(smiles)
    except:
        print(f"Could not parse smiles: {smiles}")
        return None
    molecule.sort_atoms()
    return molecule

def get_n_radius_functional_group(center_atom, molecule, all_ring_atoms, max_n_radius_neighbor=None, max_num_heavy_atoms_in_functional_group=5, min_num_heavy_atoms_in_functional_group=3):
    if center_atom in all_ring_atoms:
        return None
    
    num_heavy_atoms_in_molecule = sum(not atom.is_hydrogen() for atom in molecule.atoms)
    if max_n_radius_neighbor is None:
        max_n_radius_neighbor = max(num_heavy_atoms_in_molecule, max_num_heavy_atoms_in_functional_group)

    sampled_functional_group_smiles = None
    for n_radius_neighbor in range(1, max_n_radius_neighbor + 1):
        group = make_group(center_atom, molecule, all_ring_atoms, n_radius_neighbor)
        num_heavy_atoms_in_functional_group = sum(group_atom.atomtype[0].label!="H" for group_atom in group.atoms)
        if min_num_heavy_atoms_in_functional_group <= num_heavy_atoms_in_functional_group and num_heavy_atoms_in_functional_group <= max_num_heavy_atoms_in_functional_group:
            try:
                sampled_mol = group.make_sample_molecule()
            except:
                print(f"Could not make sample molecule from group: {group.to_adjacency_list()}")
            else:
                sampled_mol.sort_atoms()
                sampled_functional_group_smiles = sampled_mol.to_smiles()
        else:
            break
    return sampled_functional_group_smiles

def make_group(center_atom, molecule, all_ring_atoms, n_radius_neighbor=1):

    group_atoms = {}
    group_atoms[center_atom] = GroupAtom(atomtype=[center_atom.element.symbol])
    
    neighbors = list(center_atom.edges.items())
    neighbors.sort()

    group_atoms = get_neighbors(neighbors, group_atoms, all_ring_atoms, n_radius_neighbor, 1)
    
    group = Group(atoms=list(group_atoms.values()))

    group = make_bonds(molecule, group, group_atoms)
    
    group.atoms = group.sort_by_connectivity(group.atoms)
    
    group.update()

    return group


def get_neighbors(atoms, group_atoms, all_ring_atoms, n_radius_neighbor, degree):
    
    for (atom, bond) in atoms:
        if atom not in group_atoms:
            group_atoms[atom] = GroupAtom(atomtype=[atom.element.symbol])

        if atom in all_ring_atoms:
            continue

        if degree + 1 <= n_radius_neighbor:
            neighbors = list(item for item in atom.edges.items() if not item[0].is_hydrogen())
            neighbors.sort()
            group_atoms = get_neighbors(neighbors, group_atoms, all_ring_atoms, n_radius_neighbor, degree + 1)
            
    return group_atoms

def make_bonds(molecule, group, group_atoms):
    
    for atom1, group_atom1 in group_atoms.items():
        for atom2, group_atom2 in group_atoms.items():
            if molecule.has_bond(atom1, atom2):
                bond = molecule.get_bond(atom1, atom2)
                if not group.has_bond(group_atom1,group_atom2):
                    group.add_bond(GroupBond(group_atom1,group_atom2,order=[bond.order]))
            
    return group

def get_ring_functional_groups(molecule, rings=None):
    sampled_functional_group_smiles_set = set()

    if rings is None:
        sssr = molecule.get_smallest_set_of_smallest_rings()
        monorings, polyrings = molecule.get_disparate_cycles()
        rings = sssr + monorings + polyrings
    
    for ring in rings:
        group = make_ring_group(molecule, ring)

        sampled_mol = group.make_sample_molecule()
        sampled_mol.sort_atoms()
        sampled_functional_group_smiles = sampled_mol.to_smiles()
        sampled_functional_group_smiles_set.add(sampled_functional_group_smiles)

    return sampled_functional_group_smiles_set

def make_ring_group(molecule, ring):
    group_atoms = {}
    for atom in ring:
        if atom not in group_atoms:
            group_atoms[atom] = GroupAtom(atomtype=[atom.element.symbol])
            
    group = Group(atoms=list(group_atoms.values()))

    group = make_bonds(molecule, group, group_atoms)

    group.atoms = group.sort_by_connectivity(group.atoms)

    group.update()

    return group

