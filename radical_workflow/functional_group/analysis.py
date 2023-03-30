from rmgpy.molecule.molecule import Molecule
from rmgpy.molecule.group import GroupAtom, Group, GroupBond

def functional_group_analysis(smiles, max_num_heavy_atoms_in_functional_group=5):
    functional_group_smiles_set = set()

    molecule = make_rmg_mol(smiles)

    if molecule is None:
        return functional_group_smiles_set
    
    sssr = molecule.get_smallest_set_of_smallest_rings()
    monorings, polyrings = molecule.get_disparate_cycles()
    all_rings = sssr + monorings + polyrings
    
    all_ring_atoms = set()
    for ring in all_rings:
        all_ring_atoms.update(ring)

    sampled_functional_group_smiles_set = get_ring_functional_groups(molecule, rings=all_rings)
    functional_group_smiles_set.update(sampled_functional_group_smiles_set)
    
    for atom in molecule.atoms:
        if not atom.is_hydrogen():
            sampled_functional_group_smiles = get_n_radius_functional_group(atom, molecule, all_ring_atoms, max_num_heavy_atoms_in_functional_group=max_num_heavy_atoms_in_functional_group)
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

def get_n_radius_functional_group(center_atom, molecule, all_ring_atoms, max_num_heavy_atoms_in_functional_group):
    if center_atom in all_ring_atoms:
        return None
    
    sampled_functional_group_smiles = None
    num_heavy_atoms = sum(not atom.is_hydrogen() for atom in molecule.atoms)
    max_n_radius_neighbor = min(max_num_heavy_atoms_in_functional_group, num_heavy_atoms)
    for n_radius_neighbor in range(1, max_n_radius_neighbor + 1):
        group = make_group(center_atom, molecule, all_ring_atoms, n_radius_neighbor)
        num_heavy_atoms_in_functional_group = sum(group_atom.atomtype[0].label!="H" for group_atom in group.atoms)
        if num_heavy_atoms_in_functional_group <= max_num_heavy_atoms_in_functional_group:
            sampled_mol = group.make_sample_molecule()
            sampled_mol.sort_atoms()
            sampled_functional_group_smiles = sampled_mol.to_smiles()
        else:
            break
    return sampled_functional_group_smiles

def make_group(center_atom, molecule, all_ring_atoms, n_radius_neighbor=1):

    group_atoms = {}
    group_atoms[center_atom] = GroupAtom(atomtype=[center_atom.atomtype])
    
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
            group_atoms[atom] = GroupAtom(atomtype=[atom.atomtype])

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

