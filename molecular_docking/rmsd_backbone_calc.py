import os
import numpy as np
from Bio.PDB import PDBParser, PDBIO
import warnings
warnings.filterwarnings('ignore')

class RMSDCalculator:
    def __init__(self, reference_pdb, ref_ligand_chain='A', ref_ligand_resname='ACE'):
        """
        Initialize calculator with reference structure (1O86)
        
        Args:
            reference_pdb: Path to reference PDB file (1O86)
            ref_ligand_chain: Chain ID of ligand (lisinopril)
            ref_ligand_resname: Residue name of ligand (ACE, ACY, etc.)
        """
        self.parser = PDBParser(QUIET=True)
        self.reference_pdb = reference_pdb
        self.ref_ligand_chain = ref_ligand_chain
        self.ref_ligand_resname = ref_ligand_resname
        
        # Load reference structure
        self.ref_structure = self.parser.get_structure('reference', reference_pdb)
        self.ref_ligand_coords = self._extract_ligand_backbone(self.ref_structure)
        
        print(f"Reference structure loaded: {len(self.ref_ligand_coords)} backbone atoms from ligand")
    
    def _extract_ligand_backbone(self, structure, chain_id='L', ligand_resname=None):
        """
        Extract backbone-like atoms from ligand (LPR/lisinopril)
        For lisinopril, backbone = main chain carbons: C1, C2, C3, C4 + their attached N, O
        These are analogous to CA, N, C=O in peptides
        
        Args:
            structure: Bio.PDB Structure object
            chain_id: Chain ID containing the ligand
            ligand_resname: Residue name of ligand (if None, use self.ref_ligand_resname)
        
        Returns:
            coords: Nx3 array of atomic coordinates
        """
        if ligand_resname is None:
            ligand_resname = self.ref_ligand_resname
        
        # Backbone atoms for lisinopril: main chain C + attached N,O
        # C1-C2-C3-C4 is the main chain, N1/N2/N3 and O1/O2/O3 are attached
        backbone_atoms = ['C1', 'N1', 'O1', 'C2', 'N2', 'O2', 'C3', 'N3', 'O3', 'C4']
        coords = []
        
        for model in structure:
            if chain_id in model:
                chain = model[chain_id]
                for residue in chain:
                    if residue.resname == ligand_resname:
                        for atom_name in backbone_atoms:
                            if atom_name in residue:
                                atom = residue[atom_name]
                                coords.append(atom.coord)
        
        return np.array(coords)
    
    def _extract_peptide_backbone(self, structure, chain_id='B'):
        """
        Extract backbone atoms (CA, N, C) from peptide chain
        
        Args:
            structure: Bio.PDB Structure object
            chain_id: Chain ID of peptide
        
        Returns:
            coords: Nx3 array of atomic coordinates
            res_info: list of residue information for debugging
        """
        coords = []
        res_info = []
        
        for model in structure:
            if chain_id in model:
                chain = model[chain_id]
                for residue in chain:
                    res_id = residue.id[1]
                    resname = residue.resname
                    
                    # Extract backbone atoms from peptide
                    for atom_name in ['CA', 'N', 'C']:
                        if atom_name in residue:
                            atom = residue[atom_name]
                            coords.append(atom.coord)
                            res_info.append((res_id, resname, atom_name))
        
        return np.array(coords), res_info
    
    def kabsch_alignment(self, mobile_coords, reference_coords):
        """
        Align mobile coordinates to reference using Kabsch algorithm
        
        Returns:
            rmsd: Root mean square deviation
            rotated_coords: Aligned coordinates
        """
        if len(mobile_coords) == 0 or len(reference_coords) == 0:
            return None, None
        
        # Use minimum available atoms for alignment
        n_atoms = min(len(mobile_coords), len(reference_coords))
        mobile = mobile_coords[:n_atoms]
        reference = reference_coords[:n_atoms]
        
        # Center coordinates
        mobile_centered = mobile - mobile.mean(axis=0)
        reference_centered = reference - reference.mean(axis=0)
        
        # SVD for optimal rotation
        H = mobile_centered.T @ reference_centered
        U, S, Vt = np.linalg.svd(H)
        R = Vt.T @ U.T
        
        # Correct rotation matrix (ensure det = 1)
        if np.linalg.det(R) < 0:
            Vt[-1, :] *= -1
            R = Vt.T @ U.T
        
        # Apply rotation
        rotated = mobile_centered @ R.T
        
        # Calculate RMSD
        rmsd = np.sqrt(np.mean(np.sum((rotated - reference_centered) ** 2, axis=1)))
        
        return rmsd, n_atoms
    
    def calculate_rmsd(self, peptide_pdb, pep_chain_id='B', lig_chain_id='A'):
        """
        Calculate RMSD between peptide backbone and reference ligand
        
        Args:
            peptide_pdb: Path to docking output PDB
            pep_chain_id: Chain ID of peptide in docking result
            lig_chain_id: Chain ID of ligand in docking result (unused, for clarity)
        
        Returns:
            rmsd_value: RMSD in Angstroms
            n_atoms: Number of atoms used in alignment
        """
        try:
            # Load docking result structure
            peptide_structure = self.parser.get_structure('peptide', peptide_pdb)
            
            # Extract backbone from peptide chain
            pep_coords, res_info = self._extract_peptide_backbone(peptide_structure, pep_chain_id)
            
            if len(pep_coords) == 0:
                return None, 0
            
            # Calculate RMSD against reference ligand backbone
            rmsd, n_atoms = self.kabsch_alignment(pep_coords, self.ref_ligand_coords)
            
            return rmsd, n_atoms
        
        except Exception as e:
            print(f"  Error processing {peptide_pdb}: {str(e)}")
            return None, 0
    
    def batch_calculate(self, peptide_dir, output_file='rmsd_results.txt', 
                       pep_chain_id='B'):
        """
        Calculate RMSD for all docking output PDB files
        
        Args:
            peptide_dir: Directory containing docking output PDB files
            output_file: Output file for results
            pep_chain_id: Chain ID of peptide in docking outputs
        """
        results = []
        pdb_files = sorted([f for f in os.listdir(peptide_dir) 
                           if f.lower().endswith('.pdb')])
        
        print(f"\nProcessing {len(pdb_files)} docking output structures...\n")
        
        for i, pdb_file in enumerate(pdb_files, 1):
            pdb_path = os.path.join(peptide_dir, pdb_file)
            
            rmsd, n_atoms = self.calculate_rmsd(pdb_path, pep_chain_id=pep_chain_id)
            
            if rmsd is not None:
                results.append({
                    'peptide': pdb_file,
                    'rmsd': rmsd,
                    'n_atoms': n_atoms
                })
                status = f"✓ RMSD: {rmsd:.3f} Å ({n_atoms} atoms)"
            else:
                status = f"✗ No backbone atoms found in chain '{pep_chain_id}'"
            
            print(f"[{i:2d}/{len(pdb_files)}] {pdb_file:35s} {status}")
        
        # Write results to file
        if results:
            with open(output_file, 'w') as f:
                f.write("Peptide\tRMSD (Å)\tAtoms Used\n")
                for r in sorted(results, key=lambda x: x['rmsd']):
                    f.write(f"{r['peptide']}\t{r['rmsd']:.4f}\t{r['n_atoms']}\n")
            
            print(f"\n✓ Results saved to: {output_file}")
            
            # Summary statistics
            rmsd_values = [r['rmsd'] for r in results]
            print(f"\n{'='*50}")
            print(f"Summary statistics:")
            print(f"  Mean RMSD:   {np.mean(rmsd_values):.3f} Å")
            print(f"  Median RMSD: {np.median(rmsd_values):.3f} Å")
            print(f"  Min RMSD:    {np.min(rmsd_values):.3f} Å (best)")
            print(f"  Max RMSD:    {np.max(rmsd_values):.3f} Å (worst)")
            print(f"{'='*50}")
        else:
            print("\n✗ No valid results generated")
        
        return results

# Usage
if __name__ == "__main__":
    import urllib.request
    
    ref_pdb = "1O86.pdb"
    
    # Download reference if needed
    if not os.path.exists(ref_pdb):
        print("Downloading reference structure 1O86...")
        try:
            urllib.request.urlretrieve(
                f"https://files.rcsb.org/download/{ref_pdb}",
                ref_pdb
            )
            print(f"✓ Downloaded {ref_pdb}\n")
        except:
            print(f"✗ Could not download {ref_pdb}\n")
    
    # DIAGNOSTIC: Check reference structure
    print("="*60)
    print("DIAGNOSTIC: Checking reference structure 1O86")
    print("="*60)
    ref_structure = PDBParser(QUIET=True).get_structure('ref', ref_pdb)
    
    for model in ref_structure:
        chains = list(model)
        chain_ids = [c.id for c in chains]
        print(f"\nAvailable chains in 1O86: {chain_ids}")
        
        for chain in chains:
            chain_id = chain.id
            residues = list(chain)
            print(f"\n  Chain '{chain_id}': {len(residues)} residues")
            for res in residues:
                atoms = [a.name for a in res]
                print(f"    Res {res.id[1]} ({res.resname}): {len(atoms)} atoms {atoms[:8]}")
                if len(residues) > 5:
                    break  # Show only first residue per chain for brevity
    
    print("\n" + "="*60)
    print("Lisinopril (LPR) atoms - choose which correspond to backbone:")
    print("="*60)
    for model in ref_structure:
        if 'L' in [c.id for c in model]:
            chain_l = model['L']
            for res in chain_l:
                if res.resname == 'LPR':
                    all_atoms = [(a.name, a.coord) for a in res]
                    print(f"\nAll {len(all_atoms)} atoms in LPR:")
                    for i, (atom_name, coord) in enumerate(all_atoms):
                        print(f"  {atom_name:4s} -> coordinates: {coord}")
    
    print("\n" + "="*60)
    print("Choose atoms from LPR that match backbone concept")
    print("(typically C atoms on main chain)")
    print("="*60 + "\n")
    
    # Initialize calculator
    calculator = RMSDCalculator(ref_pdb, ref_ligand_chain='L', ref_ligand_resname='LPR')
    
    # DIAGNOSTIC: Check structure of first peptide file
    print("\n" + "="*60)
    print("DIAGNOSTIC: Checking structure of first peptide file")
    print("="*60)
    peptide_directory = "./peptides"
    first_file = sorted([f for f in os.listdir(peptide_directory) 
                        if f.lower().endswith('.pdb')])[0]
    
    test_structure = PDBParser(QUIET=True).get_structure('test', 
                                                          os.path.join(peptide_directory, first_file))
    
    for model in test_structure:
        chains = list(model)
        chain_ids = [c.id for c in chains]
        print(f"\nAvailable chains: {chain_ids}")
        
        for chain in chains:
            chain_id = chain.id
            residues = list(chain)
            print(f"\n  Chain '{chain_id}': {len(residues)} residues")
            for res in residues[:3]:  # Show first 3 residues
                atoms = [a.name for a in res]
                print(f"    Res {res.id[1]} ({res.resname}): atoms {atoms[:10]}")
    
    print("\n" + "="*60)
    print("Update parameters above based on actual structure")
    print("="*60 + "\n")
    
    # Process peptides
    results = calculator.batch_calculate(peptide_directory, pep_chain_id='P')