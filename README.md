<p align="right">
  <img src="biobloom_logo" alt="BioBloom Logo" width="120">
</p>

# BioBloom 🌿🫀

**Microalgae as a Source of Bioactive Compounds in Pharmacology**

---

## Overview

BioBloom is a research project exploring **microalgae proteomes** as a source of bioactive peptides with pharmacological potential (especially ACE inhibitors for cardiovascular health).  

This repository provides a collection of scripts for:
- querying **BIOPEP-UWM**,
- peptide filtering & ADMET screening,
- **SMILES/FASTA/PDB** conversions,
- peptide/receptor structure generation for docking.

---

## Repository Structure

- `biopep_uwm/` → BIOPEP-UWM automation  
- `admet/` → unify and filter ADMET outputs  
- `molecular_docking/` → peptide & receptor structure prep  
- `pLM4ACE/` → input/output formatting for pLM4ACE model  
- `smiles_conversion/` → utilities for FASTA/SMILES conversion  

---

## Workflow

Proteome FASTA
↓
BIOPEP-UWM automation → Peptide lists
↓
Unify/Deduplicate → (optional) pLM4ACE prep / SMILES↔FASTA
↓
ADMET screening (AdmetLab) → Select top peptides
↓
Peptide 3D (PyRosetta) + Receptor repack
↓
Docking preparation/output

## Documentation

See [https://github.com/martatru/BioBloom/wiki](Wiki) for **full documentation**:
- detailed script descriptions,
- requirements,
- usage examples,
- full pipeline workflow.
