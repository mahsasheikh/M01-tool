
---

# **M01 tool User Guide**

## **Introduction**

Welcome to **M01 tool**, a fully automated tool designed for generating small molecule-peptide hybrids and docking them into curated protein structures. This user manual will guide you through the installation process, hybrid generation, and docking simulation steps, ensuring smooth usage of the software.

---

## **Installation**

### **Pre-requisites**

Ensure that your system has the required dependencies installed. If you encounter issues during installation, we recommend checking out our Jupyter notebook for further assistance.

### **Manual Installation Steps**

You can manually install the required packages via `conda` and `pip`. Use the following commands to set up your environment:

#### **1. Install RDKit and PDBFixer via Conda**
```bash
conda install -c conda-forge rdkit pdbfixer
```

#### **2. Install PyTorch and Dependencies for CPU**
```bash
pip install torch==1.13.1+cpu -f https://download.pytorch.org/whl/cpu/torch_stable.html
```

#### **3. Install PyTorch Geometric Dependencies**
```bash
pip install torch-scatter -f https://data.pyg.org/whl/torch-1.13.1+cpu.html
pip install torch-sparse -f https://data.pyg.org/whl/torch-1.13.1+cpu.html
pip install torch-spline-conv torch-geometric==2.0.1 -f https://data.pyg.org/whl/torch-1.13.1+cpu.html
```

#### **4. Install Additional Tools**
```bash
pip install cairosvg svgutils molvs
pip install rdkit-pypi
pip install Bio
pip install meeko
pip install paramiko
pip install easydock
pip install -q git+https://github.com/mayrf/pkasolver.git
pip install vina
```

#### **5. Install MGLTools**
```bash
conda install bioconda::mgltools
```

---

## **Hybrid Generation**

In this section, you'll learn how to generate hybrid molecules using M01 tool. Follow the steps below to create your hybrids.

### **Steps for Hybrid Generation:**

1. **Enter the Main SMILES:**
   - Input the SMILES structure of the main molecule for the next steps.

2. **Select Connection Points and Amino Acids:**
   - Choose the connection point on your molecule.
   - Enter the connection point, amino acids, and any other amino acid-like molecules.
   - **Note:** Ensure that the "other molecule" can form a peptide bond to avoid errors or break esters in your structure.
   - Choose the N or C terminal of the peptides for the connection point, as it might cause errors or break esters if present in your structure.

3. **Molecule and Peptide Considerations:**
   - Ensure the connection site and the terminal make chemical sense. For example, if you are connecting the nitrogen (N) in your original molecule to the carboxyl group (COOH) in the amino acid, you should select the "C terminal" as in the C in the COOH of the peptide.
   - Limit your input to a total of **7 molecules** to prevent long operation times and avoid generating ligands too large for docking. **Note:** This module is not designed for docking large peptides.

---


## **3D Structure Retrieval**

### **Introduction**

In this section, we download the 3D structure of the receptor protein and make it ready for docking. 
 The first step in docking is to have a structure of a given target protein. While in some cases a high-quality comparative model will be used, most cases start with an experimentally (X-ray, NMR, cryoEM) solved three-dimensional structure.
In these cases, a given target protein structure can be downloaded with BioPython using a given accession ID. We can directly download this structure in `.pdb` file format.
In addition, we calculate the centroid of the chain in this step, in order to use in further calculations.


### **Steps for Getting the 3D Structure:**

1. **Enter the PDB ID and Chain or UniProt ID:**
   - You can either use the UniProt API to select the appropriate PDB file and chain, or input them directly into M01 tool.

2. **Run the Retrieval Process:**
   - After providing the required IDs, the system will download the 3D structure in `.pdb` format.
   - Additionally, the centroid of the chain is calculated to serve as the center for grid box definition in docking simulations.

---

## **Automated Docking**

### **Introduction**

In this section, we will walk you through the automated docking process, which involves preparing both the ligand and the receptor, and finally performing the docking simulation. The docking simulation is conducted using **Autodock Vina**, while ligand preparation is handled by **EasyDock**.

### **Steps for Automated Docking:**

1. **Preparing the Ligand with EasyDock:**
   - The ligand must be prepared for docking. You can choose to run an optional protonation step or skip it.
   - The **EasyDock** module is used to prepare the ligand, and further details are available on their [GitHub page](https://github.com/ci-lab-cz/easydock).

2. **Preparing the Receptor with MGL Tools:**
   - The receptor's PDB file is converted into PDBQT format. Hydrogens are checked, and **Gasteiger charges** are added.
   - For additional options in receptor preparation, please refer to the `prepare_receptor4.py` file from **MGLTools**.

3. **Performing the Molecular Docking Simulation:**
   - **Autodock Vina** is used to perform the docking simulation. The runtime of each simulation varies depending on the complexity of the ligand.
   - We use default settings for the grid box and docking configuration. For customization, please refer to [Vina’s documentation](https://autodock-vina.readthedocs.io/_/downloads/en/latest/pdf/).

---

## **Calculating Molecular Descriptors**

### **Introduction**

M01 tool offers the ability to calculate various molecular descriptors using the **RDKit** module. These descriptors provide insights into the chemical properties of molecules and their potential drug-likeness. Below is a list of the descriptors calculated by M01 tool, along with simple explanations of each:

### **Available Descriptors:**

1. **Lipinski's Rule of Five (lip_score):**
   - **Purpose:** Evaluates the drug-likeness of a compound based on four key factors:
     - **Molecular Weight (exact_mw):** ≤ 500 daltons
     - **LogP (logP):** ≤ 5
     - **Hydrogen Bond Donors (hbd):** ≤ 5
     - **Hydrogen Bond Acceptors (hba):** ≤ 10
   - A score is added for each violation of these rules, and the total number of violations is calculated.

2. **Exact Molecular Weight (exact_mw):**
   - The precise mass of a molecule, calculated based on the atomic masses of its constituent atoms.

3. **LogP (logP):**
   - The logarithm of the partition coefficient (P) between n-octanol and water. This measures the compound's hydrophobicity and its ability to cross cell membranes.

4. **Hydrogen Bond Donors (hbd):**
   - The count of hydrogen atoms that can participate in hydrogen bonding by donating a hydrogen atom.

5. **Hydrogen Bond Acceptors (hba):**
   - The count of atoms capable of accepting a hydrogen
  

---
**LigandBuilder Class**:

---

## **Functions in the LigandBuilder Class**

### **1. create_peptide_bond**
   - **Description:** Uses RDKit reactions to form a peptide bond by reacting a COOH group with a nitrogen.
   - **Key Point:** COOH groups are prioritized over COO to prevent ester breakage during the reaction.

### **2. add_hydrogens_and_sanitize**
   - **Description:** Ensures the chemical correctness of the molecule and adds all necessary hydrogen atoms.
   - **Usage:** This function is essential for preparing molecules for further steps, ensuring the structure's validity.

### **3. remove_extra_hydrogens**
   - **Description:** Removes any redundant hydrogen atoms to avoid errors in peptide bond generation.
   - **Usage:** This helps ensure that the structure is ready for bond formation without overhydration issues.

### **4. handle_N_heterocycles**
   - **Description:** Manages nitrogen-containing aromatic groups. If a **kekulization** error occurs, it likely indicates an unsupported heterocycle.
   - **Note:** If you encounter a kekulization error, contact us via email for troubleshooting assistance.

### **5. Input Amino Acids and Molecules**
   - **Description:** 
     - **input_amino_acids:** Use this function to input amino acids.
     - **input_molecule:** Allows for the inclusion of other molecules for generating peptide combinations.
   - **Usage:** These functions help in the generation of diverse peptides and hybrid molecules.

### **6. Find Connection Points**
   - **Description:** Identifies all NH groups within the peptide for N-terminal connections to COOH groups.
   - **Steps:**
     - **AddCs Function:** Adds a temporary atom (e.g., cesium) to the N-terminal. This atom will later be replaced by the target structure.
     - **Customization:** To replace a larger substructure, modify the function by including your desired **SMARTS** code in the designated section.
  
---

**Get_PDB Class**:

---

## **Functions in the Get_PDB Class**

### **1. SelectChains**
   - **Description:** Helper class used to accept or reject chains based on their chain IDs.
   - **Usage:** Initialized with a list of chain IDs. The `accept_chain` method checks if a given chain's ID matches the accepted list of IDs.

### **2. fetch_uniprot_data**
   - **Description:** Fetches UniProt data for a given UniProt ID.
   - **Details:** Sends a GET request to the UniProt database and returns XML data if successful.

### **3. parse_pdb_entries**
   - **Description:** Parses PDB entries from UniProt XML data.
   - **Details:** Extracts PDB IDs, methods (e.g., X-ray, NMR), resolutions, and chain details from the XML data. Returns a list of tuples containing this information.

### **4. parse_chain_details**
   - **Description:** Extracts chain details from UniProt XML data.
   - **Details:** Retrieves chain IDs and the range of amino acids for each chain. Also calculates the number of amino acids in each chain.

### **5. select_best_pdb**
   - **Description:** Selects the best PDB entry based on chain length and resolution.
   - **Details:** Sorts the PDB entries by chain length (in descending order) and resolution (in ascending order), returning the best entry for further use.

### **6. calculate_centroid**
   - **Description:** Calculates the centroid of a PDB structure.
   - **Details:** Computes the mean coordinates of all atoms in the structure to determine its centroid, used for docking and further calculations.

### **7. make_pdb**
   - **Description:** Extracts a specific chain from a PDB file and saves it to a temporary directory.
   - **Details:** 
     - Uses `PDB.PDBParser` and `PDB.PDBIO` from Biopython to parse and write PDB structures.
     - Adds a remark with the centroid coordinates to the saved PDB file.

### **8. pdbfixer**
   - **Description:** Uses PDBFixer to correct missing atoms and residues in a PDB file.
   - **Details:** 
     - Finds and adds missing residues and atoms.
     - Replaces nonstandard residues and saves the fixed PDB file, with the centroid information included as a remark.

### **9. get_pdb_from_uniprot**
   - **Description:** Processes a UniProt ID to fetch PDB data and select the best PDB entry.
   - **Details:** 
     - Retrieves the PDB file, extracts the specified chain, fixes the PDB structure, and returns the PDB ID and chain.

### **10. get_pdb_from_PDB**
   - **Description:** Processes a specified PDB ID and chain directly.
   - **Details:** Retrieves the PDB file, extracts the specified chain, fixes the PDB structure, and returns the PDB ID and chain.

---
**prepare_receptor4.py**:

---

## **Usage:**
```bash
prepare_receptor4.py -r receptor_filename
```

## **Description:**
- **`-r receptor_filename`**: Specifies the receptor file to be processed. Supported file types include:
  - `.pdb`, `.mol2`, `.pdbq`, `.pdbqs`, `.pdbqt`, `.pqr`, `.cif`.

## **Optional Parameters:**

- **`[-v]`**: Enables verbose output. The default is minimal output.
  
- **`[-o pdbqt_filename]`**: Specifies the output filename for the `.pdbqt` file (default is `molecule_name.pdbqt`).

- **`[-A]`**: Specifies the types of repairs to make:
  - `'bonds_hydrogens'`: Build bonds and add hydrogens.
  - `'bonds'`: Build a single bond from each atom with no bonds to its closest neighbor.
  - `'hydrogens'`: Add hydrogens.
  - `'checkhydrogens'`: Add hydrogens only if there are none already.
  - `'None'`: Do not make any repairs (default).

- **`[-C]`**: Preserve all input charges (default adds Gasteiger charges).

- **`[-p atom_type]`**: Preserve input charges on specific atom types, e.g., `-p Zn`, `-p Fe`.

- **`[-U]`**: Specifies cleanup actions:
  - `'nphs'`: Merge charges and remove non-polar hydrogens.
  - `'lps'`: Merge charges and remove lone pairs.
  - `'waters'`: Remove water residues.
  - `'nonstdres'`: Remove chains composed entirely of non-standard residues.
  - `'deleteAltB'`: Remove alternate location atoms labeled `@B`, and rename `@A` atoms to their original form (default is `'nphs_lps_waters_nonstdres'`).

- **`[-e]`**: Deletes every non-standard residue from any chain. A residue will be deleted if its name is not in this list:
  ```bash
  ['CYS', 'ILE', 'SER', 'VAL', 'GLN', 'LYS', 'ASN', 'PRO', 'THR', 'PHE', 'ALA', 'HIS', 'GLY', 'ASP', 'LEU', 'ARG', 'TRP', 'GLU', 'TYR', 'MET', 'HID', 'HSP', 'HIE', 'HIP', 'CYX', 'CSS']
  ```
  Note: This list includes only standard amino acids. Nucleic acids and metals are not included.  
  Default: `False` (meaning no residues are deleted).

- **`[-M]`**: Enables interactive mode (default is automatic with no further user input required).

- **`[-d dictionary_filename]`**: Saves receptor summary information to the specified file.

- **`[-w]`**: Assigns each receptor atom a unique name, where the new name is the original name plus its index (1-based).

---
