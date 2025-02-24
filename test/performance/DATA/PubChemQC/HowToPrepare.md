How to prepare data for PubChemQC2025

Required data
- 1-298189.json : a piece of raw data file from PubChemQC project
- name_chg_smi_neutral_XXX : txt file that records names (cid), charges, and smi's of data in json file
- neutral_XXX.(with|wo)RadicalXXXXXX : path to xyz files that will be used for performance test

Procedure how all files here are made are as below.

1. ./1-298189.json
Download from HuggingFace Hub.
[1-298189.json](https://huggingface.co/datasets/molssiai-hub/pubchemqc-pm6/blob/main/data/pm6opt_chnopsfclnakmgca500/train/000000001-000298189.json)

2. xyz files in ./neutral, ./anion (empty), and ./cation (empty)
run `python ./json2xyz.py`
This will parse json file into separate xyz file. We do this for parallelization.
[!CAUTION]
> This creates ~*2M* xyz files. Be cautious of I/O overload.

3. name_chg_smi_neutral
run `bash name_chg_smi_neutral.sh`
This will create the label file (which contains their names(cid), total charges, and PubChem SMILES)

4. neutral_tot.withRadicalXXXXXX
run `find neutral -type f > neutral_tot.withRadicalXXXXXX`
This simply records paths to xyz files.

5. neutral_tot.woRadicalXXXXXX
run `python ./radical_check.py > neutral_tot.woRadicalXXXXXX`
This will select xyz files of systems with even number of electrons only.

6. neutral_sel.woRadicalXXXXXX & name_chg_smi_neutral
run `python ./data_selector.py`
This will randomly choose 1M xyz files for the test.
