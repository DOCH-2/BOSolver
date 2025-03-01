# How To Reproduce the Performance Test

1. Prepare data for the test.
Follow the instruction in `../DATA/PubChemQC/HowToPrepare.md` to prepare data.

2. Run the following codes to run performance test.

- Performance test for BOSolver, xyz2mol, Openbabel, and Indigox.

```bash
>>> export SEL=$CWD/test-list/neutral_sel
>>> export EXP=$CWD/whatever_you_want/
>>> export LOG=$CWD/whatever_you_want/{"xyz2mol" or "obabel" or "indigox" or "bosolver"}/LOG.txt
>>> export ERR=$CWD/whatever_you_want/{"xyz2mol" or "obabel" or "indigox" or "bosolver"}/ERR.txt
>>> parallel -a $SEL --timeout 60 \
>>>  python -u ../test.py {} 0 \
>>>  {"--xyz2mol" or "--obabel" or "--indigo" or "--bosolver"} 1>$LOG 2>$ERR
```

Refer to `test-onlocal.sh` for an example.

[!NOTE]
> Performance test for each method can be run asynchronously,
> so if computing resource is available, we recommend you to run all test simultaneously.
> If you are using computing cluster, we suggest you to check out
> `test-batch.sh` and `tester.sh`.

- Preparation for label (SMILES from PubChem)

```bash
>>> cd label
>>> source run_parallel.sh
```

[!NOTE]
> Make sure you have proper symbolic link or copied file of
> "DATA/PubChemQC/name_chg_smi_neutral_sel.woRadicalXXXXXX"


3. Run the following codes to analyze the result.

```bash
>>> cd $CWD/whatever_you_want
>>> python ../../analyze.py
>>> python ../../compare_with_label.py
```

4. To visualize the result, use jupyter notebook file.

```bash
>>> cd ../
>>> jupyter notebook post_analyze.ipynb
```
