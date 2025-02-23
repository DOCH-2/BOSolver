Proper Installation Test
------------------------

- To check whether the package is normally installed, run the following command. (this requires *pytest* Python package)

```bash
pytest proper_install_test.py
```

- If pytest is not installed, you can run the code directly by running the following command.

```bash
python proper_install_test.py
```

The output would look like as follows:

```bash
$ python proper_install_test.py
>>> Test test_1 Passed
>>> Test test_2 Passed
>>> Test test_3 Passed
>>> Test test_4 Passed
>>> Test test_5 Passed
...
```

- This should run without any errors. If it does not, please check the installation instructions again, or versions of libraries that are installed.

Performance Test
----------------

We prepared the bond perception performance comparison test among BOSolver, RDKit, and Openbabel.
