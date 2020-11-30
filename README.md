# README

This repository contains the code used to obtain some of our results in the paper "Efficient MILP modelings for Sboxes and Linear Layers of SPN ciphers."

- `sbox` contains python (or SageMath) files that are useful for modeling 4-bit Sboxes.
- `impossible_differentials` contains python files that build Gurobi models for searching for impossible differentials for Skinny and the AES.

We have not published the code for 8-bit Sboxes yet but the models of the DDTs used in impossible differential search are given as pickle files.

## Overview of the `sbox` directory

- `arbitrary_sbox_gen.py` generates a pickle file for an arbitrary Sbox DDT (only transitions zero -> non zero are impossible).
- `check_model.sage` checks the correctness of a model of a DDT given by a pickle file.
- `convex_hull.sage` generates a big set of inequalities with the convex hull technique. Uses the SageMath Sboxes and Polyhedra tools.
- `identity_sbox_gen.py` generates a pickle file for the identity Sbox DDT.
- `minimize.py` performs step 2 given a big set of inequalities with greedy or minimization techniques.
- `utilities.py` defines some small useful functions for the other files. 

## Overview of the `impossible_differentials` directory

- `aes_equiv_sbox.pkl` is the model of the DDT of an affine equivalent AES Sbox.
- `aes.py` builds and tests the Gurobi model for the AES.
- `arbitrary_sbox_8_8.pkl` is the model of the DDT of an arbitrary 8-bit Sbox (for testing purposes).
- `identity_sbox_8.pkl` is the model of the DDT of the identity 8-bit Sbox (testing).
- `main_aes.py` launches the search for impossible differentials for the AES.
- `main_skinny.py` is the same for Skinny.
- `primitive.py` contains classes and functions common to `aes.py` and `skinny.py`.
- `skinny.py` builds and tests the Gurobi model for Skinny.
- `skinny_sbox.pkl` is the model of the DDT of the Skinny 8-bit Sbox.
- `utilities.py` defines small useful functions.

