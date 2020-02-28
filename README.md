# Folder contains codes for federated, over the air, PCA


Here are the list of files:

A: For the fully observed, static subspace setting
1. taubatch_test.m: this file considers the case when we are allowed to vary $\tau$ -- the number of iterations after which we normalize the power method output.
2. sigma_c_tes.m: this file has the codes for varying the channel noise seen at each iteration.
3. ratio_test.m: this file studies the variation in the eigen ratio.


B: For the STmiss problem which may contain missing entries
4. NORST_fed.m: this script contains the function to implement Algorithm 3. This tracks time-varying subspaces, deals with noise, and provides a "federated, over the air implementation".
5. st_miss_fed.m: this is the wrapper to generate data for st-miss problem, and implement FedSTMiss.
6. simple_evd, ccgls, calc_subspace_error, cgls, phifun.m: helper functions used inside NORST_fed.m
7. PROPACK: some linear algebra routines



