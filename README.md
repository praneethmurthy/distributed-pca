# Folder contains codes for federated, over the air, PCA

If you use the codes please cite the following paper
[1] "Federated Over-the-Air Subspace Learning from Incomplete Data", Praneeth Narayanamurthy, Namrata Vaswani, Aditya Ramamoorthy, https://arxiv.org/abs/2002.12873 

Here are the list of files:

A: For the fully observed, static subspace setting
1. taubatch_test.m: this file considers the case when we are allowed to vary $\tau$ -- the number of iterations after which we normalize the power method output.
2. sigma_c_tes.m: this file has the codes for varying the channel noise seen at each iteration.
3. ratio_test.m: this file studies the variation in the eigen ratio.


B: For the STmiss problem which may contain missing entries
1. NORST_fed.m: this script contains the function to implement Algorithm 3. This tracks time-varying subspaces, deals with noise, and provides a "federated, over the air implementation".
2. st_miss_fed.m: this is the wrapper to generate data for st-miss problem, and implement FedSTMiss.
3. simple_evd, ccgls, calc_subspace_error, cgls, phifun.m: helper functions used inside NORST_fed.m
4. PROPACK: some linear algebra routines


##need to add real data experiments, should be straightforward, but will have to figure out what to compare with. 
