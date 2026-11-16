# PU-Sinc-RBF Interpolation

MATLAB implementation of two-dimensional scattered-data interpolation using the **radial basis function partition-of-unity (RBF-PU)** method, with patchwise **leave-one-out cross-validation (LOOCV)** to select the shape parameter. Sinc, multiquadric (MQ), and Gaussian kernels are supported.

## Files

- `PU_Sinc_RBF_Interpolation.m` — main experiment script.
- `PU_SINC_RBF.m` — local RBF-PU interpolation and LOOCV selection.
- `f1.m`, `f2.m`, `f3.m` — test functions (the main script currently uses `f2`).
- `DistanceMatrix.m`, `MakeSDGrid.m` — distance and evaluation-grid utilities.
- `IntegerBasedStructure.m`, `IntegerBasedContainingQuery.m`, `IntegerBasedNeighbourhood.m`, `IntegerBasedRangeSearch.m` — spatial-search utilities.

## Run

Place all **11 `.m` files** in the same MATLAB working directory and execute:

```matlab
PU_Sinc_RBF_Interpolation
```

Change the experiment settings at the top of the main script. Results are saved to `Interpolation_Results.csv` (RMSE, maximum absolute error `LINF`, CPU time, selected shape parameters, and conditioning diagnostics).

**Requirements:** MATLAB, Statistics and Machine Learning Toolbox (`haltonset`), and Signal Processing Toolbox (`sinc`). Global Optimization Toolbox is needed only for the optional `loocv_ga` mode.

**Numerical note:** In the supplied code, the final local interpolation systems are solved without regularization, while the LOOCV objective uses a small diagonal regularization.

## References and third-party code

1. Cavoretto, R., De Rossi, A., Lancellotti, S., and Perracchione, E. (2022). *Software Implementation of the Partition of Unity Method*. **DOI:** [10.14658/PUPJ-DRNA-2022-2-4](https://doi.org/10.14658/PUPJ-DRNA-2022-2-4). Source of the `IntegerBased*.m` routines; [original implementation](https://github.com/sandro-lancellotti/PU).
2. Fasshauer, G. E. (2007). *Meshfree Approximation Methods with MATLAB*. World Scientific. **DOI:** [10.1142/6437](https://doi.org/10.1142/6437). Source credited in `DistanceMatrix.m`.
3. Rippa, S. (1999). *An algorithm for selecting a good value for the parameter c in radial basis function interpolation*. *Advances in Computational Mathematics*, 11, 193–210. **DOI:** [10.1023/A:1018975909870](https://doi.org/10.1023/A:1018975909870). LOOCV reference.

Third-party routines remain subject to their original attribution and redistribution terms.
