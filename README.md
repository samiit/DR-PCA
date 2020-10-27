# DR-PCA
Combining Data Reconciliation with Principal Component Analysis

This repository contains all the code required to reproduce the results of the paper by [Narasimhan and Bhatt (2015)](https://www.sciencedirect.com/science/article/abs/pii/S0098135415000873). The Julia code is best viewed using Pluto.

The article cited above draws a wonderful connection between the well known principal component analysis (PCA) with a less well-known but very useful technique called [Data Reconcilation (DR)](https://en.wikipedia.org/wiki/Data_validation_and_reconciliation). The example code shows how DR can be applied even when the (linear) process model is unknown.

Mostly in PCA, we concentrate on the principal components that contain the maximum variance. While these components do explain the variability of the data, the oft neglected principal components containing minimum variance actually reveal the (linear) model of the system - in short, these define the linear equations from which the data was probably generated. For non-linear cases, the same procedure can be extended using kernel-PCA, with RBF kernels able to capture any level of polynomial non-linearity. The approach works though, if you have large enough data (much greater than the number of possible constraint equations) for a steady state process.

The paper goes in depth for solving a flow process as outlined in the figure below, involving flow balance equations. But the same approach can be extended for any process with all or partial measurements of the variables.
![Schematic of a flow process](https://github.com/samiit/DR-PCA/blob/main/Flow_process.png)

Hopefully, it helps some of you!
