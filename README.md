# DR-PCA
Combining Data Reconciliation with Principal Component Analysis

This repository contains all the code required to reproduce the results of the paper by [Narasimhan and Bhatt (2015)](https://www.sciencedirect.com/science/article/abs/pii/S0098135415000873). The Julia code is best viewed using Pluto.

The article cited above draws a wonderful connection between the well known principal component analysis (PCA) with a less well-known but very useful technique called [Data Reconcilation (DR)](https://en.wikipedia.org/wiki/Data_validation_and_reconciliation). The example code shows how DR can be applied even when the (linear) process model is unknown.

Mostly in PCA, we concentrate on the principal components that contain the maximum variance. While this explains the variability of the data, the oft neglected principal components containing minimum variance actually reveal the model of the system. For non-linear cases, the same procedure can be extended using kernel-PCA. Hopefully, it helps some of you!
