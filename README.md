# BaM  <a href=""><img src="logo.png" align="right" height="138" /></a>
BaM (Bayesian Modeling) is a framework aimed at estimating a model using a bayesian approach and using it for prediction, with a particular focus on uncertainty quantification.

### Usage
The current code is a computational engine and mostly aims at generating MCMC samples from the posterior distribution associated with the model and the data, and to use these MCMC samples to make predictions. It is implemented as an executable file controlled by a set of configuration text files. BaM is still under active development and the documentation of configuration files is hence incomplete, but useful information may be found [here](https://github.com/BaM-tools/BaMdocs). It is also possible to use BaM with the R package [RBaM](https://github.com/BaM-tools/RBaM). 

### Getting BaM executable
A recent BaM executable for your system (Windows or Linux) can be downloaded [here](https://github.com/BaM-tools/RBaM/tree/main/inst/bin).

Alternatively, the executable can be recompiled from sources using the provided makefile or Code::Blocks project. Files from the following projects are needed for this purpose:

1. [BMSL](https://github.com/benRenard/BMSL)
3. [miniDMSL](https://github.com/benRenard/miniDMSL)

