# BaM  <a href=""><img src="logo.png" align="right" height="138" /></a>
BaM (Bayesian Modeling) is a framework aimed at estimating a model using a bayesian approach and using it for prediction, with a particular focus on uncertainty quantification.

### Usage
The current code is a computational engine and mostly aims at generating MCMC samples from the posterior distribution associated with the model and the data, and to use these MCMC samples to make predictions. It is implemented as an executable file controlled by a set of configuration text files. BaM is still under active development and the documentation of configuration files is hence incomplete, but useful information may be found [here](https://github.com/BaM-tools/BaMdocs). It is also possible to use BaM with the R package [RBaM](https://github.com/BaM-tools/RBaM). 

### Getting BaM executable
A recent BaM executable for your system (Windows or Linux - MacOS not available yet) can be downloaded [here](https://github.com/BaM-tools/BaM/releases/latest/).

Alternatively, the executable can be recompiled from sources using the provided [makefile](https://github.com/BaM-tools/BaM/blob/main/makefile/makefile) or [Code::Blocks project](https://github.com/BaM-tools/BaM/tree/main/CodeBlocks/BaM). Files from the following projects are needed for this purpose:

1. [BMSL](https://github.com/benRenard/BMSL)
2. [miniDMSL](https://github.com/benRenard/miniDMSL)

Here is the procedure to compile BaM using the makefile: 

1. Download-and-unzip or clone [BaM](https://github.com/BaM-tools/BaM), [BMSL](https://github.com/benRenard/BMSL) and [miniDMSL](https://github.com/benRenard/miniDMSL) into the same folder. You should end up with three folders named BaM, BMSL and miniDMSL sitting alongside (rename them if needed- warning the names are case-sensitive).
2. Open a terminal and move to the makefile directory by typing: `cd /cloning-folder/BaM/makefile`, replacing `cloning-folder` by the folder where you downloaded/cloned the three projects of step 1.
3. type `make` in the terminal.

You should see many messages and after a few minutes at most the executable file `BaM` should be created into the makefile directory.

Developper tools such as the `make` utility and the `gfortran` compiler must be installed on your machine for this procedure to work. This should be natively the case on most computers, but if not you can install the [GNU compiler selection](https://gcc.gnu.org/install/binaries.html).
