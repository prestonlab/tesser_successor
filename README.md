## Code for implementing the SR model on the Tesser fMRI data

This is the repository for implementing successor representation (SR) model on behavioral and neuroimaging data from the Tesser experiment. 
You can find more information about Tesser [here](https://github.com/prestonlab/tesser_successor/wiki).

All code libraries for reading in behavioral data and simulations should be under tesser (e.g., `from tesser import network`). Notebooks can be placed in the main directory. Notebooks should generally not have function definitions; instead, add a function to the relevant module and import that module in your notebook.

Module purposes:
* tesser
  * network.py - defines the network used for tesser
  * util.py - functions to read behavioral data (structure learning, induction)
  * sr.py - core model functions for learning/running experiment phases
  * fit.py - model fitting/parameter optimization
  * rsa.py - creating and testing model RDMs
