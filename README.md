# tesser_successor

Code for analyzing data from the Tesser study using Successor Representation modeling.

## Modules

All code for reading in behavioral data and simulations are under tesser.
* tesser
  * network.py - defines the network used for tesser
  * util.py - functions to read behavioral data (structure learning, induction)
  * sr.py - core model functions for learning/running experiment phases
  * fit.py - model fitting/parameter optimization
  * rsa.py - creating and testing model RDMs
