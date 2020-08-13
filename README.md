# Tesser Successor

Methods and code for analyzing data from the Tesser fMRI study using Successor Representation modeling.

## Stimuli

## Experiment

## Data

## Image acquisition

## Analysis software

## Image processing

## Anatomical template and ROIs

## Estimation of item-level activation processing

## Representational similarity analysis

## Successor representation analysis 

All code libraries for reading in behavioral data and simulations should be under tesser (e.g., `from tesser import network`). Notebooks can be placed in the main directory. Notebooks should generally not have function definitions; instead, add a function to the relevant module and import that module in your notebook.

Module purposes:
* tesser
  * network.py - defines the network used for tesser
  * util.py - functions to read behavioral data (structure learning, induction)
  * sr.py - core model functions for learning/running experiment phases
  * fit.py - model fitting/parameter optimization
  * rsa.py - creating and testing model RDMs
