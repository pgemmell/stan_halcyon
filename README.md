# Linking genotype and phenotype evolution using Halcyon

The Halcyon [PhyloG2P](https://doi.org/10.1016/j.tree.2020.01.005) method jointly models the evolution of a continuous trait (e.g. beak length) and a corresponding multiple sequence alignment (e.g. an enhancer) across a phylogeny.

The input for the model is a multiple sequence alignment (parameter **X**) and a background phylogeny and substitution process (parameters **Q**, **pi_distro**, **edge**, **edge_length**, and **path_to_leaf**). Halcyon reports on the association (if any) between the rate of trait evolution and the rate of nucleotide substitution (check the posterior of parameter b -- does it look different from 0?). The model uses a logistic function (parameters h, a, b) to flexibly describe the relationship between the rate of trait evolution and the rate of sequence evolution across the tree. Departures from neutral molecular evolution are modelled flexibly using per-branch substitution rate multipliers (parameter **r**). 

Halcyon is implemented using the probabilistic programming language [Stan](https://mc-stan.org/).

## Installation instructions

[Setup](https://mc-stan.org/users/interfaces/) Stan according to your OS and interface preferences. Use the [stan_halcyon_logistic.stan](https://github.com/pgemmell/stan_halcyon/blob/main/stan_halcyon_logistic.stan) model script and optionally modify the variance function, the priors, or other aspects of the model to your satisfaction. Example simulation scripts are included in the repository, and you can use these as a template for getting started with the model.

## Further information

A document describing the method is coming soon. Watch this space!
