# Using Halcyon to analyze the rate of genotypic and phenotypic evolution 

The Halcyon [PhyloG2P](https://doi.org/10.1016/j.tree.2020.01.005) method jointly models the evolution of a continuous trait (e.g. beak length) and a corresponding multiple sequence alignment (e.g. an enhancer) across a phylogeny.

The input for the model is a multiple sequence alignment (parameter **X**) and a background phylogeny and substitution process (parameters **Q**, **pi_distro**, **edge**, **edge_length**, and **path_to_leaf**). Halcyon reports on the association (if any) between the rate of trait evolution and the rate of nucleotide substitution (check the posterior of parameter b -- does it look different from 0?). The model uses a logistic function (parameters h, a, b) to flexibly describe the relationship between the rate of trait evolution and the rate of sequence evolution across the tree. Departures from neutral molecular evolution are modelled flexibly using per-branch substitution rate multipliers (parameter **r**). 

Halcyon is implemented using the probabilistic programming language [Stan](https://mc-stan.org/).

## Example model configurations

The following panels show some possible parameter values under the Halcyon logistic model. The main plot shows how the rate of trait evolution (vertical) is related to the rate of nucleotide evolution (horizontal) via the logistic function. The tree in each panel shows how the rate multipliers scale the nucleotide tree (the branch lengths themselves) and the trait tree (represented by grey triangles).

<img width="760" alt="fig_model_config" src="https://github.com/pgemmell/stan_halcyon/assets/20600944/46510ba1-3975-404a-b2b0-ae11666a7686">

## Model performance

Does the model work? Simulations to recover the direction of association between the rate of change of a trait and the rate of change of a 120bp sequence suggest parameter b can be recovered.

<img width="678" alt="image" src="https://github.com/pgemmell/stan_halcyon/assets/20600944/a5bba52f-340e-42db-926d-d5a57dee2637">

## Installation instructions

[Setup Stan](https://mc-stan.org/users/interfaces/) according to your OS and interface preferences. Use the [stan_halcyon_logistic.stan](https://github.com/pgemmell/stan_halcyon/blob/main/stan_halcyon_logistic.stan) model script and optionally modify the variance function, the priors, or other aspects of the model to your satisfaction. Example simulation scripts are included in the repository, and you can use these as a template for getting started with the model.

## Further information

A [document](https://www.biorxiv.org/content/10.1101/2024.01.23.576950v1) describing the method is available at bioRxiv. 
