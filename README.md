<div style="display: flex; align-items: center;">
  <div>
    <img src="https://github.com/KyuriP/DepressionDynamics/blob/main/figures/DDlogo.png" alt="Depression Dynamics Logo" width="100">
  </div>
  <div style="margin-left: 15px;">
    <h1 style="margin: 0;">Computational Modeling of Depression Dynamics</h1>
    <p style="margin: 0;">
      This repository contains the analysis scripts and figures used in the paper: 
      <strong>The Individual- and Population-level Mechanistic Implications of Statistical Networks of Symptoms</strong>
    </p>
  </div>
</div>





## Overview

The study explores the mechanisms driving depressive symptom dynamics using a computational model based on network theory. It links mechanistic symptom networks with empirical data and examines how resilience affects symptom interactions. The project provides insights into bistability, hysteresis, and the evolution of statistical networks under different resilience scenarios.



## Repository Structure

### 1. [`code/`](https://github.com/KyuriP/DepressionDynamics/tree/main/code)
- [`libraries.R`](https://github.com/KyuriP/DepressionDynamics/blob/main/code/libraries.R): Script to load required packages.
- [`main_simulation.R`](https://github.com/KyuriP/DepressionDynamics/blob/main/code/main_simulation.R): Runs simulations for various resilience scenarios (baseline, high, low).
- [`robustness_check.R`](https://github.com/KyuriP/DepressionDynamics/blob/main/code/robustness_check.R): Performs robustness checks (Appendix D).
- [`main_figures.R`](https://github.com/KyuriP/DepressionDynamics/blob/main/code/main_figures.R): Creates key figures.
- [`model_specification.R`](https://github.com/KyuriP/DepressionDynamics/blob/main/code/mod_specification.R): Defines the mechanistic model, including parameters and differential equations.


### 2. [`figures/`](https://github.com/KyuriP/DepressionDynamics/tree/main/figures)
- Stores generated plots for publication and analysis.



## Citation

If you use this repository or find the work helpful, please cite:

> Kyuri Park, Lourens Waldorp, and Vítor V. Vasconcelos.  
> "The Individual- and Population-level Mechanistic Implications of Statistical Networks of Symptoms (2024)"  
> [Link to Preprint](#) (update link)



## Contact

For questions, or feedback, please contact [Kyuri Park](mailto:kyurheep@gmail.com).
