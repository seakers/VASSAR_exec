Scripts to generate datasets for the EOSS problems for the paper: 
Suresh Kumar, Roshan, Srikar Srivatsa, Emilie Baker, Meredith Silberstein, and Daniel Selva. "Identifying and Leveraging Promising Design Heuristics for Multi-Objective Combinatorial Design Optimization." Journal of Mechanical Design 145, no. 12 (2023).

Important scripts:
JAVA:
MOEARun.java - Start and store results for multiple runs of either problem with different heuristic implementations
GenerateForMetricsStudyAssigning.java - Generate datasets for soft constraints screening study for the Assigning problem
GenerateForMetricsStudyPartitioning.java - Generate datasets for soft constraints screening study for the Partitioning problem
GenerateOperatorIndexDataAssigning.java - Generate datasets for repair operators screening study for the Assigning problem
GenerateOperatorIndexDataPartitioning.java - Generate datasets for repair operators screening study for the Partitioning problem

MATLAB:
metrics_study.m - Conduct soft constraints screening study
impact_indices_boxplots.m - Plot impact indices boxplots

PYTHON:
hv_satellite_heurcomp.py - Compute hypervolumes and statistics for different cases (efficacy study results)
operator_index_computation.py - Compute HDIs for repair operators 
biased_sampling_index_computation.py - Compute HDI for Instrument Count Violation biased sampling function
