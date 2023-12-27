# VASSAR_exec

Scripts to generate datasets for the EOSS problems for the paper: 
Suresh Kumar, Roshan, Srikar Srivatsa, Emilie Baker, Meredith Silberstein, and Daniel Selva. "Identifying and Leveraging Promising Design Heuristics for Multi-Objective Combinatorial Design Optimization." Journal of Mechanical Design 145, no. 12 (2023).

## Dependencies:
### JAVA (can be found in the Gradle build.gradle file):
1. MOEAFramework: https://github.com/MOEAFramework/MOEAFramework
2. mopAOS: https://github.com/seakers/mopAOS/tree/heuristics (heuristics branch)
3. System Architecture Problems: https://github.com/seakers/SystemArchitectureProblems
4. VASSAR_lib: https://github.com/seakers/VASSAR_lib/tree/heuristics (heuristics branch)
5. VASSAR_resourecs: https://github.com/seakers/VASSAR_resources/tree/heuristics (heuristics branch)
6. VASSAR_server: https://github.com/seakers/VASSAR_server
7. Jess rules engine: Separate jess.jar file
8. Older version of mopAOS for previous version of the EOSS problems: Separate mopAOS.jar
9. Old EOSS scripts (older version of satellite problems): https://github.com/seakers/EOSS/tree/AIAA (AIAA branch)

### Python:
1. PyGMO (Python Parallel Global Multiobjective Optimizer): https://esa.github.io/pygmo/
2. Scipy
3. Matplotlib

## Important scripts:
### JAVA:
1. MOEARun.java - Start and store results for multiple runs of either problem with different heuristic implementations
2. GenerateForMetricsStudyAssigning.java - Generate datasets for soft constraints screening study for the Assigning problem
3. GenerateForMetricsStudyPartitioning.java - Generate datasets for soft constraints screening study for the Partitioning problem
4. GenerateOperatorIndexDataAssigning.java - Generate datasets for repair operators screening study for the Assigning problem
5. GenerateOperatorIndexDataPartitioning.java - Generate datasets for repair operators screening study for the Partitioning problem

### MATLAB:
1. metrics_study.m - Conduct soft constraints screening study
2. impact_indices_boxplots.m - Plot impact indices boxplots (for both EOSS and Metamaterial Problems)

### PYTHON:
1. hv_satellite_heurcomp.py - Compute hypervolumes and statistics for different cases (efficacy study results)
2. operator_index_computation.py - Compute HDIs for repair operators 
3. biased_sampling_index_computation.py - Compute HDI for Instrument Count Violation biased sampling function
