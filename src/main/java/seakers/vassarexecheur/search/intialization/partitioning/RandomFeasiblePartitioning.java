package seakers.vassarexecheur.search.intialization.partitioning;

import org.moeaframework.core.Initialization;
import org.moeaframework.core.PRNG;
import org.moeaframework.core.Solution;
import seakers.architecture.util.IntegerVariable;
import seakers.vassarexecheur.search.problems.partitioning.PartitioningProblem;

public class RandomFeasiblePartitioning implements Initialization {
    private final String[] instrumentList;
    private final String[] orbitList;
    private final int populationSize;
    private final PartitioningProblem problem;

    public RandomFeasiblePartitioning (int populationSize, PartitioningProblem problem, String[] instrumentList, String[] orbitList) {
        this.populationSize = populationSize;
        this.problem = problem;
        this.instrumentList = instrumentList;
        this.orbitList = orbitList;
    }

    @Override
    public Solution[] initialize() {
        Solution[] initialPopulation = new Solution[populationSize];

        for (int i = 0; i < populationSize; i++) {
            Solution newSolution = problem.newSolution();

            int maxPartition = 1;

            IntegerVariable newVar = new IntegerVariable(0, 0, instrumentList.length-1);
            newSolution.setVariable(0, newVar);
            for (int j = 1; j < instrumentList.length; j++) {
                int currentInstrumentPartition = PRNG.nextInt(maxPartition + 1);
                newVar = new IntegerVariable(currentInstrumentPartition, 0, instrumentList.length-1);
                newSolution.setVariable(j, newVar);
                if (currentInstrumentPartition == maxPartition) {
                    maxPartition++;
                }
            }

            for (int k = 0; k < instrumentList.length; k++) {
                int currentOrbitAssignment = -1;
                if (k < maxPartition) {
                    currentOrbitAssignment = PRNG.nextInt(orbitList.length);
                }
                newVar = new IntegerVariable(currentOrbitAssignment, -1, orbitList.length-1);
                newSolution.setVariable(instrumentList.length + k, newVar);
            }
            initialPopulation[i] = newSolution;
        }
        return initialPopulation;
    }
}
