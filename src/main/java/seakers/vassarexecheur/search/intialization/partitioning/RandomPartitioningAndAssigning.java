package seakers.vassarexecheur.search.intialization.partitioning;

import org.moeaframework.core.Initialization;
import org.moeaframework.core.PRNG;
import org.moeaframework.core.Solution;
import org.apache.commons.math3.util.CombinatoricsUtils;
import seakers.architecture.util.IntegerVariable;
import seakers.vassarexecheur.search.problems.partitioning.PartitioningProblem;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

public class RandomPartitioningAndAssigning implements Initialization {
    private final String[] instrumentList;
    private final String[] orbitList;
    private final int populationSize;
    private final PartitioningProblem problem;

    public RandomPartitioningAndAssigning(int populationSize, PartitioningProblem problem, String[] instrumentList, String[] orbitList) {
        this.populationSize = populationSize;
        this.problem = problem;
        this.instrumentList = instrumentList;
        this.orbitList = orbitList;
    }

    @Override
    public Solution[] initialize() {
        Solution[] initialPopulation = new Solution[populationSize];

        int numberOfGroups = instrumentList.length;
        long[] sterlingNumbers = new long[numberOfGroups];

        for (int i = 0; i < numberOfGroups; i++) {
            sterlingNumbers[i] = CombinatoricsUtils.stirlingS2(numberOfGroups, i+1);
        }
        long sterlingSum = Arrays.stream(sterlingNumbers).sum();

        double[] cumulativeSubsetProbability = new double[numberOfGroups];
        for (int i = 0; i < numberOfGroups; i++) {
            if (i == 0) {
                cumulativeSubsetProbability[i] = (double) sterlingNumbers[i]/sterlingSum;
            } else {
                cumulativeSubsetProbability[i] = (double) sterlingNumbers[i]/sterlingSum + cumulativeSubsetProbability[i-1];
            }
        }

        int designCounter = 0;
        while (designCounter < 6) {
            Solution oneSubsetArch = problem.newSolution();
            // New 1 subset architecture
            for (int i = 0; i < instrumentList.length; i++) {
                IntegerVariable var = new IntegerVariable(0, 0, instrumentList.length-1);
                oneSubsetArch.setVariable(i, var);
            }
            int orbitIndex = PRNG.nextInt(orbitList.length);
            IntegerVariable orbitVar = new IntegerVariable(orbitIndex, -1, orbitList.length-1);
            oneSubsetArch.setVariable(instrumentList.length, orbitVar);
            initialPopulation[designCounter] = oneSubsetArch;
            designCounter += 1;

            // New n subset architecture
            Solution nSubsetArch = problem.newSolution();
            for (int i = 0; i < instrumentList.length; i++) {
                orbitIndex = PRNG.nextInt(orbitList.length);
                orbitVar = new IntegerVariable(orbitIndex, -1, orbitList.length-1);
                nSubsetArch.setVariable(instrumentList.length+i, orbitVar);
            }
            initialPopulation[designCounter] = oneSubsetArch;
            designCounter += 1;
        }

        while (designCounter < populationSize) {
            double selection = PRNG.nextDouble();
            int subsetValue = getNumberOfSubsets(selection, cumulativeSubsetProbability);
            Solution arch = problem.newSolution();
            randomlyPopulateBySubsets(arch, subsetValue);
            initialPopulation[designCounter] = arch;
            designCounter += 1;
        }

        return initialPopulation;
    }

    private void randomlyPopulateBySubsets (Solution sol, int k) {
        int maxPartition = 1;
        ArrayList<Integer> partitions = new ArrayList<>();
        IntegerVariable var = new IntegerVariable(0, 0, instrumentList.length-1);
        sol.setVariable(0, var);
        partitions.add(0);
        for (int i = 1; i < instrumentList.length; i++) {
            int instrumentPartition = PRNG.nextInt(maxPartition+1);
            partitions.add(instrumentPartition);
            var = new IntegerVariable(instrumentPartition, 0, instrumentList.length-1);
            sol.setVariable(i, var);
            maxPartition = Collections.max(partitions);
            if (maxPartition < k-2) {
                maxPartition += 1;
            }
        }

        for (int j = 0; j <= maxPartition; j++) {
            int partitionOrbit = PRNG.nextInt(orbitList.length);
            var = new IntegerVariable(partitionOrbit, -1, orbitList.length-1);
            sol.setVariable(instrumentList.length+j, var);
        }
    }

    private int getNumberOfSubsets(double randomSelection, double[] cumulativeProbabilities) {
        int subset = 0;
        for (int i = 1; i < cumulativeProbabilities.length; i++) {
            if ((randomSelection > cumulativeProbabilities[i-1]) && (randomSelection < cumulativeProbabilities[i])) {
                subset = i;
                break;
            }
        }
        return subset;
    }

}
