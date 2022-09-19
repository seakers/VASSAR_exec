package seakers.vassarexecheur.search.intialization.assignment;

import org.moeaframework.core.Initialization;
import org.moeaframework.core.PRNG;
import org.moeaframework.core.Solution;
import org.moeaframework.core.variable.BinaryVariable;
import org.moeaframework.core.variable.EncodingUtils;
import seakers.vassarexecheur.search.problems.assigning.AssigningProblem;

public class LowerInstrumentCountInitialization implements Initialization {

    private final AssigningProblem problem;
    private final double decisionProbability;
    private final int populationSize;

    public LowerInstrumentCountInitialization(AssigningProblem problem, double decisionProbability, int populationSize) {
        this.problem = problem;
        this.decisionProbability = decisionProbability;
        this.populationSize = populationSize;
    }

    @Override
    public Solution[] initialize() {
        Solution[] initialPopulation = new Solution[populationSize];

        for (int i = 0; i < populationSize; i++) {
            Solution newSolution = problem.newSolution();

            for (int j = 0; j < newSolution.getNumberOfVariables(); j++) {
                BinaryVariable var = new BinaryVariable(1);
                boolean value = PRNG.nextDouble() < decisionProbability;
                EncodingUtils.setBoolean(var, value);
                newSolution.setVariable(j, var);
            }

            initialPopulation[i] = newSolution;
        }

        return initialPopulation;
    }
}
