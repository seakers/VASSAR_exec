package seakers.vassarexecheur.search.intialization.partitioning;

import org.moeaframework.core.Initialization;
import org.moeaframework.core.Problem;
import org.moeaframework.core.Solution;
import seakers.architecture.util.IntegerVariable;
import seakers.vassarexecheur.search.problems.partitioning.PartitioningArchitecture;
import seakers.vassarheur.architecture.AbstractArchitecture;
import seakers.vassarheur.problems.PartitioningAndAssigning.Architecture;
import seakers.vassarheur.problems.PartitioningAndAssigning.ArchitectureGenerator;
import seakers.vassarheur.problems.PartitioningAndAssigning.PartitioningAndAssigningParams;
import seakers.vassarheur.BaseParams;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class PartitioningInitialization implements Initialization {

    protected final Problem problem;
    protected final int populationSize;
    protected List<Solution> injectedSolutions;
    protected PartitioningAndAssigningParams params;

    private static final boolean debug = false;

    public PartitioningInitialization(Problem problem, int populationSize, BaseParams params){
        this.problem = problem;
        this.populationSize = populationSize;
        this.injectedSolutions = new ArrayList<>();
        this.params = (PartitioningAndAssigningParams) params;
    }

    public PartitioningInitialization(Problem problem, int populationSize, List<Solution> injectedSolution, BaseParams params){
        this.problem = problem;
        this.populationSize = populationSize;
        this.injectedSolutions = injectedSolution;
        this.params = (PartitioningAndAssigningParams) params;
    }

    @Override
    public Solution[] initialize() {

        Solution[] population;

        if (this.populationSize <= this.injectedSolutions.size()) {
            population = this.injectedSolutions.toArray(new Solution[0]);

        }
        else {
            int numRandomlyGeneratedSol = this.populationSize - this.injectedSolutions.size();

            ArchitectureGenerator generator = new ArchitectureGenerator(this.params);
            List<AbstractArchitecture> randomArchs = generator.generateRandomPopulation(numRandomlyGeneratedSol);

            List<Solution> out = new ArrayList<>(injectedSolutions);

            for (AbstractArchitecture arch:randomArchs) {
                seakers.vassarheur.problems.PartitioningAndAssigning.Architecture paArch = (Architecture)arch;

                int[] partitioning = paArch.getInstrumentPartitioning();
                int[] assigning = paArch.getOrbitAssignment();

                PartitioningArchitecture sol = (PartitioningArchitecture) problem.newSolution();
                for (int i = 0; i < partitioning.length; i++) {
                    ((IntegerVariable) sol.getVariable(i)).setValue(partitioning[i]);
                }
                for (int i = 0; i < assigning.length; i++) {
                    ((IntegerVariable) sol.getVariable(i + partitioning.length)).setValue(assigning[i]);
                }
                out.add(sol);
            }
            population = out.toArray(new Solution[0]);
        }

        if (debug) {
            for (Solution sol:population){
                int[] partitioning = new int[params.getNumInstr()];
                int[] assigning = new int[params.getNumInstr()];
                for (int i = 0; i < partitioning.length; ++i) {
                    partitioning[i] = ((IntegerVariable) sol.getVariable(i)).getValue();
                }
                for (int i = 0; i < assigning.length; ++i) {
                    assigning[i] = ((IntegerVariable) sol.getVariable(i + partitioning.length)).getValue();
                }
                System.out.println(Arrays.toString(partitioning) + " | " + Arrays.toString(assigning));
            }
        }

        return population;
    }
}
