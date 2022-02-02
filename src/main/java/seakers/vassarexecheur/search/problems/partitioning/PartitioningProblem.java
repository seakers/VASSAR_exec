package seakers.vassarexecheur.search.problems.partitioning;

import org.moeaframework.core.Solution;
import org.moeaframework.problem.AbstractProblem;
import seakers.architecture.problem.SystemArchitectureProblem;
import seakers.architecture.util.IntegerVariable;
import seakers.vassarheur.Result;
import seakers.vassarheur.architecture.AbstractArchitecture;
import seakers.vassarheur.evaluation.ArchitectureEvaluationManager;
import seakers.vassarheur.BaseParams;
import seakers.vassarheur.problems.PartitioningAndAssigning.Architecture;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;
import java.util.concurrent.ExecutionException;

public class PartitioningProblem extends AbstractProblem implements SystemArchitectureProblem {

    private final String problem;
    private final ArchitectureEvaluationManager evaluationManager;
    private final BaseParams params;
    private final double dcThreshold;
    private final double massThreshold; //[kg]
    private final double packingEffThreshold;

    public PartitioningProblem(String problem, ArchitectureEvaluationManager evalManager, BaseParams params, double dcThreshold, double massThreshold, double packingEfficiencyThreshold) {
        super(2 * params.getNumInstr(),2);
        this.problem = problem;
        this.evaluationManager = evalManager;
        this.params = params;
        this.dcThreshold = dcThreshold;
        this.massThreshold = massThreshold;
        this.packingEffThreshold = packingEfficiencyThreshold;
    }

    @Override
    public void evaluate(Solution solution) {
        PartitioningArchitecture arch = (PartitioningArchitecture) solution;
        evaluateArch(arch);
    }

    public void evaluateArch (PartitioningArchitecture arch) {
        if (!arch.getAlreadyEvaluated()) {
            int numPartitioningVariables = params.getNumInstr();
            int numAssignmentVariables = params.getNumInstr();

            int[] instrumentPartitioning = new int[numPartitioningVariables];
            int[] orbitAssignment = new int[numAssignmentVariables];

            for (int i = 0; i < numPartitioningVariables; i++) {
                instrumentPartitioning[i] = ((IntegerVariable)arch.getVariable(i)).getValue();
            }

            for (int i = 0; i < numAssignmentVariables; i++) {
                orbitAssignment[i] = ((IntegerVariable) arch.getVariable(numPartitioningVariables + i)).getValue();
            }

            // Check constraint
            double constraint = 1.0;
            if (!isFeasible(instrumentPartitioning, orbitAssignment)) {
                constraint = 0.0;
            }
            arch.setConstraint(0, constraint);

            AbstractArchitecture arch_abs;
            if (problem.equalsIgnoreCase("ClimateCentric")) {
                arch_abs= new Architecture(instrumentPartitioning, orbitAssignment, 1, params);
            }
            else {
                throw new IllegalArgumentException("Unrecognizable problem type: " + problem);
            }

            try {
                Result result = evaluationManager.evaluateArchitectureSync(arch_abs, "Slow", dcThreshold, massThreshold, packingEffThreshold);
                arch.setObjective(0, -result.getScience());
                arch.setObjective(1, result.getCost());

                ArrayList<Double> archHeuristics = result.getHeuristics();
                arch.setAttribute("DCViolation",archHeuristics.get(0));
                arch.setAttribute("InstrOrbViolation",archHeuristics.get(1));
                arch.setAttribute("InterInstrViolation",archHeuristics.get(2));
                arch.setAttribute("PackEffViolation",archHeuristics.get(3));
                arch.setAttribute("SpMassViolation",archHeuristics.get(4));
                arch.setAttribute("SynergyViolation",archHeuristics.get(5));

                arch.setOperatorParameters(result.getOperatorParameters());

                arch.setSatellitePayloads(result.getSatellitePayloads());

                arch.setSatelliteOrbits(result.getSatelliteOrbits());
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        }
    }

    @Override
    public Solution newSolution() {
        return new PartitioningArchitecture(params.getNumInstr(), params.getNumOrbits(), 2);
    }

    private boolean isFeasible(int[] instrumentPartitioning, int[] orbitAssignment){

        // Check if the number of satellites matches the number of orbit assignments
        Set<Integer> satIndices = new HashSet<>();
        for (Integer instrumentPartition: instrumentPartitioning) {
            satIndices.add(instrumentPartition);
        }

        for (int i = 0; i < orbitAssignment.length; i++) {
            if (orbitAssignment[i] < 0) {
                if (satIndices.size() != i) {
                    return false;
                }
            }
        }

        // Check if the index of the new satellite is +1 of the largest index
        int max = 0;
        for (Integer instrumentPartition: instrumentPartitioning) {
            if (instrumentPartition > max) {
                if (instrumentPartition == max + 1) {
                    max = instrumentPartition;
                }
                else {
                    return false;
                }
            }
        }

        return true;
    }
}
