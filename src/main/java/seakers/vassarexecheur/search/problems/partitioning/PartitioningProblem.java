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
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import java.util.concurrent.ExecutionException;

public class PartitioningProblem extends AbstractProblem implements SystemArchitectureProblem {

    private final String problem;
    private final ArchitectureEvaluationManager evaluationManager;
    private final BaseParams params;
    private final HashMap<String, String[]> interferenceMap;
    private final HashMap<String, String[]> synergyMap;
    private final double dcThreshold;
    private final double massThreshold; //[kg]
    private final double packingEffThreshold;
    private final int numberOfHeuristicObjectives;
    private final int numberOfHeuristicConstraints;
    private final boolean[][] heuristicsConstrained;

    public PartitioningProblem(String problem, ArchitectureEvaluationManager evalManager, BaseParams params, HashMap<String, String[]> interferenceMap, HashMap<String, String[]> synergyMap, double dcThreshold, double massThreshold, double packingEfficiencyThreshold, int numberOfHeuristicObjectives, int numberOfHeuristicConstraints, boolean[][] heuristicsConstrained) {
        super(2 * params.getNumInstr(),2);
        this.problem = problem;
        this.evaluationManager = evalManager;
        this.params = params;
        this.interferenceMap = interferenceMap;
        this.synergyMap = synergyMap;
        this.dcThreshold = dcThreshold;
        this.massThreshold = massThreshold;
        this.packingEffThreshold = packingEfficiencyThreshold;
        this.numberOfHeuristicObjectives = numberOfHeuristicObjectives;
        this.numberOfHeuristicConstraints = numberOfHeuristicConstraints;
        this.heuristicsConstrained = heuristicsConstrained;
    }

    @Override
    public void evaluate(Solution solution) {
        PartitioningArchitecture arch = (PartitioningArchitecture) solution;
        evaluateArch(arch);
        //System.out.println(String.format("Arch %s Science = %10f; Cost = %10f", arch.toString(), arch.getObjective(0), arch.getObjective(1)));
    }

    public void evaluateArch (PartitioningArchitecture arch) {
        if (!arch.getAlreadyEvaluated()) {
            AbstractArchitecture abs_arch = getAbstractArchitecture(arch);

            try {
                Result result = evaluationManager.evaluateArchitectureSync(abs_arch, "Slow", interferenceMap, synergyMap, dcThreshold, massThreshold, packingEffThreshold);

                ArrayList<Double> archHeuristics = result.getHeuristics();

                // Interior Penalization
                double[] objectives = new double[2];
                objectives[0] = -result.getScience()/0.4;
                objectives[1] = result.getCost()/7250;

                // Set objective attributes for Heuristic Coevolutionary Problem
                arch.setAttribute("TrueObjective1", objectives[0]);
                arch.setAttribute("TrueObjective2", objectives[1]);

                double penaltyWeight = 1;
                double heuristicPenalty = 0;
                int numHeuristicInteriorPenalty = 0;
                for (int i = 0; i < heuristicsConstrained.length; i++) {
                    if (heuristicsConstrained[i][0]) {
                         heuristicPenalty += archHeuristics.get(i);
                         numHeuristicInteriorPenalty += 1;
                    }
                }
                if (numHeuristicInteriorPenalty > 0) {
                    heuristicPenalty /= numHeuristicInteriorPenalty;
                }

                objectives[0] += penaltyWeight*heuristicPenalty;
                objectives[1] += penaltyWeight*heuristicPenalty;

                arch.setObjective(0, objectives[0]);
                arch.setObjective(1, objectives[1]);

                for (int i = 0; i < heuristicsConstrained.length; i++) {
                    if (heuristicsConstrained[i][4]) {
                        arch.setObjective(2+i, archHeuristics.get(i));
                    }
                    if (heuristicsConstrained[i][5]) {
                        arch.setConstraint(1+i, archHeuristics.get(i));
                    }
                }

                arch.setAttribute("DCViolation",archHeuristics.get(0));
                arch.setAttribute("InstrOrbViolation",archHeuristics.get(1));
                arch.setAttribute("InterInstrViolation",archHeuristics.get(2));
                arch.setAttribute("PackEffViolation",archHeuristics.get(3));
                arch.setAttribute("SpMassViolation",archHeuristics.get(4));
                arch.setAttribute("SynergyViolation",archHeuristics.get(5));

                arch.setOperatorParameters(result.getOperatorParameters());

                arch.setSatellitePayloads(result.getSatellitePayloads());

                arch.setSatelliteOrbits(result.getSatelliteOrbits());

                arch.setAlreadyEvaluated(true);
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        }
    }

    @Override
    public Solution newSolution() {
        return new PartitioningArchitecture(params.getNumInstr(), params.getNumOrbits(), 2+numberOfHeuristicObjectives, 1+numberOfHeuristicConstraints, params);
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

    public AbstractArchitecture getAbstractArchitecture(PartitioningArchitecture arch) {
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
        return arch_abs;
    }

    public ArchitectureEvaluationManager getEvaluationManager() {
        return evaluationManager;
    }
}
