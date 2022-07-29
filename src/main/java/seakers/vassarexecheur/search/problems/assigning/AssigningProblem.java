package seakers.vassarexecheur.search.problems.assigning;

import org.moeaframework.core.Solution;
import org.moeaframework.problem.AbstractProblem;
import seakers.architecture.problem.SystemArchitectureProblem;
import seakers.vassarheur.BaseParams;
import seakers.vassarheur.Result;
import seakers.vassarheur.architecture.AbstractArchitecture;
import seakers.vassarheur.evaluation.ArchitectureEvaluationManager;
import seakers.vassarheur.problems.Assigning.Architecture;
import seakers.vassarheur.problems.Assigning.ClimateCentricAssigningParams;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.concurrent.ExecutionException;

public class AssigningProblem  extends AbstractProblem implements SystemArchitectureProblem {

    private final int[] alternativesForNumberOfSatellites;
    private final String problem;
    private final ArchitectureEvaluationManager evalManager;
    private final BaseParams params;
    private final HashMap<String, String[]> interferenceMap;
    private final HashMap<String, String[]> synergyMap;
    private final double dcThreshold;
    private final double massThreshold; //[kg]
    private final double packingEfficiencyThreshold;
    private final int numberOfHeuristicObjectives;
    private final int numberOfHeuristicConstraints;
    private final boolean[][] heuristicsConstrained;

    public AssigningProblem(int[] alternativesForNumberOfSatellites, String problem, ArchitectureEvaluationManager evalManager, BaseParams params, HashMap<String, String[]> interferenceMap, HashMap<String, String[]> synergyMap, double dcThreshold, double massThreshold, double packingEfficiencyThreshold, int numberOfHeuristicObjectives, int numberOfHeuristicConstraints, boolean[][] heuristicsConstrained) {
        super(1 + params.getNumInstr()*params.getNumOrbits(), 2+numberOfHeuristicObjectives, numberOfHeuristicConstraints);
        this.problem = problem;
        this.evalManager = evalManager;
        this.alternativesForNumberOfSatellites = alternativesForNumberOfSatellites;
        this.params = params;
        this.interferenceMap = interferenceMap;
        this.synergyMap = synergyMap;
        this.dcThreshold = dcThreshold;
        this.massThreshold = massThreshold;
        this.packingEfficiencyThreshold = packingEfficiencyThreshold;
        this.numberOfHeuristicObjectives = numberOfHeuristicObjectives;
        this.numberOfHeuristicConstraints = numberOfHeuristicConstraints;
        this.heuristicsConstrained = heuristicsConstrained;
    }

    @Override
    public void evaluate(Solution solution) {
        AssigningArchitecture arch = (AssigningArchitecture) solution;
        evaluateArch(arch);
        //System.out.println(String.format("Arch %s Science = %10f; Cost = %10f", arch.toString(), arch.getObjective(0), arch.getObjective(1)));
    }

    public void evaluateArch(AssigningArchitecture arch) {
        if (!arch.getAlreadyEvaluated()) {
            AbstractArchitecture arch_abs = getAbstractArchitecture(arch);

            try {
                Result result = evalManager.evaluateArchitectureSync(arch_abs, "Slow", interferenceMap, synergyMap, dcThreshold, massThreshold, packingEfficiencyThreshold);
                arch.setObjective(0, -result.getScience());
                arch.setObjective(1, result.getCost());

                ArrayList<Double> archHeuristics = result.getHeuristics();

                for (int i = 0; i < heuristicsConstrained.length; i++) {
                    if (heuristicsConstrained[i][4]) {
                        arch.setObjective(2+i, archHeuristics.get(i));
                    }
                    if (heuristicsConstrained[i][5]) {
                        arch.setConstraint(i, archHeuristics.get(i));
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
        return new AssigningArchitecture(alternativesForNumberOfSatellites, params.getNumInstr(), params.getNumOrbits(), 2+numberOfHeuristicObjectives, numberOfHeuristicConstraints);
    }

    public ArchitectureEvaluationManager getEvaluationManager() {
        return evalManager;
    }

    public AbstractArchitecture getAbstractArchitecture(AssigningArchitecture arch) {
        StringBuilder bitStringBuilder = new StringBuilder(this.getNumberOfVariables());
        for (int i = 1; i < this.getNumberOfVariables(); ++i) {
            bitStringBuilder.append(arch.getVariable(i).toString());
        }

        AbstractArchitecture abs_arch;
        if (problem.equalsIgnoreCase("ClimateCentric")) {
            abs_arch = new Architecture(bitStringBuilder.toString(), ((ClimateCentricAssigningParams) params).getNumSatellites()[0], (ClimateCentricAssigningParams) params);
        }
        else {
            throw new IllegalArgumentException("Unrecognizable problem type: " + problem);
        }
        return abs_arch;
    }
}
