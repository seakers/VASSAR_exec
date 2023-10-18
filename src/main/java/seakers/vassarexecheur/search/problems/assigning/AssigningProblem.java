package seakers.vassarexecheur.search.problems.assigning;

import jess.Rete;
import org.moeaframework.core.Solution;
import org.moeaframework.core.variable.EncodingUtils;
import org.moeaframework.problem.AbstractProblem;
import seakers.architecture.problem.SystemArchitectureProblem;
import seakers.vassarheur.BaseParams;
import seakers.vassarheur.QueryBuilder;
import seakers.vassarheur.Resource;
import seakers.vassarheur.Result;
import seakers.vassarheur.architecture.AbstractArchitecture;
import seakers.vassarheur.evaluation.AbstractArchitectureEvaluator;
import seakers.vassarheur.evaluation.ArchitectureEvaluationManager;
import seakers.vassarheur.problems.Assigning.Architecture;
import seakers.vassarheur.problems.Assigning.ArchitectureEvaluator;
import seakers.vassarheur.problems.Assigning.ClimateCentricAssigningParams;
import seakers.vassarheur.utils.MatlabFunctions;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.concurrent.ExecutionException;

public class AssigningProblem  extends AbstractProblem implements SystemArchitectureProblem {

    private final int[] alternativesForNumberOfSatellites;
    private final String problem;
    private final ArchitectureEvaluationManager evalManager;
    private final ArchitectureEvaluator evaluator;
    private final BaseParams params;
    private final HashMap<String, String[]> interferenceMap;
    private final HashMap<String, String[]> synergyMap;
    private final double dcThreshold;
    private final double massThreshold; //[kg]
    private final double packingEfficiencyThreshold;
    private final double instrumentCountThreshold;
    private final int numberOfHeuristicObjectives;
    private final int numberOfHeuristicConstraints;
    private final boolean[][] heuristicsConstrained;

    public AssigningProblem(int[] alternativesForNumberOfSatellites, String problem, ArchitectureEvaluationManager evalManager, ArchitectureEvaluator evaluator, BaseParams params, HashMap<String, String[]> interferenceMap, HashMap<String, String[]> synergyMap, double dcThreshold, double massThreshold, double packingEfficiencyThreshold, double instrumentCountThreshold, int numberOfHeuristicObjectives, int numberOfHeuristicConstraints, boolean[][] heuristicsConstrained) {
        super(1 + params.getNumInstr()*params.getNumOrbits(), 2+numberOfHeuristicObjectives, numberOfHeuristicConstraints);
        this.problem = problem;
        this.evalManager = evalManager;
        this.evaluator = evaluator;
        this.alternativesForNumberOfSatellites = alternativesForNumberOfSatellites;
        this.params = params;
        this.interferenceMap = interferenceMap;
        this.synergyMap = synergyMap;
        this.dcThreshold = dcThreshold;
        this.massThreshold = massThreshold;
        this.packingEfficiencyThreshold = packingEfficiencyThreshold;
        this.instrumentCountThreshold = instrumentCountThreshold;
        this.numberOfHeuristicObjectives = numberOfHeuristicObjectives;
        this.numberOfHeuristicConstraints = numberOfHeuristicConstraints;
        this.heuristicsConstrained = heuristicsConstrained;
    }

    @Override
    public void evaluate(Solution solution) {
        AssigningArchitecture arch = (AssigningArchitecture) solution;
        System.out.println("Architecture: " + arch.getBitString());
        evaluateArch(arch);
        //System.out.println(String.format("Arch %s Science = %10f; Cost = %10f", arch.toString(), arch.getObjective(0), arch.getObjective(1)));
    }

    public void evaluateArch(AssigningArchitecture arch) {
        if (!arch.getAlreadyEvaluated()) {
            AbstractArchitecture arch_abs = getAbstractArchitecture(arch);

            try {
                Result result = new Result(arch_abs, 0.0, 1.0);
                double[] objectives = new double[2];
                int numInstruments = getNumberOfInstruments(arch);
                if (numInstruments > 35) {
                    objectives[0] = 0.0; // Normalized science score (to be maximized)
                    objectives[1] = 1.0; // Normalized cost (to be minimized)

                    Resource res = evalManager.getResourcePool().getResource();
                    Rete r = res.getRete();
                    QueryBuilder qb = res.getQueryBuilder();
                    MatlabFunctions m = res.getM();

                    r.reset();

                    AbstractArchitectureEvaluator eval = evalManager.getEvaluator();
                    evaluator.assertMissions(params, r, arch_abs, m);
                    eval.evaluateHeuristicParameters(r, arch_abs, qb, m);

                    // Compute and add heuristic values to result
                    ArrayList<ArrayList<Double>> archHeuristics = eval.computeHeuristics(r, arch_abs, qb, params);
                    ArrayList<Double> archHeuristicViolations = eval.computeHeuristicsArchitecture(archHeuristics);

                    result.setHeuristics(archHeuristicViolations);

                    // Add Payloads, Orbits and operator parameters to arch
                    ArrayList<ArrayList<String>> payloads = eval.getSatellitePayloads(r, qb);
                    ArrayList<String> orbits = eval.getSatelliteOrbits(r, qb);
                    ArrayList<ArrayList<Double>> operatorParameters = eval.getOperatorParameters(r, qb);

                    arch.setSatelliteOrbits(orbits);
                    arch.setSatellitePayloads(payloads);
                    arch.setOperatorParameters(operatorParameters);

                    evalManager.getResourcePool().freeResource(res);

                } else {
                    result = evalManager.evaluateArchitectureSync(arch_abs, "Slow", interferenceMap, synergyMap, dcThreshold, massThreshold, packingEfficiencyThreshold);

                    objectives[0] = -result.getScience()/0.425;
                    objectives[1] = result.getCost()/2.5e4;

                    arch.setOperatorParameters(result.getOperatorParameters());
                    arch.setSatellitePayloads(result.getSatellitePayloads());
                    arch.setSatelliteOrbits(result.getSatelliteOrbits());
                }

                // Set objective attributes for Heuristic Coevolutionary Problem
                arch.setAttribute("TrueObjective1", objectives[0]);
                arch.setAttribute("TrueObjective2", objectives[1]);

                ArrayList<Double> archHeuristics = result.getHeuristics();

                double numInstrumentsPenalty = getInstrumentCountPenalty2(numInstruments);
                archHeuristics.add(numInstrumentsPenalty);

                // Interior Penalization
                double penaltyWeight = 1;
                double heuristicPenalty = 0;
                int numHeuristicsInteriorPenalty = 0;

                for (int i = 0; i < heuristicsConstrained.length; i++) {
                    if (heuristicsConstrained[i][0]) {
                         heuristicPenalty += archHeuristics.get(i);
                         numHeuristicsInteriorPenalty += 1;
                    }
                }
                if (numHeuristicsInteriorPenalty > 0) {
                    heuristicPenalty /= numHeuristicsInteriorPenalty;
                }

                objectives[0] += penaltyWeight*heuristicPenalty;
                objectives[1] += penaltyWeight*heuristicPenalty;

                arch.setObjective(0, objectives[0]);
                arch.setObjective(1, objectives[1]);

                // Additional Heuristic Objectives and/or Constraints
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
                arch.setAttribute("InstrCountViolation",archHeuristics.get(6));

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

    public int getNumberOfInstruments(AssigningArchitecture arch) {
        int numberOfInstruments = 0;
        for (int i = 0; i < arch.getNumberOfVariables()-1; i++) {
            if (EncodingUtils.getBoolean(arch.getVariable(i+1))) {
                numberOfInstruments++;
            }
        }
        return numberOfInstruments;
    }

    private double getInstrumentCountPenalty(int numInstruments) {
        return 1/(1 + Math.exp(-0.4*(numInstruments - 30)));
    }

    private double getInstrumentCountPenalty2(int numInstruments) {
        if (numInstruments < instrumentCountThreshold) {
            return 0.0;
        } else if (numInstruments < 40) {
            return (numInstruments - instrumentCountThreshold)/(40.0 - instrumentCountThreshold);
        } else {
            return 1.0;
        }
    }
}
