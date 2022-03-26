package seakers.vassarexecheur.search;

import org.moeaframework.core.*;
import org.moeaframework.core.operator.RandomInitialization;
import org.moeaframework.problem.AbstractProblem;
import org.moeaframework.util.TypedProperties;
import seakers.architecture.util.IntegerVariable;
import seakers.vassarexecheur.search.intialization.SynchronizedMersenneTwister;
import seakers.vassarexecheur.search.operators.assigning.*;
import seakers.vassarexecheur.search.operators.partitioning.*;
import seakers.vassarexecheur.search.problems.assigning.AssigningArchitecture;
import seakers.vassarexecheur.search.problems.assigning.AssigningProblem;
import seakers.vassarexecheur.search.problems.partitioning.PartitioningArchitecture;
import seakers.vassarexecheur.search.problems.partitioning.PartitioningProblem;
import seakers.vassarheur.BaseParams;
import seakers.vassarheur.evaluation.AbstractArchitectureEvaluator;
import seakers.vassarheur.evaluation.ArchitectureEvaluationManager;
import seakers.vassarheur.problems.Assigning.ArchitectureEvaluator;
import seakers.vassarheur.problems.Assigning.ClimateCentricAssigningParams;
import seakers.vassarheur.problems.PartitioningAndAssigning.Architecture;
import seakers.vassarheur.problems.PartitioningAndAssigning.ClimateCentricPartitioningParams;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.StringJoiner;

/**
 * Checks to make sure each heuristic operator improves the heuristic satisfaction of a design. This is done by generating random designs for each run and
 * computing heuristic function value before and after manipulation by the heuristic operator.
 */

public class CheckOperatorHeuristicImprovement {

    public static void main(String[] args) throws IOException {
        // Define problem parameters
        boolean assigningProblem = true; // True -> assigning problem, False -> partitioning problem

        double dcThreshold = 0.5;
        double massThreshold = 3000.0; // [kg]
        double packEffThreshold = 0.4;
        boolean considerFeasibility = true;

        int numCPU = 1;

        // Get time
        String timestamp = new SimpleDateFormat("yyyy-MM-dd-HH-mm").format(new Date());

        TypedProperties properties = new TypedProperties();

        String resourcesPath = "C:\\SEAK Lab\\SEAK Lab Github\\VASSAR\\VASSAR_resources-heur";
        String savePath = System.getProperty("user.dir") + File.separator + "results";

        // Heuristic Enforcement Methods
        /**
         * dutyCycleConstrained = [interior_penalty, AOS, biased_init, ACH, objective, constraint]
         * instrumentOrbitRelationsConstrained = [interior_penalty, AOS, biased_init, ACH, objective, constraint]
         * interferenceConstrained = [interior_penalty, AOS, biased_init, ACH, objective, constraint]
         * packingEfficiencyConstrained = [interior_penalty, AOS, biased_init, ACH, objective, constraint]
         * spacecraftMassConstrained = [interior_penalty, AOS, biased_init, ACH, objective, constraint]
         * synergyConstrained = [interior_penalty, AOS, biased_init, ACH, objective, constraint]
         *
         * heuristicsConstrained = [dutyCycleConstrained, instrumentOrbitRelationsConstrained, interferenceConstrained, packingEfficiencyConstrained, spacecraftMassConstrained, synergyConstrained]
         */
        boolean[] dutyCycleConstrained = {false, false, false, false, false, false};
        boolean[] instrumentOrbitRelationsConstrained = {false, false, false, false, false, false};
        boolean[] interferenceConstrained = {false, false, false, false, false, false};
        boolean[] packingEfficiencyConstrained = {false, false, false, false, false, false};
        boolean[] spacecraftMassConstrained = {false, false, false, false, false, false};
        boolean[] synergyConstrained = {false, false, false, false, false, false};

        boolean[][] heuristicsConstrained = new boolean[6][6];
        for (int i = 0; i < 6; i++) {
            heuristicsConstrained[0][i] = dutyCycleConstrained[i];
            heuristicsConstrained[1][i] = instrumentOrbitRelationsConstrained[i];
            heuristicsConstrained[2][i] = interferenceConstrained[i];
            heuristicsConstrained[3][i] = packingEfficiencyConstrained[i];
            heuristicsConstrained[4][i] =  spacecraftMassConstrained[i];
            heuristicsConstrained[5][i] = synergyConstrained[i];
        }

        int numberOfHeuristicConstraints = 0;
        int numberOfHeuristicObjectives = 0;
        for (int i = 0; i < 6; i++) {
            if (heuristicsConstrained[i][5]) {
                numberOfHeuristicConstraints++;
            }
            if (heuristicsConstrained[i][4]) {
                numberOfHeuristicObjectives++;
            }
        }

        BaseParams params;
        AbstractArchitectureEvaluator evaluator;
        if (assigningProblem) {
            params = new ClimateCentricAssigningParams(resourcesPath, "CRISP-ATTRIBUTES","test", "normal");
            evaluator = new ArchitectureEvaluator(considerFeasibility, dcThreshold, massThreshold, packEffThreshold);
        } else {
            params = new ClimateCentricPartitioningParams(resourcesPath, "CRISP-ATTRIBUTES", "test", "normal");
            evaluator = new seakers.vassarheur.problems.PartitioningAndAssigning.ArchitectureEvaluator(considerFeasibility, dcThreshold, massThreshold, packEffThreshold);
        }
        ArchitectureEvaluationManager evaluationManager = new ArchitectureEvaluationManager(params, evaluator);
        evaluationManager.init(numCPU);

        String fileSaveNameProblem  = "";
        if (assigningProblem) {
            fileSaveNameProblem = "_assigning";
        } else {
            fileSaveNameProblem = "_partitioning";
        }

        int numArchs = 100;

        PRNG.setRandom(new SynchronizedMersenneTwister());

        Initialization initialization;

        // Problem class
        AbstractProblem satelliteProblem;
        if (assigningProblem) {
            satelliteProblem = new AssigningProblem(new int[]{1}, params.getProblemName(), evaluationManager, params, dcThreshold, massThreshold, packEffThreshold, numberOfHeuristicObjectives, numberOfHeuristicConstraints, heuristicsConstrained);
        } else {
            satelliteProblem = new PartitioningProblem(params.getProblemName(), evaluationManager, params, dcThreshold, massThreshold, packEffThreshold, numberOfHeuristicObjectives, numberOfHeuristicConstraints, heuristicsConstrained);
        }

        // Initialize heuristic operators
        Variation repairDutyCycle;
        Variation repairInstrumentOrbitRelations;
        Variation repairInterference;
        Variation repairPackingEfficiency;
        Variation repairMass;
        Variation repairSynergy;

        if (assigningProblem) {
            repairDutyCycle = new RepairDutyCycleAssigning(dcThreshold, 1, params);
            repairInstrumentOrbitRelations = new RepairInstrumentOrbitAssigning(1, evaluationManager.getResourcePool(), (ArchitectureEvaluator) evaluator, params);
            repairInterference = new RepairInterferenceAssigning(1, evaluationManager.getResourcePool(), (ArchitectureEvaluator) evaluator, params);
            repairPackingEfficiency = new RepairPackingEfficiencyAssigning(packEffThreshold, 1, params);
            repairMass = new RepairMassAssigning(massThreshold, 1, params);
            repairSynergy = new RepairSynergyAssigning(1, evaluationManager.getResourcePool(), (ArchitectureEvaluator) evaluator, params);
        } else {
            repairDutyCycle = new RepairDutyCyclePartitioning(dcThreshold, 1, params);
            repairInstrumentOrbitRelations = new RepairInstrumentOrbitPartitioning(1, evaluationManager.getResourcePool(), (seakers.vassarheur.problems.PartitioningAndAssigning.ArchitectureEvaluator) evaluator, params);
            repairInterference = new RepairInterferencePartitioning(1, evaluationManager.getResourcePool(), (seakers.vassarheur.problems.PartitioningAndAssigning.ArchitectureEvaluator) evaluator, params);
            repairPackingEfficiency = new RepairPackingEfficiencyAssigning(packEffThreshold, 1, params);
            repairMass = new RepairMassPartitioning(massThreshold, 1, params);
            repairSynergy = new RepairSynergyPartitioning(1, evaluationManager.getResourcePool(), (seakers.vassarheur.problems.PartitioningAndAssigning.ArchitectureEvaluator) evaluator, params);
        }

        // Randomly generate test architectures
        initialization = new RandomInitialization(satelliteProblem, 1);

        int designCounter = 0;

        // Initialize csv file
        String csvFilename = savePath + File.separator + "operator_heuristic_satisfaction" + fileSaveNameProblem + ".csv";
        File csvFile = new File(csvFilename);

        try (FileWriter csvWriter = new FileWriter(csvFile)) {
            // CSV file headers
            StringJoiner headings = new StringJoiner(",");
            headings.add("Full Design - Initial");
            headings.add("Full Design - Duty Cycle");
            headings.add("Duty Cycle - Before");
            headings.add("Duty Cycle - After");
            headings.add("Full Design - Instrument Orbit Relations");
            headings.add("Instrument Orbit Relations - Before");
            headings.add("Instrument Orbit Relations - After");
            headings.add("Full Design - Interference");
            headings.add("Interference - Before");
            headings.add("Interference - After");
            headings.add("Full Design - Packing Efficiency");
            headings.add("Packing Efficiency - Before");
            headings.add("Packing Efficiency - After");
            headings.add("Full Design - Spacecraft Mass");
            headings.add("Spacecraft Mass - Before");
            headings.add("Spacecraft Mass - After");
            headings.add("Full Design - Synergy");
            headings.add("Synergy - Before");
            headings.add("Synergy - After");
            csvWriter.append(headings.toString());
            csvWriter.append("\n");

            while (designCounter < numArchs) {
                Solution currentSolution = initialization.initialize()[0];

                // Evaluate current solution to obtain prior heuristic values
                satelliteProblem.evaluate(currentSolution);
                double currentDutyCycleViolation = (double)currentSolution.getAttribute("DCViolation");
                double currentInstrumentOrbitRelationsViolation = (double)currentSolution.getAttribute("InstrOrbViolation");
                double currentInterferenceViolation = (double)currentSolution.getAttribute("InterInstrViolation");
                double currentPackingEfficiencyViolation = (double)currentSolution.getAttribute("PackEffViolation");
                double currentSpacecraftMassViolation = (double)currentSolution.getAttribute("SpMassViolation");
                double currentSynergyViolation = (double)currentSolution.getAttribute("SynergyViolation");

                // Only accept architectures that can be improved by the heuristic operator
                if ((currentDutyCycleViolation == 0) || (currentInstrumentOrbitRelationsViolation == 0) || (currentInterferenceViolation == 0) || (currentPackingEfficiencyViolation == 0) || (currentSpacecraftMassViolation == 0) || (currentSynergyViolation == 0)) {
                    continue;
                } else {
                    System.out.println("Viable architecture found");

                    // Pass through duty cycle violation operator and evaluate new architecture to obtain post heuristic values
                    Solution dutyCycleSolution = repairDutyCycle.evolve(new Solution[]{currentSolution})[0];
                    satelliteProblem.evaluate(dutyCycleSolution);

                    // Pass through instrument orbit relations violation operator and evaluate new architecture to obtain post heuristic values
                    Solution instrOrbSolution = repairInstrumentOrbitRelations.evolve(new Solution[]{currentSolution})[0];
                    satelliteProblem.evaluate(instrOrbSolution);

                    // Pass through interference violation operator and evaluate new architecture to obtain post heuristic values
                    Solution interferenceSolution = repairInterference.evolve(new Solution[]{currentSolution})[0];
                    satelliteProblem.evaluate(interferenceSolution);

                    // Pass through packing efficiency violation operator and evaluate new architecture to obtain post heuristic values
                    Solution packEffSolution = repairPackingEfficiency.evolve(new Solution[]{currentSolution})[0];
                    satelliteProblem.evaluate(packEffSolution);

                    // Pass through spacecraft mass violation operator and evaluate new architecture to obtain post heuristic values
                    Solution spMassSolution = repairMass.evolve(new Solution[]{currentSolution})[0];
                    satelliteProblem.evaluate(spMassSolution);

                    // Pass through synergy violation operator and evaluate new architecture to obtain post heuristic values
                    Solution synergySolution = repairSynergy.evolve(new Solution[]{currentSolution})[0];
                    satelliteProblem.evaluate(synergySolution);

                    // Populate row of CSV file
                    StringJoiner sj = new StringJoiner(",");

                    String currentSolutionString;
                    String dutyCycleSolutionString;
                    String instrOrbSolutionString;
                    String interferenceSolutionString;
                    String packEffSolutionString;
                    String spMassSolutionString;
                    String synergySolutionString;
                    if (assigningProblem) {
                        currentSolutionString = getBitStringFromAssigningArchitecture(currentSolution);
                        dutyCycleSolutionString = getBitStringFromAssigningArchitecture(dutyCycleSolution);
                        instrOrbSolutionString = getBitStringFromAssigningArchitecture(instrOrbSolution);
                        interferenceSolutionString = getBitStringFromAssigningArchitecture(interferenceSolution);
                        packEffSolutionString = getBitStringFromAssigningArchitecture(packEffSolution);
                        spMassSolutionString = getBitStringFromAssigningArchitecture(spMassSolution);
                        synergySolutionString = getBitStringFromAssigningArchitecture(synergySolution);
                    } else {
                        currentSolutionString = getBitStringFromPartitioningArchitecture(currentSolution, params);
                        dutyCycleSolutionString = getBitStringFromPartitioningArchitecture(dutyCycleSolution, params);
                        instrOrbSolutionString = getBitStringFromPartitioningArchitecture(instrOrbSolution, params);
                        interferenceSolutionString = getBitStringFromPartitioningArchitecture(interferenceSolution, params);
                        packEffSolutionString = getBitStringFromPartitioningArchitecture(packEffSolution, params);
                        spMassSolutionString = getBitStringFromPartitioningArchitecture(spMassSolution, params);
                        synergySolutionString = getBitStringFromPartitioningArchitecture(synergySolution, params);
                    }

                    sj.add(currentSolutionString);

                    sj.add(dutyCycleSolutionString);
                    sj.add(Double.toString(currentDutyCycleViolation));
                    sj.add(Double.toString((double)dutyCycleSolution.getAttribute("DCViolation")));

                    sj.add(instrOrbSolutionString);
                    sj.add(Double.toString(currentInstrumentOrbitRelationsViolation));
                    sj.add(Double.toString((double)instrOrbSolution.getAttribute("InstrOrbViolation")));

                    sj.add(interferenceSolutionString);
                    sj.add(Double.toString(currentInterferenceViolation));
                    sj.add(Double.toString((double)interferenceSolution.getAttribute("InterInstrViolation")));

                    sj.add(packEffSolutionString);
                    sj.add(Double.toString(currentPackingEfficiencyViolation));
                    sj.add(Double.toString((double)packEffSolution.getAttribute("PackEffViolation")));

                    sj.add(spMassSolutionString);
                    sj.add(Double.toString(currentSpacecraftMassViolation));
                    sj.add(Double.toString((double)spMassSolution.getAttribute("SpMassViolation")));

                    sj.add(synergySolutionString);
                    sj.add(Double.toString(currentSynergyViolation));
                    sj.add(Double.toString((double)synergySolution.getAttribute("SynergyViolation")));

                    csvWriter.append(sj.toString());
                    csvWriter.append("\n");

                    designCounter++;
                }

            }
            csvWriter.flush();
        }
    }

    private static String getBitStringFromAssigningArchitecture(Solution sol) {
        AssigningArchitecture arch = (AssigningArchitecture) sol;

        String bitString = "";
        for (int i = 1; i < arch.getNumberOfVariables(); ++i) {
            bitString += arch.getVariable(i).toString();
        }

        return bitString;
    }

    private static String getBitStringFromPartitioningArchitecture(Solution sol, BaseParams params) {
        PartitioningArchitecture arch = (PartitioningArchitecture) sol;

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

        Architecture arch_abs= new Architecture(instrumentPartitioning, orbitAssignment, 1, params);

        return arch_abs.toString("");
    }
}
