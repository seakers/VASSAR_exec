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
//import seakers.vassarheur.problems.PartitioningAndAssigning.ArchitectureEvaluator;
import seakers.vassarheur.problems.Assigning.ClimateCentricAssigningParams;
import seakers.vassarheur.problems.PartitioningAndAssigning.Architecture;
import seakers.vassarheur.problems.PartitioningAndAssigning.ClimateCentricPartitioningParams;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
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
        double packEffThreshold = 0.7;
        double instrCountThreshold = 15; // only for assigning problem
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
         * instrumentCountConstrained = [interior_penalty, AOS, biased_init, ACH, objective, constraint]
         *
         * if partitioning problem:
         * heuristicsConstrained = [dutyCycleConstrained, instrumentOrbitRelationsConstrained, interferenceConstrained, packingEfficiencyConstrained, spacecraftMassConstrained, synergyConstrained]
         * else:
         * heuristicsConstrained = [dutyCycleConstrained, instrumentOrbitRelationsConstrained, interferenceConstrained, packingEfficiencyConstrained, spacecraftMassConstrained, synergyConstrained, instrumentCountConstrained]
         */
        boolean[] dutyCycleConstrained = {false, false, false, false, false, false};
        boolean[] instrumentOrbitRelationsConstrained = {false, false, false, false, false, false};
        boolean[] interferenceConstrained = {false, false, false, false, false, false};
        boolean[] packingEfficiencyConstrained = {false, false, false, false, false, false};
        boolean[] spacecraftMassConstrained = {false, false, false, false, false, false};
        boolean[] synergyConstrained = {false, false, false, false, false, false};
        boolean[] instrumentCountConstrained = {false, false, false, false, false, false};

        boolean[][] heuristicsConstrained;
        if (assigningProblem) {
            heuristicsConstrained = new boolean[7][6];
        } else {
            heuristicsConstrained = new boolean[6][6];
        }

        for (int i = 0; i < 6; i++) {
            heuristicsConstrained[0][i] = dutyCycleConstrained[i];
            heuristicsConstrained[1][i] = instrumentOrbitRelationsConstrained[i];
            heuristicsConstrained[2][i] = interferenceConstrained[i];
            heuristicsConstrained[3][i] = packingEfficiencyConstrained[i];
            heuristicsConstrained[4][i] =  spacecraftMassConstrained[i];
            heuristicsConstrained[5][i] = synergyConstrained[i];
            if (assigningProblem) {
                heuristicsConstrained[6][i] = instrumentCountConstrained[i];
            }
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

        // Toggle whether assigning operators move on remove instruments
        BaseParams params;
        AbstractArchitectureEvaluator evaluator;
        HashMap<String, String[]> instrumentSynergyMap;
        HashMap<String, String[]> interferingInstrumentsMap;
        if (assigningProblem) {
            params = new ClimateCentricAssigningParams(resourcesPath, "FUZZY-ATTRIBUTES","test", "normal");

            instrumentSynergyMap = getInstrumentSynergyNameMap(params);
            interferingInstrumentsMap = getInstrumentInterferenceNameMap(params);

            evaluator = new ArchitectureEvaluator(considerFeasibility, interferingInstrumentsMap, instrumentSynergyMap, dcThreshold, massThreshold, packEffThreshold);
        } else {
            params = new ClimateCentricPartitioningParams(resourcesPath, "FUZZY-ATTRIBUTES", "test", "normal");

            instrumentSynergyMap = getInstrumentSynergyNameMap(params);
            interferingInstrumentsMap = getInstrumentInterferenceNameMap(params);

            evaluator = new seakers.vassarheur.problems.PartitioningAndAssigning.ArchitectureEvaluator(considerFeasibility, interferingInstrumentsMap, instrumentSynergyMap, dcThreshold, massThreshold, packEffThreshold);
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
            satelliteProblem = new AssigningProblem(new int[]{1}, params.getProblemName(), evaluationManager, (ArchitectureEvaluator) evaluator, params, interferingInstrumentsMap, instrumentSynergyMap, dcThreshold, massThreshold, packEffThreshold, instrCountThreshold, numberOfHeuristicObjectives, numberOfHeuristicConstraints, heuristicsConstrained);
        } else {
            satelliteProblem = new PartitioningProblem(params.getProblemName(), evaluationManager, params, interferingInstrumentsMap, instrumentSynergyMap, dcThreshold, massThreshold, packEffThreshold, numberOfHeuristicObjectives, numberOfHeuristicConstraints, heuristicsConstrained);
        }

        // Initialize heuristic operators
        Variation repairDutyCycle;
        Variation repairInstrumentOrbitRelations;
        Variation repairInterference;
        Variation repairPackingEfficiency;
        Variation repairMass;
        Variation repairSynergy;
        Variation repairInstrumentCount = null;

        if (assigningProblem) {
            repairDutyCycle = new RepairDutyCycleAssigning(dcThreshold, 1, params, false, (AssigningProblem) satelliteProblem, evaluationManager.getResourcePool(), (ArchitectureEvaluator) evaluator);
            repairInstrumentOrbitRelations = new RepairInstrumentOrbitAssigning(1, evaluationManager.getResourcePool(), (ArchitectureEvaluator) evaluator, params, (AssigningProblem) satelliteProblem, true);
            repairInterference = new RepairInterferenceAssigning(1, evaluationManager.getResourcePool(), (ArchitectureEvaluator) evaluator, params, (AssigningProblem) satelliteProblem, interferingInstrumentsMap, false);
            //repairPackingEfficiency = new RepairPackingEfficiencyAssigning(packEffThreshold, 1, params, moveMode, (AssigningProblem) satelliteProblem, evaluationManager.getResourcePool(), (ArchitectureEvaluator) evaluator);
            repairPackingEfficiency = new RepairPackingEfficiencyAdditionAssigning(packEffThreshold, 1, 1, params, (AssigningProblem) satelliteProblem, evaluationManager.getResourcePool(), (ArchitectureEvaluator) evaluator);
            repairMass = new RepairMassAssigning(massThreshold, 1, params, false, (AssigningProblem) satelliteProblem, evaluationManager.getResourcePool(), (ArchitectureEvaluator) evaluator);
            //repairSynergy = new RepairSynergyAssigning(1, evaluationManager.getResourcePool(), (ArchitectureEvaluator) evaluator, params, (AssigningProblem) satelliteProblem, interferingInstrumentsMap, moveMode);
            repairSynergy = new RepairSynergyAdditionAssigning(1, evaluationManager.getResourcePool(), (ArchitectureEvaluator) evaluator, params, (AssigningProblem) satelliteProblem, instrumentSynergyMap);
            repairInstrumentCount = new RepairInstrumentCountAssigning(1, 1, instrCountThreshold, (AssigningProblem) satelliteProblem, evaluationManager.getResourcePool(), (ArchitectureEvaluator) evaluator, params);
        } else {
            repairDutyCycle = new RepairDutyCyclePartitioning(dcThreshold, 1, params, (PartitioningProblem) satelliteProblem, evaluationManager.getResourcePool(), (seakers.vassarheur.problems.PartitioningAndAssigning.ArchitectureEvaluator) evaluator);
            repairInstrumentOrbitRelations = new RepairInstrumentOrbitPartitioning(1, evaluationManager.getResourcePool(), (seakers.vassarheur.problems.PartitioningAndAssigning.ArchitectureEvaluator) evaluator, params, (PartitioningProblem) satelliteProblem);
            repairInterference = new RepairInterferencePartitioning(1, evaluationManager.getResourcePool(), (seakers.vassarheur.problems.PartitioningAndAssigning.ArchitectureEvaluator) evaluator, params, (PartitioningProblem) satelliteProblem, interferingInstrumentsMap);
            repairPackingEfficiency = new RepairPackingEfficiencyPartitioning(packEffThreshold, 1, params, (PartitioningProblem) satelliteProblem, evaluationManager.getResourcePool(), (seakers.vassarheur.problems.PartitioningAndAssigning.ArchitectureEvaluator) evaluator);
            repairMass = new RepairMassPartitioning(massThreshold, 1, params, (PartitioningProblem) satelliteProblem, evaluationManager.getResourcePool(), (seakers.vassarheur.problems.PartitioningAndAssigning.ArchitectureEvaluator) evaluator);
            repairSynergy = new RepairSynergyPartitioning(1, evaluationManager.getResourcePool(), (seakers.vassarheur.problems.PartitioningAndAssigning.ArchitectureEvaluator) evaluator, params, (PartitioningProblem) satelliteProblem, instrumentSynergyMap);
        }

        // Randomly generate test architectures
        initialization = new RandomInitialization(satelliteProblem, 1);

        int archCounter = 0; // viable architectures found
        int archsEvaluated = 0; // total architectures found

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
            if (assigningProblem) {
                headings.add("Full Design - Instrument Count");
                headings.add("Instrument Count - Before");
                headings.add("Instrument Count - After");
            }
            csvWriter.append(headings.toString());
            csvWriter.append("\n");

            while (archCounter < numArchs) {
                Solution currentSolution = initialization.initialize()[0];

                // Evaluate current solution to obtain prior heuristic values
                double currentInstrumentCountViolation = 1.0;
                satelliteProblem.evaluate(currentSolution);
                double currentDutyCycleViolation = (double)currentSolution.getAttribute("DCViolation");
                double currentInstrumentOrbitRelationsViolation = (double)currentSolution.getAttribute("InstrOrbViolation");
                double currentInterferenceViolation = (double)currentSolution.getAttribute("InterInstrViolation");
                double currentPackingEfficiencyViolation = (double)currentSolution.getAttribute("PackEffViolation");
                double currentSpacecraftMassViolation = (double)currentSolution.getAttribute("SpMassViolation");
                double currentSynergyViolation = (double)currentSolution.getAttribute("SynergyViolation");
                if (assigningProblem) {
                    currentInstrumentCountViolation = (double)currentSolution.getAttribute("InstrCountViolation");
                }
                ArrayList<Double> heuristicViolations = new ArrayList<>();
                heuristicViolations.add(currentDutyCycleViolation);
                heuristicViolations.add(currentInstrumentOrbitRelationsViolation);
                heuristicViolations.add(currentInterferenceViolation);
                heuristicViolations.add(currentPackingEfficiencyViolation);
                heuristicViolations.add(currentSpacecraftMassViolation);
                heuristicViolations.add(currentSynergyViolation);
                if (assigningProblem) {
                    heuristicViolations.add(currentInstrumentCountViolation);
                }

                // Only accept architectures that can be improved by the heuristic operator
                if (!checkFullHeuristicSatisfaction(heuristicViolations)) {
                    archsEvaluated++;
                    if (archsEvaluated % 100 == 0) {
                        evaluationManager.getResourcePool().poolClean();
                        System.out.println("Rete clean initiated");
                    }
                    continue;
                } else {
                    System.out.println("Viable architecture found, no. " + (archCounter+1));

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

                    // For assigning problem, pass through instrument count operator and evaluate new architecture to obtain post heuristic values
                    Solution instrCountSolution = null;
                    if (assigningProblem) {
                        instrCountSolution = repairInstrumentCount.evolve(new Solution[]{currentSolution})[0];
                        satelliteProblem.evaluate(instrCountSolution);
                    }

                    // Populate row of CSV file
                    StringJoiner sj = new StringJoiner(",");

                    String currentSolutionString;
                    String dutyCycleSolutionString;
                    String instrOrbSolutionString;
                    String interferenceSolutionString;
                    String packEffSolutionString;
                    String spMassSolutionString;
                    String synergySolutionString;
                    String instrCountSolutionString = null;
                    if (assigningProblem) {
                        currentSolutionString = ((AssigningArchitecture) currentSolution).getBitString();
                        dutyCycleSolutionString = ((AssigningArchitecture) dutyCycleSolution).getBitString();
                        instrOrbSolutionString = ((AssigningArchitecture) instrOrbSolution).getBitString();
                        interferenceSolutionString = ((AssigningArchitecture) interferenceSolution).getBitString();
                        packEffSolutionString = ((AssigningArchitecture) packEffSolution).getBitString();
                        spMassSolutionString = ((AssigningArchitecture) spMassSolution).getBitString();
                        synergySolutionString = ((AssigningArchitecture) synergySolution).getBitString();
                        instrCountSolutionString = ((AssigningArchitecture) instrCountSolution).getBitString();
                    } else {
                        currentSolutionString = ((PartitioningArchitecture) currentSolution).getString();
                        dutyCycleSolutionString = ((PartitioningArchitecture) dutyCycleSolution).getString();
                        instrOrbSolutionString = ((PartitioningArchitecture) instrOrbSolution).getString();
                        interferenceSolutionString = ((PartitioningArchitecture) interferenceSolution).getString();
                        packEffSolutionString = ((PartitioningArchitecture) packEffSolution).getString();
                        spMassSolutionString = ((PartitioningArchitecture) spMassSolution).getString();
                        synergySolutionString = ((PartitioningArchitecture) synergySolution).getString();
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

                    if (assigningProblem) {
                        sj.add(instrCountSolutionString);
                        sj.add(Double.toString(currentInstrumentCountViolation));
                        sj.add(Double.toString((double)instrCountSolution.getAttribute("instrCountViolation")));
                    }

                    csvWriter.append(sj.toString());
                    csvWriter.append("\n");

                    archCounter++;
                    archsEvaluated++;
                    if (archsEvaluated % 100 == 0) {
                        evaluationManager.getResourcePool().poolClean();
                        System.out.println("Rete clean initiated");
                    }
                }
            }
            csvWriter.flush();
            evaluationManager.clear();
            System.out.println("DONE");
        }
    }

    /**
     * Creates instrument synergy map used to compute the instrument synergy violation heuristic (only formulated for the
     * Climate Centric problem for now) (added by roshansuresh)
     * @param params
     * @return Instrument synergy hashmap
     */
    protected static HashMap<String, String[]> getInstrumentSynergyNameMap(BaseParams params) {
        HashMap<String, String[]> synergyNameMap = new HashMap<>();
        if (params.getProblemName().equalsIgnoreCase("ClimateCentric")) {
            synergyNameMap.put("ACE_ORCA", new String[]{"DESD_LID", "GACM_VIS", "ACE_POL", "HYSP_TIR", "ACE_LID"});
            synergyNameMap.put("DESD_LID", new String[]{"ACE_ORCA", "ACE_LID", "ACE_POL"});
            synergyNameMap.put("GACM_VIS", new String[]{"ACE_ORCA", "ACE_LID"});
            synergyNameMap.put("HYSP_TIR", new String[]{"ACE_ORCA", "POSTEPS_IRS"});
            synergyNameMap.put("ACE_POL", new String[]{"ACE_ORCA", "DESD_LID"});
            synergyNameMap.put("ACE_LID", new String[]{"ACE_ORCA", "CNES_KaRIN", "DESD_LID", "GACM_VIS"});
            synergyNameMap.put("POSTEPS_IRS", new String[]{"HYSP_TIR"});
            synergyNameMap.put("CNES_KaRIN", new String[]{"ACE_LID"});
        }
        else {
            System.out.println("Synergy Map for current problem not formulated");
        }
        return synergyNameMap;
    }

    /**
     * Creates instrument interference map used to compute the instrument interference violation heuristic (only formulated for the
     * Climate Centric problem for now)
     * @param params
     * @return Instrument interference hashmap
     */
    protected static HashMap<String, String[]> getInstrumentInterferenceNameMap(BaseParams params) {
        HashMap<String, String[]> interferenceNameMap = new HashMap<>();
        if (params.getProblemName().equalsIgnoreCase("ClimateCentric")) {
            interferenceNameMap.put("ACE_LID", new String[]{"ACE_CPR", "DESD_SAR", "CLAR_ERB", "GACM_SWIR"});
            interferenceNameMap.put("ACE_CPR", new String[]{"ACE_LID", "DESD_SAR", "CNES_KaRIN", "CLAR_ERB", "ACE_POL", "ACE_ORCA", "GACM_SWIR"});
            interferenceNameMap.put("DESD_SAR", new String[]{"ACE_LID", "ACE_CPR"});
            interferenceNameMap.put("CLAR_ERB", new String[]{"ACE_LID", "ACE_CPR"});
            interferenceNameMap.put("CNES_KaRIN", new String[]{"ACE_CPR"});
            interferenceNameMap.put("ACE_POL", new String[]{"ACE_CPR"});
            interferenceNameMap.put("ACE_ORCA", new String[]{"ACE_CPR"});
            interferenceNameMap.put("GACM_SWIR", new String[]{"ACE_LID", "ACE_CPR"});
        }
        else {
            System.out.println("Interference Map fpr current problem not formulated");
        }
        return interferenceNameMap;
    }

    private static boolean checkFullHeuristicSatisfaction(ArrayList<Double> heurViolations) {
        boolean fullySatisfying = true;
        for (int i = 0; i < heurViolations.size(); i++) {
            if (heurViolations.get(i) != 0) {
                fullySatisfying = false;
                break;
            }
        }
        return fullySatisfying;
    }
}
