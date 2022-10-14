package seakers.vassarexecheur.search;

import org.moeaframework.core.Initialization;
import org.moeaframework.core.PRNG;
import org.moeaframework.core.Solution;
import org.moeaframework.core.Variation;
import org.moeaframework.core.variable.BinaryVariable;
import org.moeaframework.problem.AbstractProblem;
import org.moeaframework.util.TypedProperties;
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
import seakers.vassarheur.problems.Assigning.Architecture;
import seakers.vassarheur.problems.Assigning.ArchitectureEvaluator;
import seakers.vassarheur.problems.Assigning.AssigningParams;
import seakers.vassarheur.problems.Assigning.ClimateCentricAssigningParams;
import seakers.vassarheur.problems.PartitioningAndAssigning.ClimateCentricPartitioningParams;
import seakers.vassarheur.problems.PartitioningAndAssigning.PartitioningAndAssigningParams;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.*;

/**
 * Checks to make sure each heuristic operator improves the heuristic satisfaction of a design. The csv file with the
 * architectures that don't satisfy any heuristics is read and passed to all the repair operators. Then the heuristic
 * violations before and after the operator manipulation is compared and saved to another csv.
 *
 * @author roshansuresh
 */

public class CheckHeuristicImprovement2 {

    public static void main (String[] args) throws IOException {
        // Define problem parameters
        boolean assigningProblem = true; // True -> assigning problem, False -> partitioning problem

        double dcThreshold = 0.5;
        double massThreshold = 3000.0; // [kg]
        double packEffThreshold = 0.4;
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
         *      heuristicsConstrained = [dutyCycleConstrained, instrumentOrbitRelationsConstrained, interferenceConstrained, packingEfficiencyConstrained, spacecraftMassConstrained, synergyConstrained]
         *  else:
         *      heuristicsConstrained = [dutyCycleConstrained, instrumentOrbitRelationsConstrained, interferenceConstrained, packingEfficiencyConstrained, spacecraftMassConstrained, synergyConstrained, instrumentCountConstrained]
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

        for (int i = 0; i < heuristicsConstrained[0].length; i++) {
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
        for (int i = 0; i < heuristicsConstrained.length; i++) {
            if (heuristicsConstrained[i][5]) {
                numberOfHeuristicConstraints++;
            }
            if (heuristicsConstrained[i][4]) {
                numberOfHeuristicObjectives++;
            }
        }

        AssigningParams assigningParams = null;
        PartitioningAndAssigningParams partitionParams = null;
        AbstractArchitectureEvaluator evaluator;
        int numberOfSatellites = 0;
        ArchitectureEvaluationManager evaluationManager;
        HashMap<String, String[]> instrumentSynergyMap;
        HashMap<String, String[]> interferingInstrumentsMap;

        if (assigningProblem) {
            assigningParams = new ClimateCentricAssigningParams(resourcesPath, "CRISP-ATTRIBUTES","test", "normal");

            instrumentSynergyMap = getInstrumentSynergyNameMap(assigningParams);
            interferingInstrumentsMap = getInstrumentInterferenceNameMap(assigningParams);

            evaluator = new ArchitectureEvaluator(considerFeasibility, interferingInstrumentsMap, instrumentSynergyMap, dcThreshold, massThreshold, packEffThreshold);
            numberOfSatellites = assigningParams.getNumSatellites()[0];

            evaluationManager = new ArchitectureEvaluationManager(assigningParams, evaluator);
            evaluationManager.init(numCPU);


        } else {
            partitionParams = new ClimateCentricPartitioningParams(resourcesPath, "CRISP-ATTRIBUTES", "test", "normal");

            instrumentSynergyMap = getInstrumentSynergyNameMap(partitionParams);
            interferingInstrumentsMap = getInstrumentInterferenceNameMap(partitionParams);

            evaluator = new seakers.vassarheur.problems.PartitioningAndAssigning.ArchitectureEvaluator(considerFeasibility, interferingInstrumentsMap, instrumentSynergyMap, dcThreshold, massThreshold, packEffThreshold);
            numberOfSatellites = partitionParams.getNumSatellites()[0];

            evaluationManager = new ArchitectureEvaluationManager(partitionParams, evaluator);
            evaluationManager.init(numCPU);
        }

        String fileSaveNameProblem  = "";
        if (assigningProblem) {
            fileSaveNameProblem = "_assigning";
        } else {
            fileSaveNameProblem = "_partitioning";
        }

        int numArchs = 100;

        // Initialize csv file
        String csvFilename = savePath + File.separator + "operator_heuristic_satisfaction" + fileSaveNameProblem + ".csv";
        File csvFile = new File(csvFilename);

        PRNG.setRandom(new SynchronizedMersenneTwister());

        Initialization initialization;

        // Problem class
        AbstractProblem satelliteProblem;
        if (assigningProblem) {
            satelliteProblem = new AssigningProblem(new int[]{1}, assigningParams.getProblemName(), evaluationManager, (ArchitectureEvaluator) evaluator, assigningParams, interferingInstrumentsMap, instrumentSynergyMap, dcThreshold, massThreshold, packEffThreshold, instrCountThreshold, numberOfHeuristicObjectives, numberOfHeuristicConstraints, heuristicsConstrained);
        } else {
            satelliteProblem = new PartitioningProblem(partitionParams.getProblemName(), evaluationManager, partitionParams, interferingInstrumentsMap, instrumentSynergyMap, dcThreshold, massThreshold, packEffThreshold, numberOfHeuristicObjectives, numberOfHeuristicConstraints, heuristicsConstrained);
        }

        // Initialize heuristic operators
        Variation repairDutyCycle;
        Variation repairInstrumentOrbitRelations;
        Variation repairInterference;
        Variation repairPackingEfficiency;
        Variation repairMass;
        Variation repairSynergy;
        Variation repairInstrumentCount;

        if (assigningProblem) {
            repairDutyCycle = new RepairDutyCycleAssigning(dcThreshold, 1, assigningParams, false, (AssigningProblem) satelliteProblem, evaluationManager.getResourcePool(), (ArchitectureEvaluator) evaluator);
            repairInstrumentOrbitRelations = new RepairInstrumentOrbitAssigning(1, evaluationManager.getResourcePool(), (ArchitectureEvaluator) evaluator, assigningParams, (AssigningProblem) satelliteProblem, true);
            repairInterference = new RepairInterferenceAssigning(1, evaluationManager.getResourcePool(), (ArchitectureEvaluator) evaluator, assigningParams, (AssigningProblem) satelliteProblem, interferingInstrumentsMap, false);
            //repairPackingEfficiency = new RepairPackingEfficiencyAssigning(packEffThreshold, 1, assigningParams, moveMode, (AssigningProblem) satelliteProblem, evaluationManager.getResourcePool(), (ArchitectureEvaluator) evaluator);
            repairPackingEfficiency = new RepairPackingEfficiencyAdditionAssigning(packEffThreshold, 1, 1, assigningParams, (AssigningProblem) satelliteProblem, evaluationManager.getResourcePool(), (ArchitectureEvaluator) evaluator);
            repairMass = new RepairMassAssigning(massThreshold, 1, assigningParams, false, (AssigningProblem) satelliteProblem, evaluationManager.getResourcePool(), (ArchitectureEvaluator) evaluator);
            //repairSynergy = new RepairSynergyAssigning(1, evaluationManager.getResourcePool(), (ArchitectureEvaluator) evaluator, assigningParams, (AssigningProblem) satelliteProblem, interferingInstrumentsMap, moveMode);
            repairSynergy = new RepairSynergyAdditionAssigning(1, evaluationManager.getResourcePool(), (ArchitectureEvaluator) evaluator, assigningParams, (AssigningProblem) satelliteProblem, instrumentSynergyMap);
            repairInstrumentCount = new RepairInstrumentCountAssigning(1, 1, instrCountThreshold, (AssigningProblem) satelliteProblem, evaluationManager.getResourcePool(), (ArchitectureEvaluator) evaluator, assigningParams);
        } else {
            repairDutyCycle = new RepairDutyCyclePartitioning(dcThreshold, 1, partitionParams, (PartitioningProblem) satelliteProblem, evaluationManager.getResourcePool(), (seakers.vassarheur.problems.PartitioningAndAssigning.ArchitectureEvaluator) evaluator);
            repairInstrumentOrbitRelations = new RepairInstrumentOrbitPartitioning(1, evaluationManager.getResourcePool(), (seakers.vassarheur.problems.PartitioningAndAssigning.ArchitectureEvaluator) evaluator, partitionParams, (PartitioningProblem) satelliteProblem);
            repairInterference = new RepairInterferencePartitioning(1, evaluationManager.getResourcePool(), (seakers.vassarheur.problems.PartitioningAndAssigning.ArchitectureEvaluator) evaluator, partitionParams, (PartitioningProblem) satelliteProblem, interferingInstrumentsMap);
            repairPackingEfficiency = new RepairPackingEfficiencyPartitioning(packEffThreshold, 1, partitionParams, (PartitioningProblem) satelliteProblem, evaluationManager.getResourcePool(), (seakers.vassarheur.problems.PartitioningAndAssigning.ArchitectureEvaluator) evaluator);
            repairMass = new RepairMassPartitioning(massThreshold, 1, partitionParams, (PartitioningProblem) satelliteProblem, evaluationManager.getResourcePool(), (seakers.vassarheur.problems.PartitioningAndAssigning.ArchitectureEvaluator) evaluator);
            repairSynergy = new RepairSynergyPartitioning(1, evaluationManager.getResourcePool(), (seakers.vassarheur.problems.PartitioningAndAssigning.ArchitectureEvaluator) evaluator, partitionParams, (PartitioningProblem) satelliteProblem, instrumentSynergyMap);
        }

        // Extract architectures and heuristic violations as part of the rows of the csv file
        List<List<String>> listOfValues = new ArrayList<>();
        try (Scanner scanner = new Scanner(csvFile);) {
            while (scanner.hasNextLine()) {
                listOfValues.add(getListOfValuesFromLine(scanner.nextLine()));
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        // Pass each architecture through the different operators and store the results in a new csv file
        String newFileName = savePath + File.separator + "operator_heuristic_satisfaction" + fileSaveNameProblem + "_mod" + ".csv";
        File newCsvFile = new File(newFileName);
        int archCounter = 0;

        try(FileWriter csvFileWriter = new FileWriter(newCsvFile)) {
            // File Headers
            StringJoiner headings = new StringJoiner(",");
            for (int i = 0; i < listOfValues.get(0).size(); i++) {
                headings.add(listOfValues.get(0).get(i));
            }
            csvFileWriter.append(headings.toString());
            csvFileWriter.append("\n");

            Solution dutyCycleSolution;
            Solution instrOrbSolution;
            Solution interferenceSolution;
            Solution packEffSolution;
            Solution spMassSolution;
            Solution synergySolution;
            for (int j = 0; j < (listOfValues.size()-1); j++) {
                AssigningArchitecture assign_arch = null;
                PartitioningArchitecture part_arch = null;
                double currentDutyCycleViolation = 0;
                double currentInstrumentOrbitRelationsViolation = 0;
                double currentInterferenceViolation = 0;
                double currentPackingEfficiencyViolation = 0;
                double currentSpacecraftMassViolation = 0;
                double currentSynergyViolation = 0;
                if (assigningProblem) {
                    assign_arch = new AssigningArchitecture(new int[]{1}, assigningParams.getNumInstr(), assigningParams.getNumOrbits(), 2);
                    assign_arch.setVariablesfromString(listOfValues.get(j+1).get(0));
                    satelliteProblem.evaluate(assign_arch);

                    currentDutyCycleViolation = (double)assign_arch.getAttribute("DCViolation");
                    currentInstrumentOrbitRelationsViolation = (double)assign_arch.getAttribute("InstrOrbViolation");
                    currentInterferenceViolation = (double)assign_arch.getAttribute("InterInstrViolation");
                    currentPackingEfficiencyViolation = (double)assign_arch.getAttribute("PackEffViolation");
                    currentSpacecraftMassViolation = (double)assign_arch.getAttribute("SpMassViolation");
                    currentSynergyViolation = (double)assign_arch.getAttribute("SynergyViolation");

                    // Pass through duty cycle violation operator and evaluate new architecture to obtain post heuristic values
                    dutyCycleSolution = repairDutyCycle.evolve(new Solution[]{assign_arch})[0];
                    satelliteProblem.evaluate(dutyCycleSolution);

                    // Pass through instrument orbit relations violation operator and evaluate new architecture to obtain post heuristic values
                    instrOrbSolution = repairInstrumentOrbitRelations.evolve(new Solution[]{assign_arch})[0];
                    satelliteProblem.evaluate(instrOrbSolution);

                    // Pass through interference violation operator and evaluate new architecture to obtain post heuristic values
                    interferenceSolution = repairInterference.evolve(new Solution[]{assign_arch})[0];
                    satelliteProblem.evaluate(interferenceSolution);

                    // Pass through packing efficiency violation operator and evaluate new architecture to obtain post heuristic values
                    packEffSolution = repairPackingEfficiency.evolve(new Solution[]{assign_arch})[0];
                    satelliteProblem.evaluate(packEffSolution);

                    // Pass through spacecraft mass violation operator and evaluate new architecture to obtain post heuristic values
                    spMassSolution = repairMass.evolve(new Solution[]{assign_arch})[0];
                    satelliteProblem.evaluate(spMassSolution);

                    // Pass through synergy violation operator and evaluate new architecture to obtain post heuristic values
                    synergySolution = repairSynergy.evolve(new Solution[]{assign_arch})[0];
                    satelliteProblem.evaluate(synergySolution);

                } else {
                    part_arch = new PartitioningArchitecture(partitionParams.getNumInstr(), partitionParams.getNumOrbits(), 2, partitionParams);
                    int[] instrumentPartitions = part_arch.getInstrumentPartitionsFromString(listOfValues.get(j+1).get(0));
                    int[] orbitAssignments = part_arch.getOrbitAssignmentsFromString(listOfValues.get(j+1).get(0));

                    part_arch.setVariablesFromPartitionArrays(instrumentPartitions, orbitAssignments);
                    satelliteProblem.evaluate(part_arch);

                    currentDutyCycleViolation = (double)part_arch.getAttribute("DCViolation");
                    currentInstrumentOrbitRelationsViolation = (double)part_arch.getAttribute("InstrOrbViolation");
                    currentInterferenceViolation = (double)part_arch.getAttribute("InterInstrViolation");
                    currentPackingEfficiencyViolation = (double)part_arch.getAttribute("PackEffViolation");
                    currentSpacecraftMassViolation = (double)part_arch.getAttribute("SpMassViolation");
                    currentSynergyViolation = (double)part_arch.getAttribute("SynergyViolation");

                    // Pass through duty cycle violation operator and evaluate new architecture to obtain post heuristic values
                    dutyCycleSolution = repairDutyCycle.evolve(new Solution[]{part_arch})[0];
                    satelliteProblem.evaluate(dutyCycleSolution);

                    // Pass through instrument orbit relations violation operator and evaluate new architecture to obtain post heuristic values
                    instrOrbSolution = repairInstrumentOrbitRelations.evolve(new Solution[]{part_arch})[0];
                    satelliteProblem.evaluate(instrOrbSolution);

                    // Pass through interference violation operator and evaluate new architecture to obtain post heuristic values
                    interferenceSolution = repairInterference.evolve(new Solution[]{part_arch})[0];
                    satelliteProblem.evaluate(interferenceSolution);

                    // Pass through packing efficiency violation operator and evaluate new architecture to obtain post heuristic values
                    packEffSolution = repairPackingEfficiency.evolve(new Solution[]{part_arch})[0];
                    satelliteProblem.evaluate(packEffSolution);

                    // Pass through spacecraft mass violation operator and evaluate new architecture to obtain post heuristic values
                    spMassSolution = repairMass.evolve(new Solution[]{part_arch})[0];
                    satelliteProblem.evaluate(spMassSolution);

                    // Pass through synergy violation operator and evaluate new architecture to obtain post heuristic values
                    synergySolution = repairSynergy.evolve(new Solution[]{part_arch})[0];
                    satelliteProblem.evaluate(synergySolution);
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
                if (assigningProblem) {
                    currentSolutionString = assign_arch.getBitString();
                    dutyCycleSolutionString = ((AssigningArchitecture) dutyCycleSolution).getBitString();
                    instrOrbSolutionString = ((AssigningArchitecture) instrOrbSolution).getBitString();
                    interferenceSolutionString = ((AssigningArchitecture) interferenceSolution).getBitString();
                    packEffSolutionString = ((AssigningArchitecture) packEffSolution).getBitString();
                    spMassSolutionString = ((AssigningArchitecture) spMassSolution).getBitString();
                    synergySolutionString = ((AssigningArchitecture) synergySolution).getBitString();
                } else {
                    currentSolutionString = part_arch.getString();
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

                csvFileWriter.append(sj.toString());
                csvFileWriter.append("\n");

                archCounter++;
                System.out.println("Architecture no. " + Integer.toString(archCounter) + " evaluated");
                if (archCounter % 20 == 0) {
                    evaluationManager.getResourcePool().poolClean();
                    System.out.println("Rete clean initiated");
                }
            }
            csvFileWriter.flush();
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

    private static List<String> getListOfValuesFromLine(String line) {
        List<String> values = new ArrayList<String>();
        try (Scanner rowScanner = new Scanner(line)) {
            rowScanner.useDelimiter(",");
            while (rowScanner.hasNext()) {
                values.add(rowScanner.next());
            }
        }
        return values;
    }
}
