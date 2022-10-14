package seakers.vassarexecheur.search;

import org.moeaframework.core.Initialization;
import org.moeaframework.core.PRNG;
import org.moeaframework.core.Solution;
import org.moeaframework.core.Variation;
import org.moeaframework.core.operator.CompoundVariation;
import org.moeaframework.core.operator.RandomInitialization;
import org.moeaframework.core.operator.binary.BitFlip;
import org.moeaframework.util.TypedProperties;
import seakers.vassarexecheur.search.intialization.SynchronizedMersenneTwister;
import seakers.vassarexecheur.search.operators.assigning.*;
import seakers.vassarexecheur.search.problems.assigning.AssigningProblem;
import seakers.vassarexecheur.search.problems.assigning.AssigningArchitecture;
import seakers.vassarheur.BaseParams;
import seakers.vassarheur.evaluation.AbstractArchitectureEvaluator;
import seakers.vassarheur.evaluation.ArchitectureEvaluationManager;
import seakers.vassarheur.problems.Assigning.ArchitectureEvaluator;
import seakers.vassarheur.problems.Assigning.ClimateCentricAssigningParams;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.*;

public class GenerateOperatorIndexDataAssigning {

    public static void main(String[] args) {
        int numRuns = 10;
        RunMode runMode  = RunMode.RandomPopulation;
        int numCpus = 1;

        //boolean moveInstrument = false; // No longer used, set accordingly in assigning operator instances

        boolean finalPopulationOnly = false;

        // Get time
        String timestamp = new SimpleDateFormat("yyyy-MM-dd-HH-mm").format(new Date());

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
         * heuristicsConstrained = [dutyCycleConstrained, instrumentOrbitRelationsConstrained, interferenceConstrained, packingEfficiencyConstrained, spacecraftMassConstrained, synergyConstrained, instrumentCountConstrained]
         */
        boolean[] dutyCycleConstrained = {false, false, false, false, false, false};
        boolean[] instrumentOrbitRelationsConstrained = {false, false, false, false, false, false};
        boolean[] interferenceConstrained = {false, false, false, false, false, false};
        boolean[] packingEfficiencyConstrained = {false, false, false, false, false, false};
        boolean[] spacecraftMassConstrained = {false, false, false, false, false, false};
        boolean[] synergyConstrained = {false, false, false, false, false, false};
        boolean[] instrumentCountConstrained = {false, false, false, false, false, false};

        boolean[][] heuristicsConstrained = new boolean[7][6];
        for (int i = 0; i < 6; i++) {
            heuristicsConstrained[0][i] = dutyCycleConstrained[i];
            heuristicsConstrained[1][i] = instrumentOrbitRelationsConstrained[i];
            heuristicsConstrained[2][i] = interferenceConstrained[i];
            heuristicsConstrained[3][i] = packingEfficiencyConstrained[i];
            heuristicsConstrained[4][i] =  spacecraftMassConstrained[i];
            heuristicsConstrained[5][i] = synergyConstrained[i];
            heuristicsConstrained[6][i] = instrumentCountConstrained[i];
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

        TypedProperties properties = new TypedProperties();

        int popSize = 300;
        //int maxEvals = 5000;
        //properties.setInt("maxEvaluations", maxEvals);
        properties.setInt("populationSize", popSize);
        double crossoverProbability = 1.0;
        properties.setDouble("crossoverProbability", crossoverProbability);
        double mutationProbability = 1. / 60.;
        properties.setDouble("mutationProbability", mutationProbability);

        //setup for epsilon MOEA
        double[] epsilonDouble = new double[]{0.001, 10};

        double dcThreshold = 0.5;
        double massThreshold = 3000.0; // [kg]
        double packEffThreshold = 0.7;
        double instrCountThreshold = 15; // only for assigning problem
        boolean considerFeasibility = false; // use false only for biased random generation for random population runs

        String savePath = System.getProperty("user.dir") + File.separator + "results";

        String resourcesPath = "C:\\SEAK Lab\\SEAK Lab Github\\VASSAR\\VASSAR_resources-heur"; // for lab system
        //String resourcesPath = "C:\\Users\\rosha\\Documents\\SEAK Lab Github\\VASSAR\\VASSAR_resources-heur"; // for laptop

        ClimateCentricAssigningParams params = new ClimateCentricAssigningParams(resourcesPath, "FUZZY-ATTRIBUTES","test", "normal");

        PRNG.setRandom(new SynchronizedMersenneTwister());

        HashMap<String, String[]> instrumentSynergyMap = getInstrumentSynergyNameMap(params);
        HashMap<String, String[]> interferingInstrumentsMap = getInstrumentInterferenceNameMap(params);

        ArchitectureEvaluator evaluator = new ArchitectureEvaluator(considerFeasibility, interferingInstrumentsMap, instrumentSynergyMap, dcThreshold, massThreshold, packEffThreshold);
        ArchitectureEvaluationManager evaluationManager = new ArchitectureEvaluationManager(params, evaluator);
        evaluationManager.init(numCpus);

        AssigningProblem problem = new AssigningProblem(new int[]{1}, params.getProblemName(), evaluationManager, evaluator, params, interferingInstrumentsMap, instrumentSynergyMap, dcThreshold, massThreshold, packEffThreshold, instrCountThreshold, numberOfHeuristicObjectives, numberOfHeuristicConstraints, heuristicsConstrained);

        Initialization initialization = new RandomInitialization(problem, popSize);

        // duty Cycle, interference, mass -> remove, instrument orbit -> move, pack Eff, synergy -> Add
        Variation repairDutyCycle = new RepairDutyCycleAssigning(dcThreshold, 1, params, false, problem, evaluationManager.getResourcePool(), evaluator);
        Variation repairInstrumentOrbitRelations = new RepairInstrumentOrbitAssigning(1, evaluationManager.getResourcePool(), evaluator, params, problem, true);
        Variation repairInterference = new RepairInterferenceAssigning(1, evaluationManager.getResourcePool(), evaluator, params, problem, interferingInstrumentsMap, false);
        //Variation repairPackingEfficiency = new RepairPackingEfficiencyAssigning(packEffThreshold, 1, params, moveInstrument, problem, evaluationManager.getResourcePool(), evaluator);
        Variation repairPackingEfficiency = new RepairPackingEfficiencyAdditionAssigning(packEffThreshold, 1, 1, params, problem, evaluationManager.getResourcePool(), evaluator);
        Variation repairMass = new RepairMassAssigning(massThreshold, 1, params, false, problem, evaluationManager.getResourcePool(), evaluator);
        //Variation repairSynergy = new RepairSynergyAssigning(1, evaluationManager.getResourcePool(), evaluator, params, problem, instrumentSynergyMap, moveInstrument);
        Variation repairSynergy = new RepairSynergyAdditionAssigning(1, evaluationManager.getResourcePool(), evaluator, params, problem, instrumentSynergyMap);
        Variation repairInstrumentCount = new RepairInstrumentCountAssigning(1, 1, instrCountThreshold, problem, evaluationManager.getResourcePool(), evaluator, params);

        Variation[] heuristicOperatorsArray = {repairDutyCycle, repairInstrumentOrbitRelations, repairInterference, repairPackingEfficiency, repairMass, repairSynergy, repairInstrumentCount};

        switch (runMode) {
            case RandomPopulation:
                System.out.println("Starting random population evaluation for Assigning Problem");
                for (int i = 0; i < numRuns; i++) {
                    String runName = savePath + File.separator + "random_assign_operator_index_" + i + ".csv";

                    Solution[] population = initialization.initialize();
                    int[] nfes = new int[population.length];

                    computeAndStoreIntoCsv(population, runName, heuristicOperatorsArray, problem, true, nfes);
                    System.out.println("Data generated for run " + i);
                }
                evaluationManager.clear();
                System.out.println("DONE");
                break;

            case EpsilonMOEA: // Read completed epsilon MOEA run results, pass archs through operators and store results in new csv file
                System.out.println("Starting Epsilon MOEA data evaluation for Assigning Problem");

                String filepathData = savePath + File.separator + "Epsilon MOEA - Metrics" + File.separator + "Assigning\\";
                String filename = "EpsilonMOEA_emoea_ClimateCentric_assign_";

                for (int i = 0; i < numRuns; i++) {
                    String csvFilename = "";
                    String saveFilename = "";
                    if (finalPopulationOnly) {
                        csvFilename = filepathData + filename + "randinit_" + i + ".csv";
                        saveFilename = savePath + File.separator + filename + "operator_data_" + i + ".csv";
                    } else {
                        csvFilename = filepathData + filename + "randinit_" + i + "_fullpop.csv";
                        saveFilename = savePath + File.separator + filename + "operator_data_" + i + "_fullpop.csv";
                    }

                    ArrayList<Solution> populationList = new ArrayList<>();
                    ArrayList<Integer> nfesList = new ArrayList<>();

                    // Read appropriate csv file
                    List<List<String>> rows = new ArrayList<>();
                    try (Scanner scanner = new Scanner(new File(csvFilename))) {
                        while (scanner.hasNextLine()) {
                            rows.add(getRecordFromLine(scanner.nextLine()));
                        }
                    } catch (FileNotFoundException e) {
                        e.printStackTrace();
                    }

                    boolean header = true;
                    for (List<String> row : rows) {
                        if (header) {
                            header = false;
                            continue;
                        }
                        AssigningArchitecture newArch = (AssigningArchitecture) problem.newSolution();
                        newArch.setVariablesfromString(row.get(0));

                        populationList.add(newArch);
                        if (!finalPopulationOnly) {
                            nfesList.add(Integer.parseInt(row.get(1)));
                        }
                    }
                    Solution[] population = new Solution[populationList.size()];
                    for (int k = 0; k < populationList.size(); k++) {
                        population[k] = populationList.get(k);
                    }
                    boolean randomMode = false;
                    if (finalPopulationOnly) {
                        randomMode = true;
                    }
                    computeAndStoreIntoCsv(population, saveFilename, heuristicOperatorsArray, problem, randomMode, nfesList.stream().mapToInt(j->j).toArray());
                    System.out.println("Data generated for run " + i);
                }
                evaluationManager.clear();
                System.out.println("DONE");
                break;

            default :
                throw new IllegalStateException("Unrecognized run mode");
        }
    }

    private static void computeAndStoreIntoCsv(Solution[] archs, String fileSaveName, Variation[] heuristicOperators, AssigningProblem prob, boolean randomMode, int[] nfes) {

        Variation repairDutyCycle = heuristicOperators[0];
        Variation repairInstrumentOrbitRelations = heuristicOperators[1];
        Variation repairInterference = heuristicOperators[2];
        Variation repairPackingEfficiency = heuristicOperators[3];
        Variation repairMass = heuristicOperators[4];
        Variation repairSynergy = heuristicOperators[5];
        Variation repairInstrumentCount = heuristicOperators[6];

        File csvFile = new File(fileSaveName);
        try (FileWriter writer = new FileWriter(csvFile)) {
            StringJoiner headings = new StringJoiner(",");
            if (!randomMode) {
                headings.add("NFE");
            }
            headings.add("Architecture");
            headings.add("Science Score");
            headings.add("Cost");
            headings.add("Architecture - Instrdc");
            headings.add("Science Score - Instrdc");
            headings.add("Cost - Instrdc");
            headings.add("Architecture - Instrorb");
            headings.add("Science Score - Instrorb");
            headings.add("Cost - Instrorb");
            headings.add("Architecture - Interinstr");
            headings.add("Science Score - Interinstr");
            headings.add("Cost - Interinstr");
            headings.add("Architecture - Packeff");
            headings.add("Science Score - Packeff");
            headings.add("Cost - Packeff");
            headings.add("Architecture - Spmass");
            headings.add("Science Score - Spmass");
            headings.add("Cost - Spmass");
            headings.add("Architecture - Instrsyn");
            headings.add("Science Score - Instrsyn");
            headings.add("Cost - Instrsyn");
            headings.add("Architecture - Instrcount");
            headings.add("Science Score - Instrcount");
            headings.add("Cost - Instrcount");
            writer.append(headings.toString());
            writer.append("\n");

            for (int j = 0; j < archs.length; j++) {

                System.out.println("Architecture " + j + " is being evaluated");

                Solution currentSol = archs[j];
                AssigningArchitecture currentArch = (AssigningArchitecture) currentSol;
                prob.evaluateArch(currentArch);

                Solution instrdcSol = repairDutyCycle.evolve(new Solution[]{currentSol})[0];
                AssigningArchitecture instrdcArch = (AssigningArchitecture) instrdcSol;
                prob.evaluateArch(instrdcArch);

                Solution instrorbSol = repairInstrumentOrbitRelations.evolve(new Solution[]{currentSol})[0];
                AssigningArchitecture instrorbArch = (AssigningArchitecture) instrorbSol;
                prob.evaluateArch(instrorbArch);

                Solution interinstrSol = repairInterference.evolve(new Solution[]{currentSol})[0];
                AssigningArchitecture interinstrArch = (AssigningArchitecture) interinstrSol;
                prob.evaluateArch(interinstrArch);

                Solution packeffSol = repairPackingEfficiency.evolve(new Solution[]{currentSol})[0];
                AssigningArchitecture packeffArch = (AssigningArchitecture) packeffSol;
                prob.evaluateArch(packeffArch);

                Solution spmassSol = repairMass.evolve(new Solution[]{currentSol})[0];
                AssigningArchitecture spmassArch = (AssigningArchitecture) spmassSol;
                prob.evaluateArch(spmassArch);

                Solution instrsynSol = repairSynergy.evolve(new Solution[]{currentSol})[0];
                AssigningArchitecture instrsynArch = (AssigningArchitecture) instrsynSol;
                prob.evaluateArch(instrsynArch);

                Solution instrcountSol = repairInstrumentCount.evolve(new Solution[]{currentSol})[0];
                AssigningArchitecture instrcountArch = (AssigningArchitecture) instrcountSol;
                prob.evaluateArch(instrcountArch);

                prob.getEvaluationManager().getResourcePool().poolClean();

                StringJoiner sj = new StringJoiner(",");
                if (!randomMode) {
                    sj.add(Integer.toString(nfes[j]));
                }
                sj.add(currentArch.getBitString());
                sj.add(Double.toString(currentArch.getObjective(0)));
                sj.add(Double.toString(currentArch.getObjective(1)));
                sj.add(instrdcArch.getBitString());
                sj.add(Double.toString(instrdcArch.getObjective(0)));
                sj.add(Double.toString(instrdcArch.getObjective(1)));
                sj.add(instrorbArch.getBitString());
                sj.add(Double.toString(instrorbArch.getObjective(0)));
                sj.add(Double.toString(instrorbArch.getObjective(1)));
                sj.add(interinstrArch.getBitString());
                sj.add(Double.toString(interinstrArch.getObjective(0)));
                sj.add(Double.toString(interinstrArch.getObjective(1)));
                sj.add(packeffArch.getBitString());
                sj.add(Double.toString(packeffArch.getObjective(0)));
                sj.add(Double.toString(packeffArch.getObjective(1)));
                sj.add(spmassArch.getBitString());
                sj.add(Double.toString(spmassArch.getObjective(0)));
                sj.add(Double.toString(spmassArch.getObjective(1)));
                sj.add(instrsynArch.getBitString());
                sj.add(Double.toString(instrsynArch.getObjective(0)));
                sj.add(Double.toString(instrsynArch.getObjective(1)));
                sj.add(instrcountArch.getBitString());
                sj.add(Double.toString(instrcountArch.getObjective(0)));
                sj.add(Double.toString(instrcountArch.getObjective(1)));

                writer.append(sj.toString());
                writer.append("\n");
            }
            writer.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static List<String> getRecordFromLine(String line) {
        List<String> architectures = new ArrayList<>();
        try (Scanner rowScanner = new Scanner(line)) {
            rowScanner.useDelimiter(",");
            while (rowScanner.hasNext()) {
                architectures.add(rowScanner.next());
            }
        }
        return architectures;
    }

    public enum RunMode{
        RandomPopulation,
        EpsilonMOEA,
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

}
