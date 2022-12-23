package seakers.vassarexecheur.search;

import org.apache.commons.math3.util.CombinatoricsUtils;
import org.moeaframework.algorithm.EpsilonMOEA;
import org.moeaframework.core.*;
import org.moeaframework.core.comparator.ChainedComparator;
import org.moeaframework.core.comparator.ParetoObjectiveComparator;
import org.moeaframework.core.operator.*;
import seakers.architecture.util.IntegerVariable;
import org.moeaframework.util.TypedProperties;
import seakers.vassarexecheur.search.intialization.SynchronizedMersenneTwister;

import seakers.vassarexecheur.search.intialization.partitioning.PartitioningInitialization;
import seakers.vassarexecheur.search.operators.partitioning.PartitioningCrossover;
import seakers.vassarexecheur.search.operators.partitioning.PartitioningMutation;
import seakers.vassarexecheur.search.problems.partitioning.PartitioningArchitecture;
import seakers.vassarexecheur.search.problems.partitioning.PartitioningProblem;
import seakers.vassarheur.BaseParams;
import seakers.vassarheur.evaluation.AbstractArchitectureEvaluator;
import seakers.vassarheur.evaluation.ArchitectureEvaluationManager;
import seakers.vassarheur.problems.PartitioningAndAssigning.Architecture;
import seakers.vassarheur.problems.PartitioningAndAssigning.ArchitectureEvaluator;
import seakers.vassarheur.problems.PartitioningAndAssigning.ClimateCentricPartitioningParams;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.concurrent.*;
import java.util.stream.IntStream;

public class GenerateForMetricsStudyPartitioning {

    public static void main(String[] args) {
        int numRuns = 10;
        int numCpus = 1;

        RunMode runMode  = RunMode.RandomPopulation;
        RandomMode randomMode = RandomMode.FullyRandom;
        InitializationMode initializationMode = InitializationMode.InitializeRandom;

        ExecutorService pool = Executors.newFixedThreadPool(numCpus);
        CompletionService<Algorithm> ecs = new ExecutorCompletionService<>(pool);

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

        TypedProperties properties = new TypedProperties();

        int popSize = 300;
        int maxEvals = 5000;
        properties.setInt("maxEvaluations", maxEvals);
        properties.setInt("populationSize", popSize);
        double crossoverProbability = 1.0;
        properties.setDouble("crossoverProbability", crossoverProbability);
        double mutationProbability = 1. / 24.;
        properties.setDouble("mutationProbability", mutationProbability);
        Variation partCross;
        Variation partFlip;
        Initialization initialization;

        // Set seed for random number generator
        //PRNG.setSeed(4321);

        //setup for epsilon MOEA
        double[] epsilonDouble = new double[]{0.01, 0.01};

        double dcThreshold = 0.5;
        double massThreshold = 3000.0; // [kg]
        double packEffThreshold = 0.7;
        boolean considerFeasibility = false; // use false only for biased random generation for random population runs

        String savePath = System.getProperty("user.dir") + File.separator + "results";

        String resourcesPath = "C:\\SEAK Lab\\SEAK Lab Github\\VASSAR\\VASSAR_resources-heur"; // for lab system
        //String resourcesPath = "C:\\Users\\rosha\\Documents\\SEAK Lab Github\\VASSAR\\VASSAR_resources-heur"; // for laptop

        ClimateCentricPartitioningParams params = new ClimateCentricPartitioningParams(resourcesPath, "FUZZY-ATTRIBUTES", "test", "normal");

        HashMap<String, String[]> instrumentSynergyMap = getInstrumentSynergyNameMap(params);
        HashMap<String, String[]> interferingInstrumentsMap = getInstrumentInterferenceNameMap(params);

        AbstractArchitectureEvaluator evaluator = new ArchitectureEvaluator(considerFeasibility, interferingInstrumentsMap, instrumentSynergyMap, dcThreshold, massThreshold, packEffThreshold);
        ArchitectureEvaluationManager evaluationManager = new ArchitectureEvaluationManager(params, evaluator);
        evaluationManager.init(numCpus);

        PRNG.setRandom(new SynchronizedMersenneTwister());

        switch (runMode) {

            case EpsilonMOEA:
                System.out.println("Starting EpsilonMOEA for Partitioning Problem");
                String initializationString;

                for (int i = 0; i < numRuns; i++) {

                    PartitioningProblem problem = new PartitioningProblem(params.getProblemName(), evaluationManager, params, interferingInstrumentsMap, instrumentSynergyMap, dcThreshold, massThreshold, packEffThreshold, numberOfHeuristicObjectives, numberOfHeuristicConstraints, heuristicsConstrained);

                    switch (initializationMode) {
                        case InitializeRandom:
                            initialization = new PartitioningInitialization(problem, popSize, params);
                            initializationString = "randinit";
                            break;

                        case InitializationRandomAndInjected:
                            // Injected initialization
                            List<Solution> initial = new ArrayList<>();
                            for (int k = 0; k < popSize; k++) {
                                int[] instrumentPartitioning = new int[params.getNumInstr()];
                                int[] orbitAssignment = new int[params.getNumInstr()];

                                // There must be at least one satellite
                                int maxNumSats = PRNG.nextInt(params.getNumInstr()) + 1;

                                for(int j = 0; j < params.getNumInstr(); j++){
                                    instrumentPartitioning[j] = PRNG.nextInt(maxNumSats);
                                }

                                HashMap<Integer, Integer> map = new HashMap<>();
                                int satIndex = 0;
                                for(int m = 0; m < params.getNumInstr(); m++){
                                    int satID = instrumentPartitioning[m];
                                    if(map.keySet().contains(satID)){
                                        instrumentPartitioning[m] = map.get(satID);
                                    }else{
                                        instrumentPartitioning[m] = satIndex;
                                        map.put(satID, satIndex);
                                        satIndex++;
                                    }
                                }
                                Arrays.sort(instrumentPartitioning);

                                int numSats = map.keySet().size();
                                for(int n = 0; n < params.getNumInstr(); n++){
                                    if(n < numSats){
                                        orbitAssignment[n] = PRNG.nextInt(params.getNumOrbits());;
                                    }else{
                                        orbitAssignment[n] = -1;
                                    }
                                }

                                PartitioningArchitecture arch = createPartitioningArchitecture(instrumentPartitioning, orbitAssignment, params);

                                problem.evaluateArch(arch);
                                initial.add(arch);
                            }
                            initialization = new PartitioningInitialization(problem, initial.size(), initial, params);
                            initializationString = "randinjinit";
                            break;

                        default :
                            throw new IllegalStateException("Unrecognized initialization mode");
                    }
                    String runName = "emoea_" + params.getProblemName() + "_" + "partition" + "_" + initializationString + "_" + i;

                    //initialize population structure for algorithm
                    Population population = new Population();
                    EpsilonBoxDominanceArchive archive = new EpsilonBoxDominanceArchive(epsilonDouble);
                    ChainedComparator comp = new ChainedComparator(new ParetoObjectiveComparator());
                    TournamentSelection selection = new TournamentSelection(2, comp);

                    partCross = new PartitioningCrossover(crossoverProbability, params);
                    partFlip = new PartitioningMutation(mutationProbability, params);
                    CompoundVariation variation = new CompoundVariation(partCross, partFlip);

                    Algorithm eMOEA = new EpsilonMOEA(problem, population, archive, selection, variation, initialization);
                    ecs.submit(new PartitioningSearch(eMOEA, properties, savePath, runName, params, evaluationManager));
                }
                for (int i = 0; i < numRuns; ++i) {
                    try {
                        Algorithm alg = ecs.take().get();
                    } catch (InterruptedException | ExecutionException ex) {
                        ex.printStackTrace();
                    }
                }
                pool.shutdown();
                evaluationManager.clear();
                System.out.println("DONE");
                break;

            case RandomPopulation:
                PartitioningProblem problem = new PartitioningProblem(params.getProblemName(), evaluationManager, params, interferingInstrumentsMap, instrumentSynergyMap, dcThreshold, massThreshold, packEffThreshold, numberOfHeuristicObjectives, numberOfHeuristicConstraints, heuristicsConstrained);

                switch (randomMode) {
                    case FullyRandom:
                        System.out.println("Starting random population evaluation for Partitioning Problem");
                        for (int i = 0; i < numRuns; i++) {
                            String runName = "random_" + params.getProblemName() + "_" + "partition" + "_" + i;

                            HashSet<Solution> randomPopulation = new HashSet<>();
                            int numberOfGeneratedArchitectures = 0;
                            while(numberOfGeneratedArchitectures < popSize) {
                                int[] instrumentPartitioning = new int[params.getNumInstr()];
                                int[] orbitAssignment = new int[params.getNumInstr()];

                                // There must be at least one satellite
                                int maxNumSats = PRNG.nextInt(params.getNumInstr()) + 1;

                                for(int j = 0; j < params.getNumInstr(); j++){
                                    instrumentPartitioning[j] = PRNG.nextInt(maxNumSats);
                                }

                                HashMap<Integer, Integer> map = new HashMap<>();
                                int satIndex = 0;
                                for(int m = 0; m < params.getNumInstr(); m++){
                                    int satID = instrumentPartitioning[m];
                                    if(map.keySet().contains(satID)){
                                        instrumentPartitioning[m] = map.get(satID);
                                    }else{
                                        instrumentPartitioning[m] = satIndex;
                                        map.put(satID, satIndex);
                                        satIndex++;
                                    }
                                }
                                Arrays.sort(instrumentPartitioning);

                                int numSats = map.keySet().size();
                                for(int n = 0; n < params.getNumInstr(); n++){
                                    if(n < numSats){
                                        orbitAssignment[n] = PRNG.nextInt(params.getNumOrbits());;
                                    }else{
                                        orbitAssignment[n] = -1;
                                    }
                                }

                                PartitioningArchitecture arch = createPartitioningArchitecture(instrumentPartitioning, orbitAssignment, params);

                                problem.evaluateArch(arch);

                                if (randomPopulation.contains(arch)) {
                                    continue;
                                } else {
                                    randomPopulation.add(arch);
                                    numberOfGeneratedArchitectures++;
                                    System.out.println("Architecture " + numberOfGeneratedArchitectures + " of run " + (i+1) + " computed");
                                }
                            }
                            String filename = savePath + File.separator + runName + ".csv";
                            savePopulationCSV(randomPopulation, filename, params);
                        }
                        pool.shutdown();
                        evaluationManager.clear();
                        System.out.println("DONE");
                        break;
                    case BiasedRandom:
                        System.out.println("Starting biased random population evaluation for Partitioning Problem");
                        for (int i = 0; i < numRuns; i++) {
                            String runName = "biasedrandom_" + params.getProblemName() + "_" + "partition" + "_" + i;

                            HashSet<Solution> biasedRandomPopulation = new HashSet<>();

                            // Add all architectures with a single satellite
                            int[] instrumentPartitioning = new int[params.getNumInstr()];
                            int[] orbitAssigning = new int[params.getNumInstr()];
                            Arrays.fill(orbitAssigning, -1);
                            for (int j = 0; j < params.getNumOrbits(); j++) {
                                orbitAssigning[0] = j;

                                PartitioningArchitecture arch = createPartitioningArchitecture(instrumentPartitioning, orbitAssigning, params);

                                problem.evaluate(arch);
                                biasedRandomPopulation.add(arch);
                            }
                            System.out.println("Single satellite architectures added");

                            // Compute number of architectures with different number of satellites using Poisson Distribution
                            int lambda = 4; // lambda for Poisson's Distribution, biases higher number of architectures for lower number of satellites
                            int[] numberOfSatellites = IntStream.range(2, params.getNumInstr()).toArray();
                            int[] numberOfArchs = new int[numberOfSatellites.length];
                            for (int j = 0; j < numberOfArchs.length; j++) {
                                numberOfArchs[j] = (int) (Math.floor((Math.pow(lambda,numberOfSatellites[j])*Math.exp(-lambda)/CombinatoricsUtils.factorial(numberOfSatellites[j]))*popSize)) +1;
                            }

                            // Add remaining number of architectures with n_orbs number of satellites
                            instrumentPartitioning = IntStream.range(0, params.getNumInstr()).toArray();
                            orbitAssigning = new int[params.getNumInstr()];
                            int numberOfFullArchs = popSize - (params.getNumOrbits() + Arrays.stream(numberOfArchs).sum());
                            int j = 0;
                            while (j < numberOfFullArchs) {
                                Collections.shuffle(Arrays.asList(instrumentPartitioning));
                                for (int p = 0; p < params.getNumInstr(); p++) {
                                    orbitAssigning[p] = PRNG.nextInt(params.getNumOrbits());
                                }

                                PartitioningArchitecture arch = createPartitioningArchitecture(instrumentPartitioning, orbitAssigning, params);

                                if (biasedRandomPopulation.contains(arch)) {
                                    continue;
                                } else {
                                    problem.evaluate(arch);
                                    biasedRandomPopulation.add(arch);
                                    j++;
                                }
                            }
                            System.out.println("Full satellite architectures added");

                            // Populate architectures with different number of satellites
                            instrumentPartitioning = new int[params.getNumInstr()];
                            orbitAssigning = new int[params.getNumInstr()];
                            Arrays.fill(orbitAssigning, -1);
                            for (int k = 0; k < numberOfSatellites.length; k++) {
                                int numSatellites = numberOfSatellites[k];
                                int n = 0;
                                while (n < numberOfArchs[k]) {
                                    for (int m = 0; m < params.getNumInstr(); m++) {
                                        instrumentPartitioning[m] = PRNG.nextInt(numSatellites);
                                    }
                                    for (int m = 0; m < numSatellites; m++) {
                                        orbitAssigning[m] = PRNG.nextInt(params.getNumOrbits());
                                    }

                                    PartitioningArchitecture arch = createPartitioningArchitecture(instrumentPartitioning, orbitAssigning, params);

                                    if (biasedRandomPopulation.contains(arch)) {
                                        continue;
                                    } else {
                                        problem.evaluate(arch);
                                        biasedRandomPopulation.add(arch);
                                        System.out.println("Architecture " + biasedRandomPopulation.size() + " of run " + (i+1) + " computed");
                                        n++;
                                    }
                                }
                            }
                            String filename = savePath + File.separator + runName + ".csv";
                            savePopulationCSV(biasedRandomPopulation, filename, params);
                        }
                        pool.shutdown();
                        evaluationManager.clear();
                        System.out.println("DONE");
                        break;

                    default :
                        throw new IllegalStateException("Unrecognized random run mode");
                }

                //pool.shutdown();
                //evaluationManager.clear();
                //System.out.println("DONE");
                break;

            default :
                throw new IllegalStateException("Unrecognized run mode");
        }
    }

    public static void savePopulationCSV(HashSet<Solution> pop, String filename, ClimateCentricPartitioningParams params) {

        File results = new File(filename);
        results.getParentFile().mkdirs();

        System.out.println("Saving a population in a csv file");

        try (FileWriter writer = new FileWriter(results)) {

            StringJoiner headings = new StringJoiner(",");
            headings.add("Architecture");
            headings.add("Science Score");
            headings.add("Cost");
            headings.add("Duty Cycle Violation");
            headings.add("Instrument Orbit Assignment Violation");
            headings.add("Interference Violation");
            headings.add("Packing Efficiency Violation");
            headings.add("Spacecraft Mass Violation");
            headings.add("Instrument Synergy Violation");
            writer.append(headings.toString());
            writer.append("\n");

            Iterator<Solution> iter = pop.iterator();
            while(iter.hasNext()){

                Solution sol = iter.next();

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

                Architecture arch_abs = new Architecture(instrumentPartitioning, orbitAssignment, 1, params);

                double[] objectives = arch.getObjectives();
                double science = -objectives[0];
                double cost = objectives[1];

                double dutyCycleViolation = (double) arch.getAttribute("DCViolation");
                double instrumentOrbitAssignmentViolation = (double) arch.getAttribute("InstrOrbViolation");
                double interferenceViolation = (double) arch.getAttribute("InterInstrViolation");
                double packingEfficiencyViolation = (double) arch.getAttribute("PackEffViolation");
                double massViolation = (double) arch.getAttribute("SpMassViolation");
                double synergyViolation = (double) arch.getAttribute("SynergyViolation");

                StringJoiner sj = new StringJoiner(",");
                sj.add(arch_abs.toString(" "));
                sj.add(Double.toString(science));
                sj.add(Double.toString(cost));
                sj.add(Double.toString(dutyCycleViolation));
                sj.add(Double.toString(instrumentOrbitAssignmentViolation));
                sj.add(Double.toString(interferenceViolation));
                sj.add(Double.toString(packingEfficiencyViolation));
                sj.add(Double.toString(massViolation));
                sj.add(Double.toString(synergyViolation));

                writer.append(sj.toString());
                writer.append("\n");
            }
            writer.flush();

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static PartitioningArchitecture createPartitioningArchitecture (int[] instrumentPartitions, int[] orbitAssignments, ClimateCentricPartitioningParams params) {
        PartitioningArchitecture arch = new PartitioningArchitecture(params.getNumInstr(), params.getNumOrbits(), 2, params);

        for (int p = 0; p < params.getNumInstr(); p++) {
            IntegerVariable var = new IntegerVariable(instrumentPartitions[p], 0, params.getNumInstr());
            arch.setVariable(p, var);
        }

        for (int q = 0; q < params.getNumInstr(); q++) {
            IntegerVariable var = new IntegerVariable(orbitAssignments[q], -1, params.getNumOrbits());
            arch.setVariable(params.getNumInstr() + q, var);
        }
        return arch;
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

    public enum RunMode{
        RandomPopulation,
        EpsilonMOEA,
    }

    public enum InitializationMode{
        InitializeRandom,
        InitializationRandomAndInjected,
    }

    public enum RandomMode{
        FullyRandom,
        BiasedRandom,
    }
}
