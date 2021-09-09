package seakers.vassarexecheur.search;

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

public class GenerateForMetricsStudyPartitioning {

    public static void main(String[] args) {
        int numRuns = 1;
        int numCpus = 1;

        RunMode runMode  = RunMode.EpsilonMOEA;
        InitializationMode initializationMode = InitializationMode.InitializeRandom;

        ExecutorService pool = Executors.newFixedThreadPool(numCpus);
        CompletionService<Algorithm> ecs = new ExecutorCompletionService<>(pool);

        // Get time
        String timestamp = new SimpleDateFormat("yyyy-MM-dd-HH-mm").format(new Date());

        TypedProperties properties = new TypedProperties();

        int popSize = 300;
        int maxEvals = 5000;
        properties.setInt("maxEvaluations", maxEvals);
        properties.setInt("populationSize", popSize);
        double crossoverProbability = 1.0;
        properties.setDouble("crossoverProbability", crossoverProbability);
        double mutationProbability = 1. / 60.;
        properties.setDouble("mutationProbability", mutationProbability);
        Variation partCross;
        Variation partFlip;
        Initialization initialization;

        // Set seed for random number generator
        //PRNG.setSeed(4321);

        //setup for epsilon MOEA
        double[] epsilonDouble = new double[]{0.001, 10};

        double dcThreshold = 0.5;
        double massThreshold = 3000.0; // [kg]
        double packEffThreshold = 0.4;

        String savePath = System.getProperty("user.dir") + File.separator + "results";

        String resourcesPath = "C:\\SEAK Lab\\SEAK Lab Github\\VASSAR\\VASSAR_resources-heur";

        ClimateCentricPartitioningParams params = new ClimateCentricPartitioningParams(resourcesPath, "CRISP-ATTRIBUTES", "test", "normal");

        AbstractArchitectureEvaluator evaluator = new ArchitectureEvaluator();
        ArchitectureEvaluationManager evaluationManager = new ArchitectureEvaluationManager(params, evaluator);
        evaluationManager.init(numCpus);

        PRNG.setRandom(new SynchronizedMersenneTwister());

        switch (runMode) {

            case EpsilonMOEA:
                System.out.println("Starting EpsilonMOEA for Partitioning Problem");
                String initializationString;

                for (int i = 0; i < numRuns; i++) {

                    PartitioningProblem problem = new PartitioningProblem(params.getProblemName(), evaluationManager, params, dcThreshold, massThreshold, packEffThreshold);

                    switch (initializationMode) {
                        case InitializeRandom:
                            initialization = new PartitioningInitialization(problem, popSize, params);
                            initializationString = "randinit";
                            break;

                        case InitializationRandomAndInjected:
                            // Injected initialization
                            List<Solution> initial = new ArrayList<>();
                            for (int k = 0; k < popSize; k++) {
                                PartitioningArchitecture arch = new PartitioningArchitecture(params.getNumInstr(), params.getNumOrbits(), 2);

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

                                for (int p = 0; p < params.getNumInstr(); p++) {
                                    IntegerVariable var = new IntegerVariable(instrumentPartitioning[p], 0, params.getNumInstr());
                                    arch.setVariable(p, var);
                                }

                                for (int q = 0; q < params.getNumInstr(); q++) {
                                    IntegerVariable var = new IntegerVariable(orbitAssignment[q], -1, numSats);
                                    arch.setVariable(params.getNumInstr() + q, var);
                                }

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
                System.out.println("Starting random population evaluation for Partitioning Problem");
                PartitioningProblem problem = new PartitioningProblem(params.getProblemName(), evaluationManager, params, dcThreshold, massThreshold, packEffThreshold);

                for (int i = 0; i < numRuns; i++) {
                    String runName = "random_" + params.getProblemName() + "_" + "partition" + "_" + i;

                    List<Solution> randomPopulation = new ArrayList<>();
                    for (int k = 0; k < popSize; k++) {
                        PartitioningArchitecture arch = new PartitioningArchitecture(params.getNumInstr(), params.getNumOrbits(), 2);

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

                        for (int p = 0; p < params.getNumInstr(); p++) {
                            IntegerVariable var = new IntegerVariable(instrumentPartitioning[p], 0, params.getNumInstr());
                            arch.setVariable(p, var);
                        }

                        for (int q = 0; q < params.getNumInstr(); q++) {
                            IntegerVariable var = new IntegerVariable(orbitAssignment[q], -1, params.getNumOrbits());
                            arch.setVariable(params.getNumInstr() + q, var);
                        }

                        problem.evaluateArch(arch);
                        randomPopulation.add(arch);
                        System.out.println("Architecture " + (k+1) + " of run " + (i+1) + " computed");
                    }
                    String filename = savePath + File.separator + runName + ".csv";
                    savePopulationCSV(randomPopulation, filename, params);
                }
                pool.shutdown();
                evaluationManager.clear();
                System.out.println("DONE");
                break;

            default :
                throw new IllegalStateException("Unrecognized run mode");
        }
    }

    public static void savePopulationCSV(List<Solution> pop, String filename, ClimateCentricPartitioningParams params) {

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

    public enum RunMode{
        RandomPopulation,
        EpsilonMOEA,
    }

    public enum InitializationMode{
        InitializeRandom,
        InitializationRandomAndInjected,
    }
}
