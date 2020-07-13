package seakers.vassartest;

import org.moeaframework.algorithm.AbstractEvolutionaryAlgorithm;
import org.moeaframework.algorithm.EpsilonMOEA;
import org.moeaframework.core.*;
import org.moeaframework.core.comparator.ChainedComparator;
import org.moeaframework.core.comparator.ParetoObjectiveComparator;
import org.moeaframework.core.operator.*;
import org.moeaframework.core.operator.binary.BitFlip;
import org.moeaframework.core.variable.BinaryVariable;
import org.moeaframework.util.TypedProperties;
import seakers.orekit.util.OrekitConfig;
import seakers.vassar.Result;
import seakers.vassar.architecture.AbstractArchitecture;
import seakers.vassar.evaluation.AbstractArchitectureEvaluator;
import seakers.vassar.evaluation.ArchitectureEvaluationManager;
import seakers.vassar.problems.Assigning.ArchitectureEvaluator;
import seakers.vassar.problems.Assigning.ClimateCentricParams;
import seakers.vassartest.search.ClimateCentricProblemSearch;
import seakers.vassartest.search.problems.Assigning.AssigningArchitecture;
import seakers.vassartest.search.problems.Assigning.AssigningProblem;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.concurrent.*;


public class ClimateCentricEvaluation2019Spring {

    public static void main(String[] args){

        int numRuns = 1;
        int numCpus = 2;

        RunMode mode  = RunMode.EpsilonMOEA;

        String problem = "ClimateCentric";
        ExecutorService pool = Executors.newFixedThreadPool(numCpus);
        CompletionService<Algorithm> ecs = new ExecutorCompletionService<>(pool);

        // Get time
        String timestamp = new SimpleDateFormat("yyyy-MM-dd-HH-mm").format(new Date());

        //parameters and operators for seakers.vassar_server.search
        TypedProperties properties = new TypedProperties();

        //seakers.vassar_server.search paramaters set here
        int popSize = 500;
        int maxEvals = 10000;
        properties.setInt("maxEvaluations", maxEvals);
        properties.setInt("populationSize", popSize);
        double crossoverProbability = 1.0;
        properties.setDouble("crossoverProbability", crossoverProbability);
        double mutationProbability = 1. / 60.;
        properties.setDouble("mutationProbability", mutationProbability);
        Variation singlecross;
        Variation bitFlip;
        Initialization initialization;

        // Set seed for random number generator
        PRNG.setSeed(4321);

        //setup for epsilon MOEA
        double[] epsilonDouble = new double[]{0.001, 10};

        //setup for saving results
        properties.setBoolean("saveQuality", false);
        properties.setBoolean("saveCredits", false);
        properties.setBoolean("saveSelection", false);

        //initialize problem
        String resourcesPath = "../VASSAR_resources";
        String savePath = resourcesPath + File.separator + "results";

        ClimateCentricParams params = new ClimateCentricParams(resourcesPath, "CRISP-ATTRIBUTES",
                "test", "normal");
        AbstractArchitectureEvaluator evaluator = new ArchitectureEvaluator();
        ArchitectureEvaluationManager evaluationManager = new ArchitectureEvaluationManager(params, evaluator);
        evaluationManager.init(numCpus);
        OrekitConfig.init(numCpus, params.orekitResourcesPath);


        switch(mode) {

            case EpsilonMOEA:

                System.out.println("Starting EpsilonMOEA for binary input data");

                for (int i = 0; i < numRuns; ++i) {

                    String runName = "emoea_" + problem + "_" + timestamp + "_" + i;
                    Problem assigningProblem = new AssigningProblem(new int[]{1}, "ClimateCentric", evaluationManager, params)
                            .setSaveEvaluatedSolutions();

                    ////////////////////////////////////////////
                    // Injected initialization
                    List<Solution> initial = new ArrayList<>();
                    for(int k = 0; k < popSize; k++){
                        AssigningArchitecture arch = new AssigningArchitecture(new int[]{1},
                                params.getNumInstr(), params.getNumOrbits(), 2);

                        StringBuilder bitStringBuilder = new StringBuilder(60);
                        for (int j = 1; j < arch.getNumberOfVariables(); ++j) {
                            BinaryVariable var = new BinaryVariable(1);
                            if(PRNG.nextDouble() < 1./6.){
                                var.set(0, true);
                                bitStringBuilder.append("1");
                            }else{
                                var.set(0, false);
                                bitStringBuilder.append("0");
                            }
                            arch.setVariable(j, var);
                        }

                        try {
                            AbstractArchitecture arch_old = new seakers.vassar.problems.Assigning.Architecture(bitStringBuilder.toString(), 1, params);
                            Result result = evaluationManager.evaluateArchitectureSync(arch_old, "Slow");
                            System.out.println(k + ": " + String.format("Arch %s Science = %10f; Cost = %10f",
                                    arch.toString(), result.getScience(), result.getCost()));
                            arch.setObjective(0, -result.getScience()); //negative because MOEAFramework assumes minimization problems
                            arch.setObjective(1, result.getCost()); //normalize cost to maximum value
                            arch.setAlreadyEvaluated(true);
                        }
                        catch (Exception e) {
                            e.printStackTrace();
                        }
                        initial.add(arch);
                    }
                    initialization = new InjectedInitialization(assigningProblem, initial.size(), initial);
                    ////////////////////////////////////////////

//            initialization = new RandomInitialization(assignmentProblem, popSize);

                    //initialize population structure for algorithm
                    Population population = new Population();
                    EpsilonBoxDominanceArchive archive = new EpsilonBoxDominanceArchive(epsilonDouble);
                    ChainedComparator comp = new ChainedComparator(new ParetoObjectiveComparator());
                    TournamentSelection selection = new TournamentSelection(2, comp);

                    singlecross = new OnePointCrossover(crossoverProbability);
                    bitFlip = new BitFlip(mutationProbability);
                    CompoundVariation variation = new CompoundVariation(singlecross, bitFlip);

                    Algorithm eMOEA = new EpsilonMOEA(assigningProblem, population, archive, selection, variation, initialization);
                    ecs.submit(new ClimateCentricProblemSearch(eMOEA, properties, savePath, runName));
                }

                for (int i = 0; i < numRuns; ++i) {
                    try {
                        Algorithm alg = ecs.take().get();
                    } catch (InterruptedException | ExecutionException ex) {
                        ex.printStackTrace();
                    }
                }
                break;

            case RandomPopulation:

                System.out.println("Starting random population evaluation for binary input data");

                for (int i = 0; i < numRuns; ++i) {

                    String runName = "random_" + problem + "_" + timestamp + "_" + i;

                    List<Solution> randomPopulation = new ArrayList<>();
                    for(int k = 0; k < popSize; k++){
                        AssigningArchitecture arch = new AssigningArchitecture(new int[]{1},
                                params.getNumInstr(), params.getNumOrbits(), 2);

                        StringBuilder bitStringBuilder = new StringBuilder(60);
                        for (int j = 1; j < arch.getNumberOfVariables(); ++j) {
                            BinaryVariable var = new BinaryVariable(1);
                            if(PRNG.nextDouble() < 1./3.){
                                var.set(0, true);
                                bitStringBuilder.append("1");
                            }else{
                                var.set(0, false);
                                bitStringBuilder.append("0");
                            }
                            arch.setVariable(j, var);
                        }

                        try {
                            AbstractArchitecture arch_old = new seakers.vassar.problems.Assigning.Architecture(bitStringBuilder.toString(), 1, params);
                            Result result = evaluationManager.evaluateArchitectureSync(arch_old, "Slow");
                            System.out.println(k + ": " + String.format("Arch %s Science = %10f; Cost = %10f",
                                    arch.toString(), result.getScience(), result.getCost()));
                            arch.setObjective(0, -result.getScience()); //negative because MOEAFramework assumes minimization problems
                            arch.setObjective(1, result.getCost()); //normalize cost to maximum value
                            arch.setAlreadyEvaluated(true);
                            randomPopulation.add(arch);
                        }
                        catch (Exception e) {
                            e.printStackTrace();
                        }
                    }

                    String filename = savePath + File.separator +  runName;
                    savePopulationCSV(randomPopulation, filename);
                }
                break;

            default :
                throw new IllegalStateException("Unrecognized mode");
        }

        OrekitConfig.end();
        evaluationManager.clear();
        pool.shutdown();
        System.out.println("DONE");
    }

    public static void savePopulationCSV(List<Solution> pop, String filename) {

        File results = new File(filename);
        results.getParentFile().mkdirs();

        System.out.println("Saving a population in a csv file");

        try (FileWriter writer = new FileWriter(results)) {

            Iterator<Solution> iter = pop.iterator();
            while(iter.hasNext()){

                Solution sol = iter.next();

                AssigningArchitecture arch = (AssigningArchitecture) sol;

                String bitString = "";
                for (int i = 1; i < arch.getNumberOfVariables(); ++i) {
                    bitString += arch.getVariable(i).toString();
                }

                double[] objectives = arch.getObjectives();
                double science = -objectives[0];
                double cost = objectives[1];

                StringJoiner sj = new StringJoiner(",");
                sj.add(bitString);
                sj.add(Double.toString(science));
                sj.add(Double.toString(cost));
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
}
