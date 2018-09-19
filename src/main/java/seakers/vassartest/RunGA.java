package seakers.vassartest;

import org.moeaframework.algorithm.EpsilonMOEA;
import org.moeaframework.core.*;
import org.moeaframework.core.comparator.ChainedComparator;
import org.moeaframework.core.comparator.ParetoObjectiveComparator;
import org.moeaframework.core.operator.*;
import org.moeaframework.core.operator.binary.BitFlip;
import org.moeaframework.core.variable.BinaryVariable;
import org.moeaframework.util.TypedProperties;
import seakers.vassar.evaluation.AbstractArchitectureEvaluator;
import seakers.vassar.evaluation.ArchitectureEvaluationManager;
import seakers.vassartest.search.problems.Assigning.AssigningArchitecture;
import seakers.vassartest.search.problems.Assigning.AssigningProblem;
import seakers.vassartest.search.TimedSearch;
import seakers.vassar.problems.Assigning.ClimateCentricParams;
import seakers.architecture.operators.IntegerUM;
import seakers.orekit.util.OrekitConfig;
import seakers.vassar.problems.Assigning.ArchitectureEvaluator;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.*;

public class RunGA {
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        System.out.println("Starting GA for binary input data");

        int numRuns = 15;
        int numCpus = 5;
        int startOn = 15;
        String problem = "weather";
        String problemCap = "Weather";

        ExecutorService pool = Executors.newFixedThreadPool(numCpus);
        CompletionService<Algorithm> ecs = new ExecutorCompletionService<>(pool);

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
        Variation intergerMutation;
        Initialization initialization;

        //setup for epsilon MOEA
        double[] epsilonDouble = new double[]{0.001, 1};

        //setup for saving results
        properties.setBoolean("saveQuality", true);
        properties.setBoolean("saveCredits", true);
        properties.setBoolean("saveSelection", true);

        //initialize problem
        String path = "../VASSAR_resources";
        ClimateCentricParams params = new ClimateCentricParams(path, "CRISP-ATTRIBUTES",
                "test", "normal", "search_heuristic_rules_smap_127");
        params.aggregationXls = params.problemPath + "/xls/Aggregation Rules-" + problemCap + ".xls";
        AbstractArchitectureEvaluator evaluator = new ArchitectureEvaluator(params);
        ArchitectureEvaluationManager AEM = new ArchitectureEvaluationManager(params, evaluator);
        AEM.init(numCpus);
        OrekitConfig.init(numCpus, params.orekitResourcesPath);

        for (int i = 0; i < numRuns; ++i) {

            Problem assignmentProblem = new AssigningProblem(new int[]{1}, "ClimateCentric", AEM, params);

            // Create a solution for each input arch in the dataset
            String csvFile = params.pathSaveResults + "/start_" + problem + ".csv";
            String line = "";
            String cvsSplitBy = ",";

            List<Solution> initial = new ArrayList<>();
            try (BufferedReader br = new BufferedReader(new FileReader(csvFile))) {
                while ((line = br.readLine()) != null) {
                    // use comma as separator
                    String[] csvArch = line.split(cvsSplitBy);
                    AssigningArchitecture arch = new AssigningArchitecture(new int[]{1},
                            params.getNumInstr(), params.getNumOrbits(), 2);
                    for (int j = 1; j < arch.getNumberOfVariables(); ++j) {
                        BinaryVariable var = new BinaryVariable(1);
                        var.set(0, csvArch[0].charAt(j-1) == '1');
                        arch.setVariable(j, var);
                    }
                    arch.setObjective(0, -Double.valueOf(csvArch[1]));
                    arch.setObjective(1, Double.valueOf(csvArch[2]));
                    arch.setAlreadyEvaluated(true);
                    initial.add(arch);
                }
            }
            catch (IOException e) {
                e.printStackTrace();
            }

            initialization = new InjectedInitialization(assignmentProblem, initial.size(), initial);

            //initialize population structure for algorithm
            Population population = new Population();
            EpsilonBoxDominanceArchive archive = new EpsilonBoxDominanceArchive(epsilonDouble);
            ChainedComparator comp = new ChainedComparator(new ParetoObjectiveComparator());
            TournamentSelection selection = new TournamentSelection(2, comp);

            singlecross = new OnePointCrossover(crossoverProbability);
            bitFlip = new BitFlip(mutationProbability);
            intergerMutation = new IntegerUM(mutationProbability);
            CompoundVariation var = new CompoundVariation(singlecross, bitFlip, intergerMutation);

            Algorithm eMOEA = new EpsilonMOEA(assignmentProblem, population, archive, selection, var, initialization);
            ecs.submit(new TimedSearch(eMOEA, properties, path + File.separator + "results", "emoea_" + problem + (i+startOn)));
        }

        for (int i = 0; i < numRuns; ++i) {
            try {
                Algorithm alg = ecs.take().get();
            } catch (InterruptedException | ExecutionException ex) {
                ex.printStackTrace();
            }
        }

        OrekitConfig.end();
        AEM.clear();
        pool.shutdown();

        System.out.println("DONE");
    }
}
