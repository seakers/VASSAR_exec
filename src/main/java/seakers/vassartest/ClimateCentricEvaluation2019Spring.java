package seakers.vassartest;

import org.moeaframework.algorithm.EpsilonMOEA;
import org.moeaframework.core.*;
import org.moeaframework.core.comparator.ChainedComparator;
import org.moeaframework.core.comparator.ParetoObjectiveComparator;
import org.moeaframework.core.operator.*;
import org.moeaframework.core.operator.binary.BitFlip;
import org.moeaframework.util.TypedProperties;
import seakers.orekit.util.OrekitConfig;
import seakers.vassar.evaluation.AbstractArchitectureEvaluator;
import seakers.vassar.evaluation.ArchitectureEvaluationManager;
import seakers.vassar.problems.Assigning.ArchitectureEvaluator;
import seakers.vassar.problems.Assigning.ClimateCentricParams;
import seakers.vassartest.search.ClimateCentricProblemSearch;
import seakers.vassartest.search.problems.Assigning.AssigningProblem;

import java.io.File;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.concurrent.*;


public class ClimateCentricEvaluation2019Spring {

    public static void main(String[] args){

        System.out.println("Starting GA for binary input data");

        int numRuns = 1;
        int numCpus = 1;

        String problem = "ClimateCentric";
        ExecutorService pool = Executors.newFixedThreadPool(numCpus);
        CompletionService<Algorithm> ecs = new ExecutorCompletionService<>(pool);

        // Get time
        String timestamp = new SimpleDateFormat("yyyy-MM-dd-HH-mm").format(new Date());

        //parameters and operators for seakers.vassar_server.search
        TypedProperties properties = new TypedProperties();

        //seakers.vassar_server.search paramaters set here
        int popSize = 500;
        int maxEvals = 7000;
        properties.setInt("maxEvaluations", maxEvals);
        properties.setInt("populationSize", popSize);
        double crossoverProbability = 1.0;
        properties.setDouble("crossoverProbability", crossoverProbability);
        double mutationProbability = 1. / 60.;
        properties.setDouble("mutationProbability", mutationProbability);
        Variation singlecross;
        Variation bitFlip;
        Initialization initialization;

        //setup for epsilon MOEA
        double[] epsilonDouble = new double[]{0.001, 10};

        //setup for saving results
        properties.setBoolean("saveQuality", false);
        properties.setBoolean("saveCredits", false);
        properties.setBoolean("saveSelection", false);

        //initialize problem
        String resourcesPath = "../VASSAR_resources";
        ClimateCentricParams params = new ClimateCentricParams(resourcesPath, "CRISP-ATTRIBUTES",
                "test", "normal");
        AbstractArchitectureEvaluator evaluator = new ArchitectureEvaluator();
        ArchitectureEvaluationManager AEM = new ArchitectureEvaluationManager(params, evaluator);
        AEM.init(numCpus);
        OrekitConfig.init(numCpus, params.orekitResourcesPath);

        for (int i = 0; i < numRuns; ++i) {

            Problem assignmentProblem = new AssigningProblem(new int[]{1}, "ClimateCentric", AEM, params);
            initialization = new RandomInitialization(assignmentProblem, popSize);

            //initialize population structure for algorithm
            Population population = new Population();
            EpsilonBoxDominanceArchive archive = new EpsilonBoxDominanceArchive(epsilonDouble);
            ChainedComparator comp = new ChainedComparator(new ParetoObjectiveComparator());
            TournamentSelection selection = new TournamentSelection(2, comp);

            singlecross = new OnePointCrossover(crossoverProbability);
            bitFlip = new BitFlip(mutationProbability);
            CompoundVariation variation = new CompoundVariation(singlecross, bitFlip);

            Algorithm eMOEA = new EpsilonMOEA(assignmentProblem, population, archive, selection, variation, initialization);

            String runName = "emoea_" + problem + "_" + timestamp + "_" + i;
            ecs.submit(new ClimateCentricProblemSearch(eMOEA, properties, resourcesPath + File.separator + "results", runName));
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
