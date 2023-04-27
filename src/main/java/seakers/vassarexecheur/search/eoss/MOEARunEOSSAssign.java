package seakers.vassarexecheur.search.eoss;

import eoss.problem.EOSSDatabase;
import org.moeaframework.core.Algorithm;
import org.moeaframework.core.Initialization;
import org.moeaframework.core.PRNG;
import org.moeaframework.core.Variation;
import org.moeaframework.util.TypedProperties;
import seakers.orekit.util.OrekitConfig;
import seakers.vassarheur.evaluation.AbstractArchitectureEvaluator;
import seakers.vassarheur.evaluation.ArchitectureEvaluationManager;
import seakers.vassarheur.problems.Assigning.ArchitectureEvaluator;
import seakers.vassarheur.problems.Assigning.ClimateCentricAssigningParams;
import seakers.vassarexecheur.search.intialization.SynchronizedMersenneTwister;

import java.io.File;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public class MOEARunEOSSAssign {

    public static void main(String[] args) {

        int numRuns = 2;
        int numCpus = 1;

        //// Select enforcement methods for different heuristics
        // boolean_array = [int_pen, AOS, bias_init, ACH, objective, constraint]
        boolean[] dutyCycleEnforcement = {false, false, false, false, false, false};
        boolean[] instrumentOrbitEnforcement = {false, false, false, false, false, false};
        boolean[] instrumentInterferenceEnforcement = {false, false, false, false, false, false};
        boolean[] packingEfficiencyEnforcement = {false, false, false, false, false, false};
        boolean[] spacecraftMassEnforcement = {false, false, false, false, false, false};
        boolean[] instrumentSynergyEnforcement = {false, false, false, false, false, false};

        InitializationMode initializationMode = InitializationMode.InitializeRandom;

        String problem = "ClimateCentric";
        ExecutorService pool = Executors.newFixedThreadPool(numCpus);
        CompletionService<Algorithm> ecs = new ExecutorCompletionService<>(pool);

        //parameters and operators for seakers.vassar_server.search
        TypedProperties properties = new TypedProperties();

        //seakers.vassar_server.search paramaters set here
        int popSize = 1000;
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

        //setup for epsilon MOEA
        double[] epsilonDouble = new double[]{0.001, 10};

        //initialize problem
        String resourcesPath = "C:\\SEAK Lab\\SEAK Lab Github\\VASSAR\\VASSAR_resources-master";
        String savePath = System.getProperty("user.dir") + File.separator + "results";

        ClimateCentricAssigningParams params = new ClimateCentricAssigningParams(resourcesPath, "CRISP-ATTRIBUTES","test", "normal");
        AbstractArchitectureEvaluator evaluator = new ArchitectureEvaluator();
        ArchitectureEvaluationManager evaluationManager = new ArchitectureEvaluationManager(params, evaluator);
        evaluationManager.init(numCpus);
        OrekitConfig.init(numCpus, params.orekitResourcesPath);

        //initialize EOSS database
        String problemPath = "C:\\SEAK Lab\\SEAK Lab Github\\EOSS\\EOSS-dev\\problems\\climateCentric";
        EOSSDatabase.getInstance();
        EOSSDatabase.loadBuses(new File(problemPath + File.separator + "config" + File.separator + "candidateBuses.xml"));
        EOSSDatabase.loadInstruments(new File(problemPath + File.separator + "xls" + File.separator + "Instrument Capability Definition.xls"));
        EOSSDatabase.loadOrbits(new File(problemPath + File.separator + "config" + File.separator + "candidateOrbits5.xml"));
        EOSSDatabase.loadLaunchVehicles(new File(problemPath + File.separator + "config" + File.separator + "candidateLaunchVehicles.xml"));

        PRNG.setRandom(new SynchronizedMersenneTwister());

        StringBuilder fileSaveNameConstraint = new StringBuilder();



    }

    public enum InitializationMode{
        InitializeRandom,
        InitializationRandomAndInjected,
    }
}
