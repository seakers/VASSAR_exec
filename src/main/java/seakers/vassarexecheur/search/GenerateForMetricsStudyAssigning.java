package seakers.vassarexecheur.search;

import org.moeaframework.algorithm.EpsilonMOEA;
import org.moeaframework.core.*;
import org.moeaframework.core.comparator.ChainedComparator;
import org.moeaframework.core.comparator.ParetoObjectiveComparator;
import org.moeaframework.core.operator.*;
import org.moeaframework.core.operator.binary.BitFlip;
import org.moeaframework.core.variable.BinaryVariable;
import org.moeaframework.util.TypedProperties;
import seakers.vassarexecheur.search.intialization.SynchronizedMersenneTwister;
import seakers.vassarexecheur.search.problems.assigning.AssigningArchitecture;
import seakers.vassarexecheur.search.problems.assigning.AssigningProblem;
import seakers.vassarheur.evaluation.AbstractArchitectureEvaluator;
import seakers.vassarheur.evaluation.ArchitectureEvaluationManager;
import seakers.vassarheur.problems.Assigning.ArchitectureEvaluator;
import seakers.vassarheur.problems.Assigning.ClimateCentricAssigningParams;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.concurrent.*;

public class GenerateForMetricsStudyAssigning {

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
        Variation singlecross;
        Variation bitFlip;
        Initialization initialization;

        // Set seed for random number generator
        //PRNG.setSeed(4321);

        //setup for epsilon MOEA
        double[] epsilonDouble = new double[]{0.001, 10};

        double dcThreshold = 0.5;
        double massThreshold = 3000.0; // [kg]
        double packEffThreshold = 0.4;
        boolean considerFeasibility = true;

        String savePath = System.getProperty("user.dir") + File.separator + "results";

        //String resourcesPath = "C:\\SEAK Lab\\SEAK Lab Github\\VASSAR\\VASSAR_resources-heur"; // for lab system
        String resourcesPath = "C:\\Users\\rosha\\Documents\\SEAK Lab Github\\VASSAR\\VASSAR_resources-heur"; // for laptop

        ClimateCentricAssigningParams params = new ClimateCentricAssigningParams(resourcesPath, "CRISP-ATTRIBUTES","test", "normal");

        PRNG.setRandom(new SynchronizedMersenneTwister());

        AbstractArchitectureEvaluator evaluator = new ArchitectureEvaluator(considerFeasibility, dcThreshold, massThreshold, packEffThreshold);
        ArchitectureEvaluationManager evaluationManager = new ArchitectureEvaluationManager(params, evaluator);
        evaluationManager.init(numCpus);

        switch (runMode) {

            case EpsilonMOEA:
                System.out.println("Starting EpsilonMOEA for Assigning Problem");
                String initializationString;

                for (int i = 0; i < numRuns; i++) {

                    AssigningProblem problem = new AssigningProblem(new int[]{1}, params.getProblemName(), evaluationManager, params, dcThreshold, massThreshold, packEffThreshold);

                    switch (initializationMode) {
                        case InitializeRandom:
                            initialization = new RandomInitialization(problem, popSize);
                            initializationString = "randinit";
                            break;

                        case InitializationRandomAndInjected:
                            // Injected initialization
                            List<Solution> initial = new ArrayList<>();
                            for(int k = 0; k < popSize; k++){
                                AssigningArchitecture arch = new AssigningArchitecture(new int[]{1}, params.getNumInstr(), params.getNumOrbits(), 2);

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
                                problem.evaluateArch(arch);

                                initial.add(arch);
                            }

                            initialization = new InjectedInitialization(problem, initial.size(), initial);
                            initializationString = "randinjinit";
                            break;

                        default :
                            throw new IllegalStateException("Unrecognized initialization mode");
                    }

                    String runName = "emoea_" + params.getProblemName() + "_" + "assign" + "_" + initializationString + "_" + i;

                    //initialize population structure for algorithm
                    Population population = new Population();
                    EpsilonBoxDominanceArchive archive = new EpsilonBoxDominanceArchive(epsilonDouble);
                    ChainedComparator comp = new ChainedComparator(new ParetoObjectiveComparator());
                    TournamentSelection selection = new TournamentSelection(2, comp);

                    singlecross = new OnePointCrossover(crossoverProbability);
                    bitFlip = new BitFlip(mutationProbability);
                    CompoundVariation variation = new CompoundVariation(singlecross, bitFlip);

                    Algorithm eMOEA = new EpsilonMOEA(problem, population, archive, selection, variation, initialization);
                    ecs.submit(new AssigningSearch(eMOEA, properties, savePath, runName, evaluationManager));
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
                System.out.println("Starting random population evaluation for Assigning Problem");

                for (int i = 0; i < numRuns; i++) {
                    AssigningProblem problem = new AssigningProblem(new int[]{1}, params.getProblemName(), evaluationManager, params, dcThreshold, massThreshold, packEffThreshold);
                    String runName = "random_" + params.getProblemName() + "_" + "assign" + "_" + i;

                    List<Solution> randomPopulation = new ArrayList<>();
                    for (int k = 0; k < popSize; k++) {
                        AssigningArchitecture arch = new AssigningArchitecture(new int[]{1}, params.getNumInstr(), params.getNumOrbits(), 2);

                        StringBuilder bitStringBuilder = new StringBuilder(60);
                        for (int j = 1; j < arch.getNumberOfVariables(); ++j) {
                            BinaryVariable var = new BinaryVariable(1);
                            if (PRNG.nextDouble() < 1. / 2.) {
                                var.set(0, true);
                                bitStringBuilder.append("1");
                            } else {
                                var.set(0, false);
                                bitStringBuilder.append("0");
                            }
                            arch.setVariable(j, var);
                        }
                        problem.evaluateArch(arch);
                        randomPopulation.add(arch);
                        System.out.println("Architecture " + (k+1) + " of run " + (i+1) + " computed");
                    }
                    String filename = savePath + File.separator + runName + ".csv";
                    savePopulationCSV(randomPopulation, filename);
                }
                pool.shutdown();
                evaluationManager.clear();
                System.out.println("DONE");
                break;
            default :
                throw new IllegalStateException("Unrecognized run mode");
        }
    }

    public static void savePopulationCSV(List<Solution> pop, String filename) {

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

                AssigningArchitecture arch = (AssigningArchitecture) sol;

                String bitString = "";
                for (int i = 1; i < arch.getNumberOfVariables(); ++i) {
                    bitString += arch.getVariable(i).toString();
                }

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
                sj.add(bitString);
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
