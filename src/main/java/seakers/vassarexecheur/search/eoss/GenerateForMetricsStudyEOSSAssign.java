package seakers.vassarexecheur.search.eoss;

import eoss.problem.EOSSDatabase;
import eoss.problem.assignment.InstrumentAssignment;
import eoss.problem.assignment.InstrumentAssignmentArchitecture;
import eoss.problem.evaluation.RequirementMode;
import org.moeaframework.algorithm.EpsilonMOEA;
import org.moeaframework.core.*;
import org.moeaframework.core.comparator.ChainedComparator;
import org.moeaframework.core.comparator.ParetoObjectiveComparator;
import org.moeaframework.core.operator.*;
import org.moeaframework.core.operator.binary.BitFlip;
import org.moeaframework.core.variable.BinaryVariable;
import org.moeaframework.util.TypedProperties;
import seakers.orekit.util.OrekitConfig;
import seakers.vassarexecheur.search.intialization.SynchronizedMersenneTwister;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.concurrent.*;

public class GenerateForMetricsStudyEOSSAssign {

    public static void main(String[] args){

        int numRuns = 1;
        int numCpus = 1;

        RunMode runMode  = RunMode.EpsilonMOEA;
        InitializationMode initializationMode = InitializationMode.InitializeRandom;

        String problemString = "ClimateCentric";
        ExecutorService pool = Executors.newFixedThreadPool(numCpus);
        CompletionService<Algorithm> ecs = new ExecutorCompletionService<>(pool);

        // Get time
        String timestamp = new SimpleDateFormat("yyyy-MM-dd-HH-mm").format(new Date());

        //parameters and operators for seakers.vassar_server.search
        TypedProperties properties = new TypedProperties();

        //seakers.vassar_server.search paramaters set here
        int popSize = 300;
        int maxEvals = 3000;
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

        //setup for saving results (not used)
        //properties.setBoolean("saveQuality", false);
        //properties.setBoolean("saveCredits", false);
        //properties.setBoolean("saveSelection", false);

        //initialize problem

        String savePath = System.getProperty("user.dir") + File.separator + "results";

        //initialize EOSS database
        // String problemPath = "C:\\SEAK Lab\\SEAK Lab Github\\EOSS\\EOSS-AIAA\\problems\\climateCentric"; // for lab system
        String problemPath = "C:\\Users\\rosha\\Documents\\SEAK Lab Github\\EOSS\\problems\\climateCentric"; // for laptop
        EOSSDatabase.getInstance();
        EOSSDatabase.loadBuses(new File(problemPath + File.separator + "config" + File.separator + "candidateBuses.xml"));
        EOSSDatabase.loadInstruments(new File(problemPath + File.separator + "xls" + File.separator + "Instrument Capability Definition.xls"));
        EOSSDatabase.loadOrbits(new File(problemPath + File.separator + "config" + File.separator + "candidateOrbits5.xml"));
        EOSSDatabase.loadLaunchVehicles(new File(problemPath + File.separator + "config" + File.separator + "candidateLaunchVehicles.xml"));

        // OrekitConfig.init(numCpus,"C:\\SEAK Lab\\SEAK Lab Github\\EOSS\\EOSS-AIAA"); \\ for lab system
        OrekitConfig.init(numCpus,"C:\\Users\\rosha\\Documents\\SEAK Lab Github\\EOSS"); // for laptop

        PRNG.setRandom(new SynchronizedMersenneTwister());
        switch(runMode) {

            case EpsilonMOEA:
                System.out.println("Starting EpsilonMOEA for Assigning Problem");
                String initializationString;

                for (int i = 0; i < numRuns; ++i) {

                    InstrumentAssignment problem = getAssignmentProblem(problemPath, RequirementMode.FUZZYATTRIBUTE);

                    switch(initializationMode) {
                        case InitializeRandom:
                            initialization = new RandomInitialization(problem, popSize);
                            initializationString = "randinit";
                            break;

                        case InitializationRandomAndInjected:
                            // Injected initialization
                            List<Solution> initial = new ArrayList<>();
                            for(int k = 0; k < popSize; k++){
                                InstrumentAssignmentArchitecture arch = new InstrumentAssignmentArchitecture(new int[]{1}, EOSSDatabase.getNumberOfInstruments(), EOSSDatabase.getNumberOfOrbits(), 2);

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
                                arch.setMissions();
                                problem.evaluateArch(arch);

                                initial.add(arch);
                            }

                            initialization = new InjectedInitialization(problem, initial.size(), initial);
                            initializationString = "randinjinit";
                            break;

                        default :
                            throw new IllegalStateException("Unrecognized initialization mode");

                    }
                    String runName = "emoea_" + problemString + "_" + "assign" + "_" + initializationString + "_" + i;

                    //initialize population structure for algorithm
                    Population population = new Population();
                    EpsilonBoxDominanceArchive archive = new EpsilonBoxDominanceArchive(epsilonDouble);
                    ChainedComparator comp = new ChainedComparator(new ParetoObjectiveComparator());
                    TournamentSelection selection = new TournamentSelection(2, comp);

                    singlecross = new OnePointCrossover(crossoverProbability);
                    bitFlip = new BitFlip(mutationProbability);
                    CompoundVariation variation = new CompoundVariation(singlecross, bitFlip);

                    Algorithm eMOEA = new EpsilonMOEA(problem, population, archive, selection, variation, initialization);
                    ecs.submit(new EOSSSearchAssign(eMOEA, properties, savePath, runName));
                }

                for (int i = 0; i < numRuns; ++i) {
                    try {
                        Algorithm alg = ecs.take().get();
                    } catch (InterruptedException | ExecutionException ex) {
                        ex.printStackTrace();
                    }
                }
                OrekitConfig.end();
                pool.shutdown();
                System.out.println("DONE");
                break;

            case RandomPopulation:

                System.out.println("Starting random population evaluation for Assigning Problem");
                InstrumentAssignment problem = getAssignmentProblem(problemPath, RequirementMode.FUZZYATTRIBUTE);

                for (int i = 0; i < numRuns; ++i) {

                    String runName = "random_" + problemString + "_" + "assign" + "_" + i;

                    //// **1
                    //Initialization randomInitialization = new RandomInitialization(problem, popSize);
                    //List<Solution> randomPopulation = Arrays.asList(randomInitialization.initialize());
                    //List<InstrumentAssignmentArchitecture> randomPopulationArchitectures = new ArrayList<>();

                    //for (Solution sol : randomPopulation) {
                        //InstrumentAssignmentArchitecture arch = new InstrumentAssignmentArchitecture(new int[]{1}, EOSSDatabase.getNumberOfInstruments(), EOSSDatabase.getNumberOfOrbits(), 2);
                        //for (int j = 1; j < arch.getNumberOfVariables(); ++j) {
                            //arch.setVariable(j, sol.getVariable(j));
                        //}

                        //arch.setMissions();
                        //problem.evaluateArch(arch);
                        //randomPopulationArchitectures.add(arch);
                    //}

                    //// **2
                    List<Solution> randomPopulation = new ArrayList<>();
                    for (int k = 0; k < popSize; k++) {
                        InstrumentAssignmentArchitecture arch = new InstrumentAssignmentArchitecture(new int[]{1}, EOSSDatabase.getNumberOfInstruments(), EOSSDatabase.getNumberOfOrbits(), 2);

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

                        arch.setMissions();
                        problem.evaluateArch(arch);

                        randomPopulation.add(arch);
                    }

                    String filename = savePath + File.separator + runName + ".csv";
                    savePopulationCSV(randomPopulation, filename); //// **2
                    //savePopulationCSV2(randomPopulationArchitectures, filename); //// **1
                }
                break;

            default :
                throw new IllegalStateException("Unrecognized run mode");
        }

        OrekitConfig.end();
        pool.shutdown();
        System.out.println("DONE");
    }

    public static void savePopulationCSV2(List<InstrumentAssignmentArchitecture> pop, String filename) {

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

            Iterator<InstrumentAssignmentArchitecture> iter = pop.iterator();
            while(iter.hasNext()){

                InstrumentAssignmentArchitecture arch = iter.next();

                String bitString = "";
                for (int i = 1; i < arch.getNumberOfVariables(); ++i) {
                    bitString += arch.getVariable(i).toString();
                }

                double[] objectives = arch.getObjectives();
                double science = -objectives[0];
                double cost = objectives[1];

                double dutyCycleViolation = (double) arch.getAttribute("dcViolationSum");
                double instrumentOrbitAssignmentViolation = (double) arch.getAttribute("instrumentOrbitAssignmentViolationSum");
                double interferenceViolation = (double) arch.getAttribute("interferenceViolationSum");
                double packingEfficiencyViolation = (double) arch.getAttribute("packingEfficiencyViolationSum");
                double massViolation = (double) arch.getAttribute("massViolationSum");
                double synergyViolation = (double) arch.getAttribute("synergyViolationSum");

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

                InstrumentAssignmentArchitecture arch = (InstrumentAssignmentArchitecture) sol;

                String bitString = "";
                for (int i = 1; i < arch.getNumberOfVariables(); ++i) {
                    bitString += arch.getVariable(i).toString();
                }

                double[] objectives = arch.getObjectives();
                double science = -objectives[0];
                double cost = objectives[1];

                double dutyCycleViolation = (double) arch.getAttribute("dcViolationSum");
                double instrumentOrbitAssignmentViolation = (double) arch.getAttribute("instrumentOrbitAssignmentViolationSum");
                double interferenceViolation = (double) arch.getAttribute("interferenceViolationSum");
                double packingEfficiencyViolation = (double) arch.getAttribute("packingEfficiencyViolationSum");
                double massViolation = (double) arch.getAttribute("massViolationSum");
                double synergyViolation = (double) arch.getAttribute("synergyViolationSum");

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

    private static InstrumentAssignment getAssignmentProblem(String path, RequirementMode mode) {
        return new InstrumentAssignment(path, mode, new int[]{1}, true);
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
