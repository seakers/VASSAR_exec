package seakers.vassarexecheur.search;

import org.moeaframework.core.*;
import org.moeaframework.core.comparator.ChainedComparator;
import org.moeaframework.core.comparator.DominanceComparator;
import org.moeaframework.core.comparator.ParetoObjectiveComparator;
import org.moeaframework.core.operator.*;
import org.moeaframework.core.operator.binary.BitFlip;
import org.moeaframework.problem.AbstractProblem;
import org.moeaframework.util.TypedProperties;
import seakers.aos.aos.AOSMOEA;
import seakers.aos.creditassignment.offspringparent.OffspringParentDomination;
import seakers.aos.creditassignment.setimprovement.SetImprovementDominance;
import seakers.aos.operator.AOSVariation;
import seakers.aos.operator.AOSVariationOP;
import java.io.File;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.concurrent.*;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.moeaframework.algorithm.EpsilonMOEA;
import seakers.aos.operator.AOSVariationSI;
import seakers.aos.operatorselectors.AdaptivePursuit;
import seakers.aos.operatorselectors.OperatorSelector;
import seakers.aos.operatorselectors.ProbabilityMatching;
import seakers.vassarexecheur.search.constrainthandling.KnowledgeStochasticRanking;
import seakers.vassarexecheur.search.intialization.SynchronizedMersenneTwister;
import seakers.vassarexecheur.search.operators.assigning.*;
import seakers.vassarexecheur.search.operators.partitioning.*;
import seakers.vassarexecheur.search.problems.assigning.AssigningProblem;
import seakers.vassarexecheur.search.problems.partitioning.PartitioningProblem;
import seakers.vassarheur.BaseParams;
import seakers.vassarheur.evaluation.AbstractArchitectureEvaluator;
import seakers.vassarheur.evaluation.ArchitectureEvaluationManager;
import seakers.vassarheur.problems.Assigning.ArchitectureEvaluator;
import seakers.vassarheur.problems.Assigning.ClimateCentricAssigningParams;
import seakers.vassarheur.problems.PartitioningAndAssigning.ClimateCentricPartitioningParams;

/**
 * Executable class for different eMOEA run experiments for the EOSS Optimization problems.
 *
 * @author roshan94
 */

public class MOEARun {
    /**
     * pool of resources
     */
    private static ExecutorService pool;

    /**
     * Executor completion services helps remove completed tasks
     */
    private static CompletionService<Algorithm> ecs;

    public static void main (String[] args) {

        // Define problem parameters
        String csvPath = System.getProperty("user.dir");

        boolean assigningProblem = true; // True -> assigning problem, False -> partitioning problem

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

        int numCPU = 4;
        int numRuns = 4;
        pool = Executors.newFixedThreadPool(numCPU);
        ecs = new ExecutorCompletionService<>(pool);

        //setup for epsilon MOEA
        double[] epsilonDouble = new double[]{0.001, 10};

        double dcThreshold = 0.5;
        double massThreshold = 3000.0; // [kg]
        double packEffThreshold = 0.4;
        boolean considerFeasibility = true;

        // Get time
        String timestamp = new SimpleDateFormat("yyyy-MM-dd-HH-mm").format(new Date());

        TypedProperties properties = new TypedProperties();

        int popSize = 300;
        int maxEvals = 5000;
        properties.setInt("maxEvaluations", maxEvals);
        properties.setInt("populationSize", popSize);
        double crossoverProbability = 1.0;
        properties.setDouble("crossoverProbability", crossoverProbability);

        String resourcesPath = "C:\\SEAK Lab\\SEAK Lab Github\\VASSAR\\VASSAR_resources-heur";
        String savePath = System.getProperty("user.dir") + File.separator + "results";

        double mutationProbability;
        BaseParams params;
        AbstractArchitectureEvaluator evaluator;
        if (assigningProblem) {
            mutationProbability = 1. / 60.;
            params = new ClimateCentricAssigningParams(resourcesPath, "CRISP-ATTRIBUTES","test", "normal");
            evaluator = new ArchitectureEvaluator(considerFeasibility, dcThreshold, massThreshold, packEffThreshold);
        } else {
            mutationProbability = 1. / 24.; // Based on the 12 instruments for the ClimateCentric Problem
            params = new ClimateCentricPartitioningParams(resourcesPath, "CRISP-ATTRIBUTES", "test", "normal");
            evaluator = new seakers.vassarheur.problems.PartitioningAndAssigning.ArchitectureEvaluator(considerFeasibility, dcThreshold, massThreshold, packEffThreshold);
        }
        ArchitectureEvaluationManager evaluationManager = new ArchitectureEvaluationManager(params, evaluator);
        evaluationManager.init(numCPU);
        properties.setDouble("mutationProbability", mutationProbability);

        Initialization initialization = null;

        String fileSaveNameProblem  = "";
        if (assigningProblem) {
            fileSaveNameProblem = "_assigning";
        } else {
            fileSaveNameProblem = "_partitioning";
        }

        String[] heuristicAbbreviations = {"d","o","i","p","m","s"};
        StringBuilder fileSaveNameConstraint = new StringBuilder();
        for (int i = 0; i < heuristicsConstrained[0].length; i++) {
            StringBuilder enforcedHeuristics = new StringBuilder();
            int heuristicCount = 0;
            for (int j = 0; j < heuristicsConstrained.length; j++) {
                if (heuristicsConstrained[j][i]) {
                    enforcedHeuristics.append(heuristicAbbreviations[j]);
                    heuristicCount++;
                }
            }
            if (heuristicCount > 0) {
                fileSaveNameConstraint.append(enforcedHeuristics.toString()).append("con").append(Integer.toString(i)).append("_");
            }
        }

        PRNG.setRandom(new SynchronizedMersenneTwister());

        for (int i = 0; i < numRuns; i++) {
            // Problem class
            AbstractProblem satelliteProblem;
            if (assigningProblem) {
                satelliteProblem = new AssigningProblem(new int[]{1}, params.getProblemName(), evaluationManager, params, dcThreshold, massThreshold, packEffThreshold, numberOfHeuristicObjectives, numberOfHeuristicConstraints, heuristicsConstrained);
            } else {
                satelliteProblem = new PartitioningProblem(params.getProblemName(), evaluationManager, params, dcThreshold, massThreshold, packEffThreshold, numberOfHeuristicObjectives, numberOfHeuristicConstraints, heuristicsConstrained);
            }

            // Initial population
            if (dutyCycleConstrained[2] || instrumentOrbitRelationsConstrained[2] || interferenceConstrained[2] || packingEfficiencyConstrained[2] || spacecraftMassConstrained[2] || synergyConstrained[2]) {
                System.out.println("Biased Initialization not programmed yet");
            } else {
                initialization = new RandomInitialization(satelliteProblem, popSize);
            }

            // Initialize the population
            Population population = new Population();

            EpsilonBoxDominanceArchive archive = new EpsilonBoxDominanceArchive(epsilonDouble);

            Algorithm moeaObj;
            DominanceComparator comp = null;
            TournamentSelection selection;
            CompoundVariation var = null;
            AOSVariation aosStrategy = null;

            // Initialize heuristic operators
            Variation repairDutyCycle;
            Variation repairInstrumentOrbitRelations;
            Variation repairInterference;
            Variation repairPackingEfficiency;
            Variation repairMass;
            Variation repairSynergy;

            if (assigningProblem) {
                repairDutyCycle = new RepairDutyCycleAssigning(dcThreshold, 1, params);
                repairInstrumentOrbitRelations = new RepairInstrumentOrbitAssigning(1, evaluationManager.getResourcePool(), (ArchitectureEvaluator) evaluator, params);
                repairInterference = new RepairInterferenceAssigning(1, evaluationManager.getResourcePool(), (ArchitectureEvaluator) evaluator, params);
                repairPackingEfficiency = new RepairPackingEfficiencyAssigning(packEffThreshold, 1, params);
                repairMass = new RepairMassAssigning(massThreshold, 1, params);
                repairSynergy = new RepairSynergyAssigning(1, evaluationManager.getResourcePool(), (ArchitectureEvaluator) evaluator, params);
            } else {
                repairDutyCycle = new RepairDutyCyclePartitioning(dcThreshold, 1, params);
                repairInstrumentOrbitRelations = new RepairInstrumentOrbitPartitioning(1, evaluationManager.getResourcePool(), (seakers.vassarheur.problems.PartitioningAndAssigning.ArchitectureEvaluator) evaluator, params);
                repairInterference = new RepairInterferencePartitioning(1, evaluationManager.getResourcePool(), (seakers.vassarheur.problems.PartitioningAndAssigning.ArchitectureEvaluator) evaluator, params);
                repairPackingEfficiency = new RepairPackingEfficiencyAssigning(packEffThreshold, 1, params);
                repairMass = new RepairMassPartitioning(massThreshold, 1, params);
                repairSynergy = new RepairSynergyPartitioning(1, evaluationManager.getResourcePool(), (seakers.vassarheur.problems.PartitioningAndAssigning.ArchitectureEvaluator) evaluator, params);
            }

            Variation[] heuristicOperators ={repairDutyCycle, repairInstrumentOrbitRelations, repairInterference, repairPackingEfficiency, repairMass, repairSynergy};
            String[] heuristicAttributes = {"DCViolation","InstrOrbViolation","InterInstrViolation","PackEffViolation","SpMassViolation","SynergyViolation"};

            if (dutyCycleConstrained[3] || instrumentOrbitRelationsConstrained[3] || interferenceConstrained[3] || packingEfficiencyConstrained[3] || spacecraftMassConstrained[3] || synergyConstrained[3]) { // Adaptive Constraint Handling objects

                KnowledgeStochasticRanking ksr;

                HashMap<Variation, String> constraintOperatorMap = new HashMap<>();

                for (int j = 0; j < heuristicsConstrained.length; j++) {
                    if (heuristicsConstrained[j][3]) {
                        constraintOperatorMap.put(heuristicOperators[j], heuristicAttributes[j]);
                    }
                }
                ksr = new KnowledgeStochasticRanking(constraintOperatorMap.size(), constraintOperatorMap.values(), archive);
                comp = new ChainedComparator(ksr, new ParetoObjectiveComparator());
            } else { // Epsilon MOEA objects

                comp = new ChainedComparator(new ParetoObjectiveComparator());
            }
            selection = new TournamentSelection(2, comp);

            if (dutyCycleConstrained[1] || instrumentOrbitRelationsConstrained[1] || interferenceConstrained[1] || packingEfficiencyConstrained[1] || spacecraftMassConstrained[1] || synergyConstrained[1]) { // AOS objects

                // Setup for saving results
                properties.setBoolean("saveQuality", true);
                properties.setBoolean("saveCredits", true);
                properties.setBoolean("saveSelection", true);

                // IMPLEMENTATION WITH ACTUAL REPAIR OPERATORS

                ArrayList<Variation> operators = new ArrayList<>();
                for (int k = 0; k < heuristicsConstrained.length; k++) {
                    if (heuristicsConstrained[k][1]) {
                        operators.add(heuristicOperators[k]);
                    }
                }

                if (assigningProblem) {
                    operators.add(new CompoundVariation(new OnePointCrossover(crossoverProbability), new BitFlip(mutationProbability)));
                } else {
                    operators.add(new CompoundVariation(new PartitioningCrossover(crossoverProbability, params), new PartitioningMutation(mutationProbability, params)));
                }
                properties.setDouble("pmin", 0.03);

                // Create operator selector
                OperatorSelector operatorSelector = new AdaptivePursuit(operators, 0.8, 0.8, 0.1);
                //OperatorSelector operatorSelector = new ProbabilityMatching(operators, 0.6, 0.03);

                // Create credit assignment
                SetImprovementDominance creditAssignment = new SetImprovementDominance(archive, 1, 0);
                //OffspringParentDomination creditAssignment = new OffspringParentDomination(1.0, 0.5, 0.0, comp);

                // Create AOS
                aosStrategy = new AOSVariationSI(operatorSelector, creditAssignment, popSize);
                //aosStrategy = new AOSVariationOP(operatorSelector, creditAssignment, popSize);

                // Creating AOS MOEA object
                EpsilonMOEA emoea = new EpsilonMOEA(satelliteProblem, population, archive, selection, aosStrategy, initialization, comp);
                moeaObj = new AOSMOEA(emoea, aosStrategy, true);

            } else { // Epsilon MOEA objects
                if (assigningProblem) {
                    var = new CompoundVariation(new OnePointCrossover(crossoverProbability), new BitFlip(mutationProbability));
                } else {
                    var = new CompoundVariation(new PartitioningCrossover(crossoverProbability, params), new PartitioningMutation(mutationProbability, params));
                }
                moeaObj = new EpsilonMOEA(satelliteProblem, population, archive, selection, var, initialization, comp);
            }

            // Submit to ECS
            if (assigningProblem) {
                ecs.submit(new AssigningSearch(moeaObj, properties, savePath, "emoea_" + String.valueOf(i) + fileSaveNameConstraint + fileSaveNameProblem, evaluationManager));
            } else {
                ecs.submit(new PartitioningSearch(moeaObj, properties, savePath, "emoea_" + String.valueOf(i) + fileSaveNameConstraint + fileSaveNameProblem, params, evaluationManager));
            }
        }

        for (int i = 0; i < numRuns; i++) {
            try {
                Algorithm alg = ecs.take().get();
            } catch (InterruptedException | ExecutionException ex) {
                Logger.getLogger(MOEARun.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        evaluationManager.clear();
        pool.shutdown();
    }

}
