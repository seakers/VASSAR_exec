package seakers.vassarexecheur.search;

import org.moeaframework.core.Solution;
import org.moeaframework.core.Variable;
import org.moeaframework.core.Variation;
import org.moeaframework.core.variable.BinaryVariable;
import org.moeaframework.core.variable.IntegerVariable;
import org.moeaframework.problem.AbstractProblem;
import org.moeaframework.util.TypedProperties;
import seakers.vassarexecheur.search.operators.assigning.*;
import seakers.vassarexecheur.search.operators.partitioning.*;
import seakers.vassarexecheur.search.problems.assigning.AssigningArchitecture;
import seakers.vassarexecheur.search.problems.assigning.AssigningProblem;
import seakers.vassarexecheur.search.problems.partitioning.PartitioningArchitecture;
import seakers.vassarexecheur.search.problems.partitioning.PartitioningProblem;
import seakers.vassarheur.BaseParams;
import seakers.vassarheur.architecture.AbstractArchitecture;
import seakers.vassarheur.evaluation.AbstractArchitectureEvaluator;
import seakers.vassarheur.evaluation.ArchitectureEvaluationManager;
import seakers.vassarheur.problems.Assigning.Architecture;
import seakers.vassarheur.problems.Assigning.ArchitectureEvaluator;
import seakers.vassarheur.problems.Assigning.ClimateCentricAssigningParams;
import seakers.vassarheur.problems.PartitioningAndAssigning.ClimateCentricPartitioningParams;

import java.io.File;
import java.text.SimpleDateFormat;
import java.util.Date;

/**
 * Used to check the working of each operator
 */

public class CheckOperators {

    public static void main(String[] args) {
        // Define problem parameters
        boolean assigningProblem = true; // True -> assigning problem, False -> partitioning problem

        double dcThreshold = 0.5;
        double massThreshold = 3000.0; // [kg]
        double packEffThreshold = 0.4;
        boolean considerFeasibility = true;

        int numCPU = 1;

        String resourcesPath = "C:\\SEAK Lab\\SEAK Lab Github\\VASSAR\\VASSAR_resources-heur";

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

        BaseParams params;
        AbstractArchitectureEvaluator evaluator;
        if (assigningProblem) {
            params = new ClimateCentricAssigningParams(resourcesPath, "CRISP-ATTRIBUTES","test", "normal");
            evaluator = new ArchitectureEvaluator(considerFeasibility, dcThreshold, massThreshold, packEffThreshold);

        } else {
            params = new ClimateCentricPartitioningParams(resourcesPath, "CRISP-ATTRIBUTES", "test", "normal");
            evaluator = new seakers.vassarheur.problems.PartitioningAndAssigning.ArchitectureEvaluator(considerFeasibility, dcThreshold, massThreshold, packEffThreshold);
        }
        ArchitectureEvaluationManager evaluationManager = new ArchitectureEvaluationManager(params, evaluator);
        evaluationManager.init(numCPU);

        // Problem class
        AbstractProblem satelliteProblem;
        if (assigningProblem) {
            satelliteProblem = new AssigningProblem(new int[]{1}, params.getProblemName(), evaluationManager, params, dcThreshold, massThreshold, packEffThreshold, numberOfHeuristicObjectives, numberOfHeuristicConstraints, heuristicsConstrained);
        } else {
            satelliteProblem = new PartitioningProblem(params.getProblemName(), evaluationManager, params, dcThreshold, massThreshold, packEffThreshold, numberOfHeuristicObjectives, numberOfHeuristicConstraints, heuristicsConstrained);
        }

        // Initialize heuristic operator
        String operatorChoice = "dutyCycle"; // can be dutyCycle, instrOrbit, interInstr, packEff, spMass or instrSyn
        Variation operator = null;

        if (assigningProblem) {
            switch (operatorChoice) {
                case "dutyCycle": operator = new RepairDutyCycleAssigning(dcThreshold, 1, params);
                break;
                case "instrOrbit": operator = new RepairInstrumentOrbitAssigning(1, evaluationManager.getResourcePool(), (ArchitectureEvaluator) evaluator, params);
                break;
                case "interInstr": operator = new RepairInterferenceAssigning(1, evaluationManager.getResourcePool(), (ArchitectureEvaluator) evaluator, params);
                break;
                case "packEff": operator = new RepairPackingEfficiencyAssigning(packEffThreshold, 1, params);
                break;
                case "spMass": operator = new RepairMassAssigning(massThreshold, 1, params);
                break;
                case "instrSyn": operator = new RepairSynergyAssigning(1, evaluationManager.getResourcePool(), (ArchitectureEvaluator) evaluator, params);
                break;
                default: System.out.println("Invalid operator choice");
                break;
            }
        } else {
            switch (operatorChoice) {
                case "dutyCycle": operator = new RepairDutyCyclePartitioning(dcThreshold, 1, params);
                break;
                case "instrOrbit": operator = new RepairInstrumentOrbitPartitioning(1, evaluationManager.getResourcePool(), (seakers.vassarheur.problems.PartitioningAndAssigning.ArchitectureEvaluator) evaluator, params);
                break;
                case "interInstr": operator = new RepairInterferencePartitioning(1, evaluationManager.getResourcePool(), (seakers.vassarheur.problems.PartitioningAndAssigning.ArchitectureEvaluator) evaluator, params);
                break;
                case "packEff": operator = new RepairPackingEfficiencyPartitioning(packEffThreshold, 1, params);
                break;
                case "spMass": operator = new RepairMassPartitioning(massThreshold, 1, params);
                break;
                case "instrSyn": operator = new RepairSynergyPartitioning(1, evaluationManager.getResourcePool(), (seakers.vassarheur.problems.PartitioningAndAssigning.ArchitectureEvaluator) evaluator, params);
                break;
                default: System.out.println("Invalid operator choice");
                break;
            }
        }

        // Initialize architecture
        if (assigningProblem) {
            String architectureString = "110110101101101101110110010111011011111010110101101101011110";
            AssigningArchitecture arch = new AssigningArchitecture(new int[]{1}, params.getNumInstr(), params.getNumOrbits(), 2);
            // Populate arch with bits from architectureString
            for (int i = 0; i < architectureString.length(); i++) {
                BinaryVariable var = new BinaryVariable(1);
                String bit = architectureString.substring(i,i+1);
                if (bit.equalsIgnoreCase("1")) {
                    var.set(0, true);
                } else {
                    var.set(0, false);
                }
                arch.setVariable(i, var);
            }

            // Evaluate and evolve architecture using operator for testing
            satelliteProblem.evaluate(arch);
            operator.evolve(new Solution[]{arch});

        } else {
            int[] instrumentPartitioning = new int[]{2, 2, 2, 1, 4, 0, 3, 1, 4, 3, 0, 2};
            int[] orbitAssignment = new int[]{1, 1, 0, 2, 4, -1, -1, -1, -1, -1, -1, -1};

            PartitioningArchitecture arch = new PartitioningArchitecture(params.getNumInstr(), params.getNumOrbits(), 2, 1);
            // Populate arch with assigning and partitioning values
            for (int i = 0; i < 2*instrumentPartitioning.length; i++) {
                if (i < instrumentPartitioning.length) {
                    IntegerVariable var = new IntegerVariable(instrumentPartitioning[i], 0, instrumentPartitioning.length);
                    arch.setVariable(i, var);
                }
                if (i >= instrumentPartitioning.length) {
                    IntegerVariable var = new IntegerVariable(orbitAssignment[instrumentPartitioning.length - i], -1, params.getNumOrbits());
                    arch.setVariable(i, var);
                }
            }

            // Evaluate and evolve architecture using operator for testing
            satelliteProblem.evaluate(arch);
            operator.evolve(new Solution[]{arch});
        }




    }

}
