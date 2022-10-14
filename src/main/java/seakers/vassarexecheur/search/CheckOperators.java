package seakers.vassarexecheur.search;

import org.moeaframework.core.Solution;
import org.moeaframework.core.Variable;
import org.moeaframework.core.Variation;
import org.moeaframework.core.variable.BinaryVariable;
import org.moeaframework.core.variable.RealVariable;
import org.moeaframework.problem.AbstractProblem;
import org.moeaframework.util.TypedProperties;
import seakers.architecture.util.IntegerVariable;
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
import java.util.HashMap;

/**
 * Used to check the working of each operator
 */

public class CheckOperators {

    public static void main(String[] args) {
        // Define problem parameters
        boolean assigningProblem = false; // True -> assigning problem, False -> partitioning problem

        double dcThreshold = 0.5;
        double massThreshold = 3000.0; // [kg]
        double packEffThreshold = 0.7;
        double instrCountThreshold = 15; // only for assigning problem
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
         * instrumentCountConstrained = [interior_penalty, AOS, biased_init, ACH, objective, constraint]
         *
         * if partitioning problem:
         * heuristicsConstrained = [dutyCycleConstrained, instrumentOrbitRelationsConstrained, interferenceConstrained, packingEfficiencyConstrained, spacecraftMassConstrained, synergyConstrained]
         * else:
         * heuristicsConstrained = [dutyCycleConstrained, instrumentOrbitRelationsConstrained, interferenceConstrained, packingEfficiencyConstrained, spacecraftMassConstrained, synergyConstrained, instrumentCountConstrained]
         */
        boolean[] dutyCycleConstrained = {true, false, false, false, false, false};
        boolean[] instrumentOrbitRelationsConstrained = {true, false, false, false, false, false};
        boolean[] interferenceConstrained = {true, false, false, false, false, false};
        boolean[] packingEfficiencyConstrained = {false, false, false, false, false, false};
        boolean[] spacecraftMassConstrained = {true, false, false, false, false, false};
        boolean[] synergyConstrained = {true, false, false, false, false, false};
        boolean[] instrumentCountConstrained = {false, false, false, false, false, false};

        boolean[][] heuristicsConstrained;
        if (assigningProblem) {
            heuristicsConstrained = new boolean[7][6];
        } else {
            heuristicsConstrained = new boolean[6][6];
        }

        for (int i = 0; i < heuristicsConstrained[0].length; i++) {
            heuristicsConstrained[0][i] = dutyCycleConstrained[i];
            heuristicsConstrained[1][i] = instrumentOrbitRelationsConstrained[i];
            heuristicsConstrained[2][i] = interferenceConstrained[i];
            heuristicsConstrained[3][i] = packingEfficiencyConstrained[i];
            heuristicsConstrained[4][i] =  spacecraftMassConstrained[i];
            heuristicsConstrained[5][i] = synergyConstrained[i];
            if (assigningProblem) {
                heuristicsConstrained[6][i] = instrumentCountConstrained[i];
            }
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
        HashMap<String, String[]> instrumentSynergyMap;
        HashMap<String, String[]> interferingInstrumentsMap;
        if (assigningProblem) {
            params = new ClimateCentricAssigningParams(resourcesPath, "FUZZY-ATTRIBUTES","test", "normal");

            instrumentSynergyMap = getInstrumentSynergyNameMap(params);
            interferingInstrumentsMap = getInstrumentInterferenceNameMap(params);

            evaluator = new ArchitectureEvaluator(considerFeasibility, interferingInstrumentsMap, instrumentSynergyMap, dcThreshold, massThreshold, packEffThreshold);

        } else {
            params = new ClimateCentricPartitioningParams(resourcesPath, "FUZZY-ATTRIBUTES", "test", "normal");

            instrumentSynergyMap = getInstrumentSynergyNameMap(params);
            interferingInstrumentsMap = getInstrumentInterferenceNameMap(params);

            evaluator = new seakers.vassarheur.problems.PartitioningAndAssigning.ArchitectureEvaluator(considerFeasibility, interferingInstrumentsMap, instrumentSynergyMap, dcThreshold, massThreshold, packEffThreshold);
        }
        ArchitectureEvaluationManager evaluationManager = new ArchitectureEvaluationManager(params, evaluator);
        evaluationManager.init(numCPU);

        // Problem class
        AbstractProblem satelliteProblem;
        if (assigningProblem) {
            satelliteProblem = new AssigningProblem(new int[]{1}, params.getProblemName(), evaluationManager, (ArchitectureEvaluator) evaluator, params, interferingInstrumentsMap, instrumentSynergyMap, dcThreshold, massThreshold, packEffThreshold, instrCountThreshold, numberOfHeuristicObjectives, numberOfHeuristicConstraints, heuristicsConstrained);
        } else {
            satelliteProblem = new PartitioningProblem(params.getProblemName(), evaluationManager, params, interferingInstrumentsMap, instrumentSynergyMap, dcThreshold, massThreshold, packEffThreshold, numberOfHeuristicObjectives, numberOfHeuristicConstraints, heuristicsConstrained);
        }

        // Initialize heuristic operator
        String operatorChoice = "instrCount"; // can be dutyCycle, instrOrbit, interInstr, packEff, spMass or instrSyn (or instrCount for assigning problem)
        Variation operator = null;
        String attribute = "";

        if (assigningProblem) {
            switch (operatorChoice) {
                case "dutyCycle":
                    operator = new RepairDutyCycleAssigning(dcThreshold, 1, params, false, (AssigningProblem) satelliteProblem, evaluationManager.getResourcePool(), (ArchitectureEvaluator) evaluator);
                    attribute = "DCViolation";
                    break;
                case "instrOrbit":
                    operator = new RepairInstrumentOrbitAssigning(1, evaluationManager.getResourcePool(), (ArchitectureEvaluator) evaluator, params, (AssigningProblem) satelliteProblem, true);
                    attribute = "InstrOrbViolation";
                    break;
                case "interInstr":
                    operator = new RepairInterferenceAssigning(1, evaluationManager.getResourcePool(), (ArchitectureEvaluator) evaluator, params, (AssigningProblem) satelliteProblem, interferingInstrumentsMap, false);
                    attribute = "InterInstrViolation";
                    break;
                case "packEff":
                    //operator = new RepairPackingEfficiencyAssigning(packEffThreshold, 1, params, moveInstrument, (AssigningProblem) satelliteProblem, evaluationManager.getResourcePool(), (ArchitectureEvaluator) evaluator);
                    operator = new RepairPackingEfficiencyAdditionAssigning(packEffThreshold, 1, 1, params, (AssigningProblem) satelliteProblem, evaluationManager.getResourcePool(), (ArchitectureEvaluator) evaluator);
                    attribute = "PackEffViolation";
                    break;
                case "spMass":
                    operator = new RepairMassAssigning(massThreshold, 1, params, false, (AssigningProblem) satelliteProblem, evaluationManager.getResourcePool(), (ArchitectureEvaluator) evaluator);
                    attribute = "SpMassViolation";
                    break;
                case "instrSyn":
                    //operator = new RepairSynergyAssigning(1, evaluationManager.getResourcePool(), (ArchitectureEvaluator) evaluator, params, (AssigningProblem) satelliteProblem, instrumentSynergyMap, moveInstrument);
                    operator = new RepairSynergyAdditionAssigning(1, evaluationManager.getResourcePool(), (ArchitectureEvaluator) evaluator, params, (AssigningProblem) satelliteProblem, instrumentSynergyMap);
                    attribute = "SynergyViolation";
                    break;
                case "instrCount":
                    operator = new RepairInstrumentCountAssigning(1, 1, instrCountThreshold, (AssigningProblem) satelliteProblem, evaluationManager.getResourcePool(), (ArchitectureEvaluator) evaluator, params);
                    attribute = "InstrCountViolation";
                    break;
                default: System.out.println("Invalid operator choice");
                    break;
            }
        } else {
            switch (operatorChoice) {
                case "dutyCycle":
                    operator = new RepairDutyCyclePartitioning(dcThreshold, 1, params, (PartitioningProblem) satelliteProblem, evaluationManager.getResourcePool(), (seakers.vassarheur.problems.PartitioningAndAssigning.ArchitectureEvaluator) evaluator);
                    attribute = "DCViolation";
                    break;
                case "instrOrbit":
                    operator = new RepairInstrumentOrbitPartitioning(1, evaluationManager.getResourcePool(), (seakers.vassarheur.problems.PartitioningAndAssigning.ArchitectureEvaluator) evaluator, params, (PartitioningProblem) satelliteProblem);
                    attribute = "InstrOrbViolation";
                    break;
                case "interInstr":
                    operator = new RepairInterferencePartitioning(1, evaluationManager.getResourcePool(), (seakers.vassarheur.problems.PartitioningAndAssigning.ArchitectureEvaluator) evaluator, params, (PartitioningProblem) satelliteProblem, interferingInstrumentsMap);
                    attribute = "InterInstrViolation";
                    break;
                case "packEff":
                    operator = new RepairPackingEfficiencyPartitioning(packEffThreshold, 1, params, (PartitioningProblem) satelliteProblem, evaluationManager.getResourcePool(), (seakers.vassarheur.problems.PartitioningAndAssigning.ArchitectureEvaluator) evaluator);
                    attribute = "PackEffViolation";
                    break;
                case "spMass":
                    operator = new RepairMassPartitioning(massThreshold, 1, params, (PartitioningProblem) satelliteProblem, evaluationManager.getResourcePool(), (seakers.vassarheur.problems.PartitioningAndAssigning.ArchitectureEvaluator) evaluator);
                    attribute = "SpMassViolation";
                    break;
                case "instrSyn":
                    operator = new RepairSynergyPartitioning(1, evaluationManager.getResourcePool(), (seakers.vassarheur.problems.PartitioningAndAssigning.ArchitectureEvaluator) evaluator, params, (PartitioningProblem) satelliteProblem, instrumentSynergyMap);
                    attribute = "SynergyViolation";
                    break;
                default: System.out.println("Invalid operator choice");
                    break;
            }
        }

        // Initialize architecture
        Solution child = null;
        if (assigningProblem) {
            String architectureString = "101110110101010101111110100101010000111001110110101101101101";
            AssigningArchitecture arch = new AssigningArchitecture(new int[]{1}, params.getNumInstr(), params.getNumOrbits(), 2);
            arch.setVariablesfromString(architectureString);

            // Evaluate and evolve architecture using operator for testing
            satelliteProblem.evaluate(arch);
            child = operator.evolve(new Solution[]{arch})[0];
            satelliteProblem.evaluate(child);

            evaluationManager.clear();
            System.out.println("Old Architecture Heuristic: " + arch.getAttribute(attribute));
            System.out.println("New Architecture Heuristic: " + child.getAttribute(attribute));
            System.out.println("DONE");

        } else {
            //int[] instrumentPartitioning = new int[]{0, 0, 0, 1, 1, 2, 3, 3, 3, 3, 4, 4};
            int[] instrumentPartitioning = new int[]{0, 0, 1, 2, 3, 4, 5, 6, 7, 7, 8, 9};
            int[] orbitAssignment = new int[]{3, 1, 0, 3, 2, 4, 3, 2, 0, 4, -1, -1};

            PartitioningArchitecture arch = new PartitioningArchitecture(params.getNumInstr(), params.getNumOrbits(), 2, params);
            arch.setVariablesFromPartitionArrays(instrumentPartitioning, orbitAssignment);

            // Evaluate and evolve architecture using operator for testing
            satelliteProblem.evaluate(arch);
            child = operator.evolve(new Solution[]{arch})[0];
            satelliteProblem.evaluate(child);

            evaluationManager.clear();
            System.out.println("Old Architecture Heuristic: " + arch.getAttribute(attribute));
            System.out.println("New Architecture Heuristic: " + child.getAttribute(attribute));
            System.out.println("DONE");
        }

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

}
