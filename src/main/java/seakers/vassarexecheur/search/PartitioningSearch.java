package seakers.vassarexecheur.search;

import org.moeaframework.algorithm.AbstractEvolutionaryAlgorithm;
import org.moeaframework.core.Algorithm;
import org.moeaframework.core.Population;
import org.moeaframework.core.Solution;
import org.moeaframework.util.TypedProperties;
import seakers.aos.aos.AOS;
import seakers.aos.history.AOSHistoryIO;
import seakers.architecture.util.IntegerVariable;
import seakers.vassarheur.BaseParams;
import seakers.vassarexecheur.search.problems.partitioning.PartitioningArchitecture;
import seakers.vassarheur.evaluation.ArchitectureEvaluationManager;
import seakers.vassarheur.problems.PartitioningAndAssigning.Architecture;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashSet;
import java.util.Iterator;
import java.util.StringJoiner;
import java.util.concurrent.Callable;

public class PartitioningSearch implements Callable<Algorithm> {

    private final String savePath;
    private final String name;
    private final Algorithm alg;
    private final TypedProperties properties;
    private final ArchitectureEvaluationManager AEM;
    private final BaseParams params;

    public PartitioningSearch (Algorithm alg, TypedProperties properties, String savePath, String name, BaseParams params, ArchitectureEvaluationManager evalManager) {
        this.alg = alg;
        this.properties = properties;
        this.savePath = savePath;
        this.name = name;
        this.AEM = evalManager;
        this.params = params;
    }

    @Override
    public Algorithm call() throws Exception {
        int populationSize = (int) properties.getDouble("populationSize", 600);
        int maxEvaluations = (int) properties.getDouble("maxEvaluations", 10000);

        // run the executor using the listener to collect results
        System.out.println("Starting " + alg.getClass().getSimpleName() + " on " + alg.getProblem().getName() + " with pop size: " + populationSize);
        alg.step();
        long startTime = System.currentTimeMillis();

        HashSet<Solution> allSolutions = new HashSet<>();
        Population initPop = ((AbstractEvolutionaryAlgorithm) alg).getPopulation();
        for (int i = 0; i < initPop.size(); i++) {
            initPop.get(i).setAttribute("NFE", 0);
            allSolutions.add( initPop.get(i));
        }

        while (!alg.isTerminated() && (alg.getNumberOfEvaluations() < maxEvaluations)) {
            if (alg.getNumberOfEvaluations() % 100 == 0) {
                System.out.println("NFE: " + alg.getNumberOfEvaluations());
                System.out.print("Popsize: " + ((AbstractEvolutionaryAlgorithm) alg).getPopulation().size());
                System.out.println("  Archivesize: " + ((AbstractEvolutionaryAlgorithm) alg).getArchive().size());
            }
            alg.step();
            Population pop = ((AbstractEvolutionaryAlgorithm) alg).getPopulation();
            for(int i=1; i<3; i++){
                Solution s = pop.get(pop.size() - i);
                s.setAttribute("NFE", alg.getNumberOfEvaluations());
                allSolutions.add(s);
            }
            if (alg.getNumberOfEvaluations() % 100 == 0) {
                AEM.getResourcePool().poolClean();
                System.out.println("Rete clean initiated");
            }
        }

        Population pop = ((AbstractEvolutionaryAlgorithm) alg).getPopulation();
        String filename = savePath + File.separator + alg.getClass().getSimpleName() + "_" + name + ".csv";
        savePopulationCSV(pop, filename);

        String fullPopulationFilename = savePath + File.separator + alg.getClass().getSimpleName() + "_" + name + "_fullpop.csv";
        savePopulationHistory(allSolutions, fullPopulationFilename);

        alg.terminate();
        long finishTime = System.currentTimeMillis();
        System.out.println("Done with optimization. Execution time: " + ((finishTime - startTime) / 1000) + "s");

        if (alg instanceof AOS) {
            AOS algAOS = (AOS) alg;
            if (properties.getBoolean("saveQuality", false)) {
                AOSHistoryIO.saveQualityHistory(algAOS.getQualityHistory(), new File(savePath + File.separator + name + "_qual" + ".csv"), ",");
            }
            if (properties.getBoolean("saveCredits", false)) {
                AOSHistoryIO.saveCreditHistory(algAOS.getCreditHistory(), new File(savePath + File.separator + name + "_credit" + ".csv"), ",");
            }
            if (properties.getBoolean("saveSelection", false)) {
                AOSHistoryIO.saveSelectionHistory(algAOS.getSelectionHistory(), new File(savePath + File.separator + name + "_hist" + ".csv"), ",");
            }
        }
        return alg;
    }

    private void savePopulationCSV(Population pop, String filename) {

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
                sj.add(arch_abs.toString(""));
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

    private void savePopulationHistory(HashSet<Solution> solutionSet, String filename) {
        File results = new File(filename);
        results.getParentFile().mkdirs();

        try (FileWriter writer = new FileWriter(results)) {

            StringJoiner headings = new StringJoiner(",");
            headings.add("Architecture");
            headings.add("NFE");
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

            Iterator<Solution> iter = solutionSet.iterator();
            while(iter.hasNext()){

                Solution sol = iter.next();
                int numberOfFunctionEvaluations = (int) sol.getAttribute("NFE");

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

                Architecture arch_abs= new Architecture(instrumentPartitioning, orbitAssignment, 1, params);

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
                sj.add(arch_abs.toString(""));
                sj.add(Integer.toString(numberOfFunctionEvaluations));
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
}
