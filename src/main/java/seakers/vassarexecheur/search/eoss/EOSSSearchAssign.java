package seakers.vassarexecheur.search.eoss;

import eoss.problem.assignment.InstrumentAssignmentArchitecture;
import org.moeaframework.algorithm.AbstractEvolutionaryAlgorithm;
import org.moeaframework.core.Algorithm;
import org.moeaframework.core.Population;
import org.moeaframework.core.Solution;
import org.moeaframework.util.TypedProperties;
//import seakers.aos.aos.AOS;
//import seakers.aos.history.AOSHistoryIO;
//import seakers.architecture.io.ResultIO;
//import seakers.vassarheur.search.problems.Assigning.AssigningArchitecture;
import aos.IO.IOCreditHistory;
import aos.IO.IOQualityHistory;
import aos.IO.IOSelectionHistory;
import aos.aos.IAOS;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashSet;
import java.util.Iterator;
import java.util.StringJoiner;
import java.util.concurrent.Callable;

public class EOSSSearchAssign implements Callable<Algorithm> {

    private final String savePath;
    private final String name;
    private final Algorithm alg;
    private final TypedProperties properties;

    public EOSSSearchAssign(Algorithm alg, TypedProperties properties, String savePath, String name) {
        this.alg = alg;
        this.properties = properties;
        this.savePath = savePath;
        this.name = name;
    }

    @Override
    public Algorithm call() throws IOException {

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
            if (alg.getNumberOfEvaluations() % 500 == 0) {
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
        }

        Population pop = ((AbstractEvolutionaryAlgorithm) alg).getPopulation();
        String filename = savePath + File.separator + alg.getClass().getSimpleName() + "_" + name + ".csv";
        savePopulationCSV(pop, filename);

        String fullPopulationFilename = savePath + File.separator + alg.getClass().getSimpleName() + "_" + name + "_fullpop.csv";
        savePopulationHistory(allSolutions, fullPopulationFilename);

        alg.terminate();
        long finishTime = System.currentTimeMillis();
        System.out.println("Done with optimization. Execution time: " + ((finishTime - startTime) / 1000) + "s");

        //String filename = savePath + File.separator + alg.getClass().getSimpleName() + "_" + name;
        //ResultIO.savePopulation(((AbstractEvolutionaryAlgorithm) alg).getPopulation(), filename);
        //ResultIO.savePopulation(new Population(allSolutions), filename + "_all");
        //ResultIO.saveObjectives(alg.getResult(), filename);

        //if (alg instanceof AOS) {
            //AOS algAOS = (AOS) alg;
            //if (properties.getBoolean("saveQuality", false)) {
                //AOSHistoryIO.saveQualityHistory(algAOS.getQualityHistory(), new File(savePath + File.separator + name + ".qual"), ",");
            //}
            //if (properties.getBoolean("saveCredits", false)) {
                //AOSHistoryIO.saveCreditHistory(algAOS.getCreditHistory(), new File(savePath + File.separator + name + ".credit"), ",");
            //}
            //if (properties.getBoolean("saveSelection", false)) {
                //AOSHistoryIO.saveSelectionHistory(algAOS.getSelectionHistory(), new File(savePath + File.separator + name + ".hist"), ",");
            //}
        //}

        if (alg instanceof IAOS) {
            IAOS algAOS = (IAOS) alg;
            if (properties.getBoolean("saveQuality", false)) {
                IOQualityHistory ioqh = new IOQualityHistory();
                ioqh.saveHistory(algAOS.getQualityHistory(), savePath + File.separator + name + ".qual", ",");
            }
            if (properties.getBoolean("saveCredits", false)) {
                IOCreditHistory ioch = new IOCreditHistory();
                ioch.saveHistory(algAOS.getCreditHistory(), savePath + File.separator + name + ".credit", ",");
            }
            if (properties.getBoolean("saveSelection", false)) {
                IOSelectionHistory iosh = new IOSelectionHistory();
                iosh.saveHistory(algAOS.getSelectionHistory(), savePath + File.separator + name + ".hist", ",");
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
