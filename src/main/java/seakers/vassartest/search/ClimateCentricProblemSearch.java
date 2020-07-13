/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package seakers.vassartest.search;


import org.moeaframework.algorithm.AbstractEvolutionaryAlgorithm;
import org.moeaframework.core.Algorithm;
import org.moeaframework.core.Population;
import org.moeaframework.core.Solution;
import org.moeaframework.util.TypedProperties;
import seakers.vassartest.search.problems.Assigning.AssigningArchitecture;
import seakers.vassartest.search.problems.Assigning.AssigningProblem;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.concurrent.Callable;

/**
 *
 * @author nozomihitomi
 */
public class ClimateCentricProblemSearch implements Callable<Algorithm> {

    private final String savePath;
    private final String name;
    private Algorithm alg;
    private final TypedProperties properties;

    public ClimateCentricProblemSearch(Algorithm alg, TypedProperties properties, String savePath, String name) {
        this.alg = alg;
        this.properties = properties;
        this.savePath = savePath;
        this.name = name;
    }

    @Override
    public Algorithm call() throws IOException  {

        int populationSize = (int) properties.getDouble("populationSize", 600);
        int maxEvaluations = (int) properties.getDouble("maxEvaluations", 10000);

        // run the executor using the listener to collect results
        System.out.println("Starting " + alg.getClass().getSimpleName() + " on " + alg.getProblem().getName() + " with pop size: " + populationSize);
        alg.step();
        long startTime = System.currentTimeMillis();

        Population initPop = ((AbstractEvolutionaryAlgorithm) alg).getPopulation();
        for (int i = 0; i < initPop.size(); i++) {
            initPop.get(i).setAttribute("NFE", 0);
        }

        int cnt = 0;
        while (!alg.isTerminated() && (alg.getNumberOfEvaluations() < maxEvaluations)) {
            if (alg.getNumberOfEvaluations() % 500 == 0) {
                System.out.println("NFE: " + alg.getNumberOfEvaluations());
                System.out.print("Popsize: " + ((AbstractEvolutionaryAlgorithm) alg).getPopulation().size());
                System.out.println("  Archivesize: " + ((AbstractEvolutionaryAlgorithm) alg).getArchive().size());

                Population pop = ((AbstractEvolutionaryAlgorithm) alg).getPopulation();
                String filename = savePath + File.separator + alg.getClass().getSimpleName() + "_" + name + "_" + cnt;
                savePopulationCSV(pop, filename);
                cnt++;

                ((AssigningProblem)alg.getProblem()).resetEvaluationManager();
            }
            alg.step();
        }

        Population pop = ((AbstractEvolutionaryAlgorithm) alg).getPopulation();
        String filename = savePath + File.separator + alg.getClass().getSimpleName() + "_" + name + "_" + cnt;
        savePopulationCSV(pop, filename);

        alg.terminate();
        long finishTime = System.currentTimeMillis();
        System.out.println("Done with optimization. Execution time: " + ((finishTime - startTime) / 1000) + "s");
        return alg;
    }


    public void savePopulationCSV(Population pop, String filename) {

        File results = new File(filename);
        results.getParentFile().mkdirs();

        System.out.println("Saving a population in a csv file");

        try (FileWriter writer = new FileWriter(results)) {

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

                StringJoiner sj = new StringJoiner(",");
                sj.add(bitString);
                sj.add(Double.toString(science));
                sj.add(Double.toString(cost));
                writer.append(sj.toString());
                writer.append("\n");
            }
            writer.flush();

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}
