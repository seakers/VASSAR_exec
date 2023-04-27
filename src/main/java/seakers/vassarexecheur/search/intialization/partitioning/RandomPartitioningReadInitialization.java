package seakers.vassarexecheur.search.intialization.partitioning;

import org.moeaframework.core.Initialization;
import org.moeaframework.core.Solution;
import seakers.vassarexecheur.search.problems.partitioning.PartitioningArchitecture;
import seakers.vassarexecheur.search.problems.partitioning.PartitioningProblem;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

/**
 * This class creates the initial population for a Partitioning Problem run by reading
 * the appropriate csv file of pre-generated architectures
 *
 * @author roshansuresh
 */

public class RandomPartitioningReadInitialization implements Initialization {
    private final String readPath;
    private final int runNumber;
    private final String[] instrumentList;
    private final String[] orbitList;
    private final int populationSize;
    private final PartitioningProblem problem;

    public RandomPartitioningReadInitialization(String readPath, int runNumber, int populationSize, PartitioningProblem problem, String[] instrumentList, String[] orbitList) {
        this.readPath = readPath;
        this.runNumber = runNumber;
        this.instrumentList = instrumentList;
        this.orbitList = orbitList;
        this.populationSize = populationSize;
        this.problem = problem;
    }

    @Override
    public Solution[] initialize() {
        String csvFileName = readPath + File.separator + "Partitioning" + File.separator + "Injected Initialization" + File.separator + "random_feasible_partitioning" + Integer.toString(runNumber) + ".csv"; // Change accordingly

        Solution[] initialPopulation = new Solution[populationSize];

        // Read appropriate csv file
        List<List<String>> rows = new ArrayList<>();
        try (Scanner scanner = new Scanner(new File(csvFileName))) {
            while (scanner.hasNextLine()) {
                rows.add(getRecordFromLine(scanner.nextLine()));
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        int populationCounter = 0;
        boolean header = true;

        for (List<String> row : rows) {
            if (header) {
                header = false;
                continue;
            }
            PartitioningArchitecture newArch = (PartitioningArchitecture) problem.newSolution();
            int[] archInstrumentPartitioning = newArch.getInstrumentPartitionsFromString(row.get(0));
            int[] archOrbitAssigning = newArch.getOrbitAssignmentsFromString(row.get(0));
            newArch.setVariablesFromPartitionArrays(archInstrumentPartitioning, archOrbitAssigning);

            initialPopulation[populationCounter] = newArch;
            populationCounter++;
            if (populationCounter == populationSize) {
                break;
            }
        }

        return initialPopulation;
    }

    private List<String> getRecordFromLine(String line) {
        List<String> architectures = new ArrayList<>();
        try (Scanner rowScanner = new Scanner(line)) {
            rowScanner.useDelimiter(",");
            while (rowScanner.hasNext()) {
                architectures.add(rowScanner.next());
            }
        }
        return architectures;
    }


}
