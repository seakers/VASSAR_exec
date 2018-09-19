package seakers.vassartest;

import org.moeaframework.core.Population;
import org.moeaframework.core.Solution;
import seakers.architecture.io.ResultIO;

import java.io.*;

public class PrintData {
    public static void main(String[] args) {
        int numRuns = 30;

        String path = "./problems/ClimateCentric/results/";
        String problemType = "weather";

        for (int run = 0; run < numRuns; ++run) {
            Population population;
            try {
                population = ResultIO.loadPopulation(path + "EpsilonMOEA_emoea_" + problemType + run + "_all.pop");
                String csvFile = path + problemType + run + ".csv";

                BufferedWriter outputWriter = new BufferedWriter(new FileWriter(csvFile));
                for (Solution solution: population) {
                    String bitString = "";
                    for (int i = 1; i < solution.getNumberOfVariables(); ++i) {
                        bitString += solution.getVariable(i).toString();
                    }
                    outputWriter.write(bitString);
                    outputWriter.write(",");
                    outputWriter.write(Double.toString(-solution.getObjectives()[0]));
                    outputWriter.write(",");
                    outputWriter.write(Double.toString(solution.getObjectives()[1]));
                    outputWriter.newLine();
                }

                outputWriter.close();
            }
            catch (IOException e) {
                e.printStackTrace();
            }
        }
    }
}
