package seakers.vassartest;

import org.moeaframework.core.Solution;
import org.moeaframework.core.variable.BinaryVariable;
import seakers.orekit.util.OrekitConfig;
import seakers.vassar.ResultCollection;
import seakers.vassar.architecture.AbstractArchitecture;
import seakers.vassar.evaluation.AbstractArchitectureEvaluator;
import seakers.vassar.evaluation.ArchitectureEvaluationManager;
import seakers.vassar.io.ResultCollectionRecorder;
import seakers.vassar.problems.Assigning.*;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public class ReevaluateDataset {

    public static void main(String[] args) {
        SMAPJPL2Params params;

        String resourcesPath = "../VASSAR_resources";
        params = new SMAPJPL2Params(resourcesPath, "CRISP-ATTRIBUTES","test","normal");
        AbstractArchitectureEvaluator eval = new ArchitectureEvaluator();
        ArchitectureEvaluationManager AE = new ArchitectureEvaluationManager(params, eval);
        OrekitConfig.init(6, params.orekitResourcesPath);

        String csvFile = resourcesPath + "/problems/SMAP/results/2019-03-06_23-13-37_test.csv";
        String line = "";
        String cvsSplitBy = ",";

        ArrayList<AbstractArchitecture> archSet = new ArrayList<>();
        boolean header = true;
        try (BufferedReader br = new BufferedReader(new FileReader(csvFile))) {
            while ((line = br.readLine()) != null) {
                if (header) {
                    header = false;
                    continue;
                }
                // use comma as separator
                String[] csvArch = line.split(cvsSplitBy);
                Architecture arch = new Architecture(csvArch[0], 1, params);
                archSet.add(arch);
            }
        } catch (
                IOException e) {
            e.printStackTrace();
        }

        AE.init(6);
//        AE.setPopulation(archSet);
//        AE.evaluatePopulation();
        ResultCollection c = new ResultCollection(params, AE.getResults());

        ResultCollectionRecorder writer = new ResultCollectionRecorder(params);
        writer.write(c);

        OrekitConfig.end();
        AE.clear();
        System.out.println("DONE");

    }
}
