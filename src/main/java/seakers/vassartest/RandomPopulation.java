package seakers.vassartest;

import seakers.vassar.ResultCollection;
import seakers.vassar.io.ResultCollectionRecorder;
import seakers.vassar.architecture.AbstractArchitecture;
import seakers.vassar.evaluation.AbstractArchitectureEvaluator;
import seakers.vassar.evaluation.ArchitectureEvaluationManager;
import seakers.vassar.architecture.AbstractArchitectureGenerator;
import seakers.vassar.problems.Assigning.*;
import seakers.orekit.util.OrekitConfig;

import java.io.File;
import java.util.ArrayList;

public class RandomPopulation {
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        int POP_SIZE = 500;

        SMAPParams params;

        String resourcesPath = "../VASSAR_resources";
        params = new SMAPParams(resourcesPath,
                "CRISP-ATTRIBUTES","test","normal");
        AbstractArchitectureEvaluator eval = new ArchitectureEvaluator();
        ArchitectureEvaluationManager AE = new ArchitectureEvaluationManager(params, eval);
        AbstractArchitectureGenerator archGenerator = new ArchitectureGenerator(params);
        OrekitConfig.init(6, params.orekitResourcesPath);


        ArrayList<AbstractArchitecture> initialPopulation = archGenerator.generateBiasedRandomPopulation(POP_SIZE, 0.15);
        AE.init(6);
//        AE.setPopulation(initialPopulation);
//        AE.evaluatePopulation();
        ResultCollection c = new ResultCollection(params, AE.getResults());

        ResultCollectionRecorder writer = new ResultCollectionRecorder(params);
        writer.write(c);

        OrekitConfig.end();
        AE.clear();
        System.out.println("DONE");
    }
}
