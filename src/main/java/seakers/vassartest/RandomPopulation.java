package seakers.vassartest;

import seakers.vassar.ResultCollection;
import seakers.vassar.io.ResultCollectionRecorder;
import seakers.vassar.architecture.AbstractArchitecture;
import seakers.vassar.evaluation.AbstractArchitectureEvaluator;
import seakers.vassar.evaluation.ArchitectureEvaluationManager;
import seakers.vassar.architecture.AbstractArchitectureGenerator;
import seakers.vassar.problems.Assigning.ArchitectureGenerator;
import seakers.vassar.problems.Assigning.ClimateCentricParams;
import seak.orekit.util.OrekitConfig;
import seakers.vassar.problems.Assigning.ArchitectureEvaluator;

import java.io.File;
import java.util.ArrayList;

public class RandomPopulation {
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        int POP_SIZE = 200;

        ClimateCentricParams params;

        String path = "." + File.separator + "problems" + File.separator + "ClimateCentric";
        params = new ClimateCentricParams(path,
                "CRISP-ATTRIBUTES","test","normal","search_heuristic_rules_smap_127");
        params.aggregationXls = params.path + "/xls/Aggregation Rules-Weather.xls";
        AbstractArchitectureEvaluator eval = new ArchitectureEvaluator(params);
        ArchitectureEvaluationManager AE = new ArchitectureEvaluationManager(params, eval);
        AbstractArchitectureGenerator archGenerator = new ArchitectureGenerator(params);
        OrekitConfig.init(6);

        ArrayList<AbstractArchitecture> initialPopulation = archGenerator.generateBiasedRandomPopulation(POP_SIZE, 0.25);
        AE.init(6);
        AE.setPopulation(initialPopulation);
        AE.evaluatePopulation();
        ResultCollection c = new ResultCollection(params, AE.getResults());

        ResultCollectionRecorder writer = new ResultCollectionRecorder(params);
        writer.write(c);

        OrekitConfig.end();
        AE.clear();
        System.out.println("DONE");
    }
}