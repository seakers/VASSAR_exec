package seakers.vassartest;

import com.google.gson.JsonObject;
import seakers.vassar.evaluation.ArchitectureEvaluationManager;
import seakers.vassar.problems.Assigning.ArchitectureEvaluator;
import seakers.vassar.problems.Assigning.GigaAssigningParams;
import seakers.vassartest.search.problems.Assigning.AssigningProblem;
import seakers.vassartest.search.problems.Assigning.DesignBuilder;
import seakers.vassartest.search.problems.Assigning.GigaArchitecture;

import com.google.gson.JsonParser;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

public class DebugRun {

    public static void main (String[] args) {
        String dir_path = "/home/ec2-user/designs2_100k/designs_all";




        // Problem Parameters
        int orekit_threads = 1;
        String resourcesPath = "/home/ec2-user/vassar/giga/VASSAR_resources";
        GigaAssigningParams params = new GigaAssigningParams(resourcesPath, "FUZZY-CASES", "test", "normal", orekit_threads);
        ArchitectureEvaluator evaluator = new ArchitectureEvaluator();

        // Manager
        int numCpus = 1;
        ArchitectureEvaluationManager evaluationManager = new ArchitectureEvaluationManager(params, evaluator);
        evaluationManager.init(numCpus);

        // Problem
        AssigningProblem problem = new AssigningProblem(new int[]{1}, "GigaProblem", evaluationManager, params);






        JsonObject error_design = getErrorDesign(dir_path);
        System.out.println(error_design.get("science").getAsInt());
        String bit_string = error_design.get("design").getAsString();
        GigaArchitecture arch = new GigaArchitecture(bit_string);

        // Evaluate
        problem.evaluate(arch);






        System.exit(0);
    }

    public static JsonObject getErrorDesign(String path_str){
        File directory = new File(path_str);
        File[] files = directory.listFiles((dir, name) -> name.endsWith(".json"));
        if (files != null) {
            JsonParser parser = new JsonParser();
            for (File file : files) {
                try (FileReader reader = new FileReader(file)) {
                    JsonObject jsonObject = parser.parse(reader).getAsJsonObject();
                    if (jsonObject.has("science") && jsonObject.get("science").getAsInt() == -1) {
                        return jsonObject;
                    }
                } catch (IOException e) {
                    System.err.println("Error reading file: " + file.getName());
                }
            }
        }
        return null;
    }



}
