package seakers.vassartest;

import seakers.orekit.util.OrekitConfig;
import seakers.vassar.Result;
import seakers.vassar.architecture.AbstractArchitecture;
import seakers.vassar.evaluation.AbstractArchitectureEvaluator;
import seakers.vassar.evaluation.ArchitectureEvaluationManager;
import seakers.vassar.problems.Assigning.Architecture;
import seakers.vassar.problems.Assigning.ArchitectureEvaluator;
import seakers.vassar.problems.Assigning.AssigningParams;
import seakers.vassar.problems.Assigning.ClimateCentricParams;

import py4j.GatewayServer;
import org.json.JSONObject;
import org.json.JSONArray;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.concurrent.Future;

public class IntegrationManager {


    public static void main(String[] args) {

        GatewayServer gatewayServer = new GatewayServer(new IntegrationManager());
        gatewayServer.start();
        System.out.println("Gateway Server Started");

    }

    public void getArchitectureScience(String archPath){
        try {
            // Getting params
            String resourcesPath = "C:\\Users\\demagall\\Documents\\IntelliJ\\VASSAR\\VASSAR_resources_dev\\VASSAR_resources";
            AssigningParams params = new ClimateCentricParams(resourcesPath, "CRISP-ATTRIBUTES",
                    "test", "normal");

            String content = new String(Files.readAllBytes(Paths.get(archPath)));
            JSONObject arch = new JSONObject(content);
            JSONArray constellationList = arch.getJSONArray("spaceSegment");
            ArchitectureEvaluator evaluator = new ArchitectureEvaluator();

            for (int i = 0; i < constellationList.length(); i++) {
                JSONObject constellation = constellationList.getJSONObject(i);
                Result result = evaluator.evaluatePerformanceFromJSON(constellation, params);
                System.out.println(result.toString());
            }

            System.out.println("DONE");

        } catch(Exception e) {
            e.printStackTrace();
        }

    }
}
