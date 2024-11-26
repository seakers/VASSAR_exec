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

//        String archPath = "C:\\Users\\demagall\\Documents\\VS Code\\Research\\SpaDes\\adjArch.json";
//        Double science = getArchitectureScience(archPath, 10.);
//        System.out.println("Science: " + science);
    }

    public static Double getArchitectureScience(String archPath, Double revisit){
        Double science = 0.0;
        try {
            // Getting params
            String resourcesPath = "C:\\Users\\demagall\\Documents\\IntelliJ\\VASSAR\\VASSAR_resources_dev\\VASSAR_resources";
            AssigningParams params = new ClimateCentricParams(resourcesPath, "CRISP-ATTRIBUTES",
                    "test", "normal");

            String content = new String(Files.readAllBytes(Paths.get(archPath)));
            JSONObject arch = new JSONObject(content);
//            JSONArray constList = arch.getJSONArray("spaceSegment");
//            for (int i = 0; i < constList.length(); i++) {
//                JSONObject constellation = constList.getJSONObject(i);
//                JSONArray satList = constellation.getJSONArray("satellites");
//                for (int j = 0; j < satList.length(); j++) {
//                    JSONObject satellite = satList.getJSONObject(j);
//                    JSONArray instrList = satellite.getJSONArray("instruments");
//                    for (int k = 0; k < instrList.length(); k++) {
//                        JSONObject instrument = instrList.getJSONObject(k);
//                        String instrName = instrument.getString("name");
//                        String instrOrbit = satellite.getString("orbit");
//                        System.out.println(instrName + " " + instrOrbit);
//                    }
//                }
//            }
//            JSONObject mission = arch.getJSONObject("spaceSegment");
            ArchitectureEvaluator evaluator = new ArchitectureEvaluator();

//            for (int i = 0; i < constellationList.length(); i++) {
//                JSONObject constellation = constellationList.getJSONObject(i);
//                Result result = evaluator.evaluatePerformanceFromJSON(constellation, params);
//                System.out.println(result.toString());
//            }
//            JSONObject constellation = constellationList.getJSONObject;
            Result result = evaluator.evaluatePerformanceFromJSON(arch, revisit, params);
//            System.out.println(result.toString());
            science = result.getScience();
            System.out.print("Science: " + science);

            System.out.println("DONE");

        } catch(Exception e) {
            e.printStackTrace();
        }

        return science;

    }
}
