package seakers.vassartest;

import seakers.vassar.BaseParams;
import seakers.vassar.Result;
import seakers.vassar.architecture.AbstractArchitecture;
import seakers.vassar.evaluation.AbstractArchitectureEvaluator;
import seakers.vassar.evaluation.ArchitectureEvaluationManager;
import seakers.vassar.problems.PartitioningAndAssigning.Architecture;
import seakers.vassar.problems.PartitioningAndAssigning.ArchitectureEvaluator;
import seakers.vassar.problems.PartitioningAndAssigning.Decadal2017AerosolsParams;

import java.util.*;


public class DecadalSurveyEvaluationTest {

    public static void main(String[] args){

        //PATH
        String path = "./problems/SMAP";

        BaseParams params = new Decadal2017AerosolsParams(path,"CRISP-ATTRIBUTES","test","normal","");
        AbstractArchitectureEvaluator eval = new ArchitectureEvaluator(params);
        ArchitectureEvaluationManager evalManager = new ArchitectureEvaluationManager(params, eval);
        AbstractArchitecture testArch;

        // List of instruments and orbits
//      {"BIOMASS","SMAP_RAD","SMAP_MWR","CMIS","VIIRS"};
//      {"LEO-600-polar-NA","SSO-600-SSO-AM","SSO-600-SSO-DD","SSO-800-SSO-AM","SSO-800-SSO-DD"};

        List<Set<String>> instrumentPartitioning = new ArrayList<>();
        Map<Set<String>, String> orbitAssignment = new HashMap<>();

        Set<String> set1 = new HashSet<>();
        set1.add("BIOMASS");
        instrumentPartitioning.add(set1);
        orbitAssignment.put(set1,"SSO-800-SSO-DD");

        Set<String> set2 = new HashSet<>();
        set2.add("SMAP_RAD");
        set2.add("SMAP_MWR");
        instrumentPartitioning.add(set2);
        orbitAssignment.put(set2,"LEO-600-polar-NA");

        Set<String> set3 = new HashSet<>();
        set3.add("CMIS");
        set3.add("VIIRS");
        instrumentPartitioning.add(set3);
        orbitAssignment.put(set3,"LEO-600-polar-NA");
        testArch = new Architecture(instrumentPartitioning, orbitAssignment, 1, params);

        System.out.println(testArch.ppString());

        evalManager.reset();
        evalManager.init(1);
        Result resu = evalManager.evaluateArchitectureSync(testArch, "Slow");

        System.out.println(resu.getScience());
        System.out.println(resu.getCost());
        System.out.println("DONE");
    }
}
