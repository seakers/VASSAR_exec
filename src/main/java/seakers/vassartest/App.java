package seakers.vassartest;



import org.moeaframework.core.*;
import org.moeaframework.core.variable.BinaryVariable;
import seakers.vassartest.search.problems.Assigning.AssigningArchitecture;
import seakers.vassartest.search.problems.Assigning.AssigningProblem;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.*;
import java.util.ArrayList;
import seakers.vassar.BaseParams;
import seakers.vassar.evaluation.ArchitectureEvaluationManager;
import seakers.vassar.problems.Assigning.ArchitectureEvaluator;
//import seakers.vassar.problems.Assigning.ClimateCentricAssigningParams;
import seakers.vassar.problems.Assigning.GigaAssigningParams;

import seakers.vassar.Result;
import seakers.vassartest.search.problems.Assigning.DesignBuilder;
import seakers.vassartest.search.problems.Assigning.GigaArchitecture;


public class App {


    public static void main (String[] args) {
        System.out.println("--> TESTING");
        int orekit_threads = 1;



//        String resourcesPath = "/home/ec2-user/vassar/VASSAR_resources";
//        ClimateCentricAssigningParams l_params = new ClimateCentricAssigningParams(resourcesPath, "FUZZY-ATTRIBUTES", "test", "normal");

        // Problem Parameters
        String resourcesPath = "/home/ec2-user/vassar/giga/VASSAR_resources";
        GigaAssigningParams params = new GigaAssigningParams(resourcesPath, "FUZZY-CASES", "test", "normal", orekit_threads);

        // Evaluator
        ArchitectureEvaluator evaluator = new ArchitectureEvaluator();

        // Manager
        int numCpus = 1;
        ArchitectureEvaluationManager evaluationManager = new ArchitectureEvaluationManager(params, evaluator);
        evaluationManager.init(numCpus);

        // Problem
        AssigningProblem problem = new AssigningProblem(new int[]{1}, "GigaProblem", evaluationManager, params);





        // Create Testing Design: 7175 Variables
        int num_insts = params.getNumInstr();
        int num_orbs = params.getNumOrbits();


//        String bit_string = DesignBuilder.randomBitString(num_insts, num_orbs, 5);


//        ArrayList<Integer> indices = new ArrayList<>();
//        int inst_count = 35;
//        int offset = 41 * 100 + 0;
//        inst_count += offset;
//        for(int x = offset; x < inst_count; x++){
//            indices.add(x);
//        }
//        String bit_string = DesignBuilder.specificIdxBitString(num_insts, num_orbs, indices);
////        String bit_string = DesignBuilder.specificOrbInstBitString(num_insts, num_orbs, 41, 26);
//
//        // Architecture
//        GigaArchitecture arch = new GigaArchitecture(bit_string);
//
//        // Evaluate
//        problem.evaluate(arch);



        // Eval in loop
        int num_designs = 10;
        for(int x = 0; x < num_designs; x++){

            String bit_string = DesignBuilder.randomBitStringOrbs(num_insts, num_orbs, 10);

            // Architecture
            GigaArchitecture arch = new GigaArchitecture(bit_string);

            // Evaluate
            problem.evaluate(arch);
        }






        System.exit(0);
    }


}