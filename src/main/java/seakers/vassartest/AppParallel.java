package seakers.vassartest;

import seakers.vassar.Result;
import seakers.vassar.evaluation.ArchitectureEvaluationManager;
import seakers.vassar.problems.Assigning.ArchitectureEvaluator;
import seakers.vassar.problems.Assigning.GigaAssigningParams;
import seakers.vassartest.search.problems.Assigning.AssigningProblem;
import seakers.vassartest.search.problems.Assigning.DesignBuilder;
import seakers.vassartest.search.problems.Assigning.GigaArchitecture;

import java.util.ArrayList;
import java.util.concurrent.Future;

public class AppParallel {

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
        int numCpus = 5;
        ArchitectureEvaluationManager evaluationManager = new ArchitectureEvaluationManager(params, evaluator);
        evaluationManager.init(numCpus);

        // Problem
        AssigningProblem problem = new AssigningProblem(new int[]{1}, "GigaProblem", evaluationManager, params);



        // Create Testing Design: 7175 Variables
        int num_insts = params.getNumInstr();
        int num_orbs = params.getNumOrbits();




        // --> LINEAR
//        String bit_string = DesignBuilder.randomBitString(num_insts, num_orbs, 40);
//        System.out.println(bit_string);
//        GigaArchitecture arch = new GigaArchitecture(bit_string);
//        Future<Result> result = problem.evaluateGigaArchAsync(arch);
//        try{
//            result.get();
//        }
//        catch (Exception ex){
//            ex.printStackTrace();
//        }

        // --> PARALLEL
        int num_designs = 5;
        ArrayList<Future<Result>> futures = new ArrayList<Future<Result>>();
        for(int x = 0; x < num_designs; x++){
            String bit_string = DesignBuilder.randomBitString(num_insts, num_orbs, 40);
            GigaArchitecture arch = new GigaArchitecture(bit_string);
            Future<Result> future = problem.evaluateGigaArchAsync(arch);
            futures.add(future);
        }
        int future_count = 0;
        for(Future<Result> future: futures){
            try{
                System.out.println("--> WAITING ON FUTURE: " + future_count);
                future.get();
            }
            catch (Exception ex){
                ex.printStackTrace();
            }
            future_count++;
        }

        System.exit(0);
    }


}
