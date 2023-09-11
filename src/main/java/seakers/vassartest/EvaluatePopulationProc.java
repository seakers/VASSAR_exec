package seakers.vassartest;

import seakers.vassar.Result;
import seakers.vassar.evaluation.ArchitectureEvaluationManager;
import seakers.vassar.problems.Assigning.ArchitectureEvaluator;
import seakers.vassar.problems.Assigning.GigaAssigningParams;
import seakers.vassartest.search.problems.Assigning.AssigningProblem;
import seakers.vassartest.search.problems.Assigning.GigaArchitecture;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Future;

public class EvaluatePopulationProc {


    public static void main(String[] args) {
        String proc_num = args[0];
        int pop_lb = Integer.parseInt(args[1]);
        int pop_ub = Integer.parseInt(args[2]);

        // Problem Parameters
        int orekit_threads = 1;
        String resourcesPath = "/home/ec2-user/vassar/giga/VASSAR_resources" + proc_num;
        GigaAssigningParams params = new GigaAssigningParams(resourcesPath, "FUZZY-CASES", "test", "normal", orekit_threads);

        // Manager
        int numCpus = 10;
        ArchitectureEvaluator evaluator = new ArchitectureEvaluator();
        ArchitectureEvaluationManager evaluationManager = new ArchitectureEvaluationManager(params, evaluator);
        evaluationManager.init(numCpus);

        // Problem
        AssigningProblem problem = new AssigningProblem(new int[]{1}, "GigaProblem", evaluationManager, params);

        // --> PARALLEL
        List<String> population = GeneratePopulation.getPopulationChunkPruned(pop_lb, pop_ub);
        ArrayList<Future<Result>> futures = new ArrayList<Future<Result>>();
        for(String design: population){
            GigaArchitecture arch = new GigaArchitecture(design);
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
