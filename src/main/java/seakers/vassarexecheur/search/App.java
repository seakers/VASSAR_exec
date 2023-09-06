package seakers.vassarexecheur.search;



import org.moeaframework.core.*;
import org.moeaframework.core.variable.BinaryVariable;
import seakers.vassarexecheur.search.problems.assigning.AssigningArchitecture;
import seakers.vassarexecheur.search.problems.assigning.AssigningProblem;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.*;
import java.util.ArrayList;
import seakers.vassarheur.BaseParams;
import seakers.vassarheur.evaluation.ArchitectureEvaluationManager;
import seakers.vassarheur.problems.Assigning.ArchitectureEvaluator;
import seakers.vassarheur.problems.Assigning.ClimateCentricAssigningParams;
import seakers.vassarheur.problems.Assigning.GigaAssigningParams;

import seakers.vassarheur.Result;





public class App {


    public static void main (String[] args) {
        System.out.println("--> TESTING");



//        String resourcesPath = "/home/ec2-user/vassar/VASSAR_resources";
//        ClimateCentricAssigningParams l_params = new ClimateCentricAssigningParams(resourcesPath, "FUZZY-ATTRIBUTES", "test", "normal");

        // Problem Parameters
        String resourcesPath = "/home/ec2-user/vassar/giga/VASSAR_resources";
        GigaAssigningParams params = new GigaAssigningParams(resourcesPath, "FUZZY-ATTRIBUTES", "test", "normal");




//        ArchitectureEvaluator evaluator = new ArchitectureEvaluator(considerFeasibility, interferingInstrumentsMap, instrumentSynergyMap, dcThreshold, massThreshold, packEffThreshold);



    }


}
