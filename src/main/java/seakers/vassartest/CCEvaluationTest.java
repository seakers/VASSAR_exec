package seakers.vassartest;

import seakers.vassar.Result;
import seakers.vassar.architecture.AbstractArchitecture;
import seakers.vassar.evaluation.AbstractArchitectureEvaluator;
import seakers.vassar.evaluation.ArchitectureEvaluationManager;
import seakers.vassar.problems.Assigning.AssigningParams;
import seakers.vassar.problems.Assigning.ClimateCentricParams;
import seakers.vassar.problems.Assigning.Architecture;
import seakers.vassar.problems.Assigning.ArchitectureEvaluator;


public class CCEvaluationTest {

    public static void main(String[] args){

        //PATH
        String path = "./problems/ClimateCentric";

        AssigningParams params = new ClimateCentricParams(path, "CRISP-ATTRIBUTES",
                "test", "normal", "search_heuristic_rules_smap_127");
        AbstractArchitectureEvaluator eval = new ArchitectureEvaluator(params);
        ArchitectureEvaluationManager AE = new ArchitectureEvaluationManager(params, eval);
        AbstractArchitecture testArch;

        // List of instruments and orbits
//      {"ACE_ORCA", "ACE_POL", "ACE_LID", "CLAR_ERB", "ACE_CPR", "DESD_SAR", "DESD_LID", "GACM_VIS", "GACM_SWIR", "HYSP_TIR", "POSTEPS_IRS", "CNES_KaRIN"};
//      {"LEO-600-polar-NA", "SSO-600-SSO-AM", "SSO-600-SSO-DD", "SSO-800-SSO-DD", "SSO-800-SSO-PM"};

        //Reference rbsa.eoss.architecture #1
        //testArch = new Architecture("001001011110000001110101101000100101000101110010010011101101", 1, params);
        testArch = new Architecture("111100000000000000000000000000000000000000000000000000000000", 1, params);

        AE.init(1);
        Result result = AE.evaluateArchitectureSync(testArch, "Slow");
        AE.clear();

        System.out.println("DONE");
    }
}
