package seakers.vassarexecheur.search.operators.assigning;

import jess.Fact;
import jess.JessException;
import jess.ValueVector;
import org.moeaframework.core.ParallelPRNG;
import org.moeaframework.core.Solution;
import org.moeaframework.core.Variation;
import jess.Rete;
import org.moeaframework.core.variable.BinaryVariable;
import seakers.vassarheur.BaseParams;
import seakers.vassarheur.QueryBuilder;
import seakers.vassarexecheur.search.problems.assigning.AssigningArchitecture;
import seakers.vassarheur.problems.Assigning.ClimateCentricAssigningParams;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;


/**
 * Checks that the data rate duty cycle of the spacecraft fall within the
 * acceptable bounds. If not one or more random instruments are removed to try
 * to alleviate the situation. User can define whether to change one or multiple
 * spacecraft. Class modified from the original RepairDutyCycle.java in the EOSS
 * repo authored by nozomihitomi
 *
 * @author roshansuresh
 */

public class RepairDutyCycle implements Variation {

    /**
     * The duty cycle that a spacecraft must be at or higher
     */
    private final double threshold;

    /**
     * The number of instruments to remove from each satellite that does not
     * meet the threshold
     */
    private final int xInstruments;

    /**
     * The number of satellites to modify
     */
    private final int ySatellites;

    private final ParallelPRNG pprng;

    private final BaseParams params;

    public RepairDutyCycle(double threshold, int xInstruments, int ySatellites, BaseParams params) {
        this.xInstruments = xInstruments;
        this.ySatellites = ySatellites;
        this.threshold = threshold;
        this.params = params;
        this.pprng = new ParallelPRNG();
    }

    /**
     * removes x number of instruments from the payload of y number of satellite
     * that does not meet the data rate duty cycle threshold
     *
     * @param solutions
     * @return
     */
    @Override
    public Solution[] evolve(Solution[] solutions) {
        AssigningArchitecture parent = (AssigningArchitecture) solutions[0];
        ArrayList<ArrayList<String>> payloads = parent.getSatellitePayloads();
        ArrayList<String> orbits = parent.getSatelliteOrbits();

        ArrayList<ArrayList<String>> childPayloads = new ArrayList<>(payloads);
        ArrayList<String> childOrbits = new ArrayList<>(orbits);

        ArrayList<ArrayList<Double>> operatorParameters = parent.getOperatorParameters(); //{duty cycle, wet mass, packing efficiency} for each satellite
        ArrayList<Integer> candidateSatellites = new ArrayList<>();

        for (int i = 0; i < payloads.size(); i++) {
            if ((operatorParameters.get(i).get(0)) < threshold) {
                candidateSatellites.add(i);
            }
        }

        for (int j = 0; j < ySatellites; j++) {
            if (j > childPayloads.size() || j >= candidateSatellites.size()) {
                break;
            }
            int satelliteIndex = pprng.nextInt(candidateSatellites.size());
            for (int k = 0; k < xInstruments; k++) {
                ArrayList<String> currentPayloads = new ArrayList<>();
                if (childPayloads.isEmpty()) {
                    childOrbits.remove(satelliteIndex);
                    break;
                } else {
                    int payloadIndex = pprng.nextInt(childPayloads.get(satelliteIndex).size());
                    for (int m = 0; m < childPayloads.get(satelliteIndex).size(); m++) {
                        if (m != payloadIndex) {
                            currentPayloads.add(childPayloads.get(satelliteIndex).get(m));
                        }
                    }
                    childPayloads.set(satelliteIndex, currentPayloads);
                }
            }
        }
        AssigningArchitecture child = getArchitectureFromPayloadsAndOrbits(childPayloads, childOrbits);

        return new Solution[]{child};
    }

    @Override
    public int getArity() {
        return 1;
    }

    private AssigningArchitecture getArchitectureFromPayloadsAndOrbits(ArrayList<ArrayList<String>> currentPayloads, ArrayList<String> currentOrbits) {
        // ORDER ORBITS AND INSTRUMENTS BASED ON CLIMATECENTRICPARAMS
        ArrayList<String> instrumentList = new ArrayList<>(Arrays.asList(params.getInstrumentList()));
        ArrayList<String> orbitList = new ArrayList<>(Arrays.asList(params.getOrbitList()));

        AssigningArchitecture arch = new AssigningArchitecture(new int[]{1}, params.getNumInstr(), params.getNumOrbits(), 2);

        for (int i = 0; i < params.getNumOrbits(); i++) {
            if (!currentOrbits.contains(orbitList.get(i))) {
                for (int j = 0; j < params.getNumInstr(); j++) {
                    BinaryVariable var = new BinaryVariable(1);
                    var.set(0, false);
                    int decisionIndex = params.getNumInstr()*i + j;
                    arch.setVariable(decisionIndex, var);
                }
            } else {
                ArrayList<String> currentInstruments = currentPayloads.get(i);
                for (int j = 0; j < params.getNumInstr(); j++) {
                    BinaryVariable var = new BinaryVariable(1);
                    int decisionIndex = params.getNumInstr()*i + j;
                    if (!currentInstruments.contains(instrumentList.get(j))) {
                        var.set(0, false);
                    } else {
                        var.set(0, true);
                    }
                    arch.setVariable(decisionIndex, var);
                }
            }
        }
        return arch;
    }


}
