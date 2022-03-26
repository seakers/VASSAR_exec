package seakers.vassarexecheur.search.operators.assigning;

import jess.Fact;
import jess.JessException;
import jess.ValueVector;
import org.moeaframework.core.PRNG;
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
 * Checks to see if packing efficiency of the architecture falls within the acceptable bounds.
 * If not, one or more instruments are moved from one or more orbits to try to alleviate the situation
 *
 * @author roshansuresh
 */

public class RepairPackingEfficiencyAssigning implements Variation {
    /**
     * The packing efficiency that a spacecraft must be at or higher
     */
    private final double threshold;

    /**
     * The number of instruments to remove from each satellite that does not
     * meet the threshold
     */
    private final int xInstruments;

    //private final ParallelPRNG pprng;

    private final BaseParams params;

    public RepairPackingEfficiencyAssigning(double threshold, int xInstruments, BaseParams params) {
        this.xInstruments = xInstruments;
        this.threshold = threshold;
        this.params = params;
        //this.pprng = new ParallelPRNG();
    }

    @Override
    public int getArity() {
        return 1;
    }

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
            if ((operatorParameters.get(i).get(2)) < threshold) {
                candidateSatellites.add(i);
            }
        }

        if (candidateSatellites.size() > 0) {
            for (int k = 0; k < xInstruments; k++) {
                int satelliteIndex = PRNG.nextInt(candidateSatellites.size());
                ArrayList<String> currentPayloads = new ArrayList<>();
                if (childPayloads.isEmpty()) {
                    childOrbits.remove(satelliteIndex);
                    break;
                } else {
                    // Remove payload from current satellite
                    int payloadIndex = PRNG.nextInt(childPayloads.get(satelliteIndex).size());
                    String payload = childPayloads.get(satelliteIndex).get(payloadIndex);
                    ArrayList<Integer> candidatePayloadSatellites = getCandidateSatellitesForPayload(childPayloads, payload);
                    if (candidatePayloadSatellites.size() == 0) {
                        k--;
                        continue;
                    }
                    //for (int m = 0; m < childPayloads.get(satelliteIndex).size(); m++) {
                    //if (m != payloadIndex) {
                    //currentPayloads.add(childPayloads.get(satelliteIndex).get(m));
                    //}
                    //}
                    childPayloads.get(satelliteIndex).remove(payload);

                    // Add payload to different satellite
                    int candidatePayloadSatelliteIndex = PRNG.nextInt(candidatePayloadSatellites.size());
                    //ArrayList<String> candidateSatellitePayload = childPayloads.get(candidatePayloadSatelliteIndex);
                    //candidateSatellitePayload.add(payload);
                    childPayloads.get(candidatePayloadSatelliteIndex).add(payload);
                    childPayloads.set(satelliteIndex, currentPayloads);
                }
            }
        }
        AssigningArchitecture child = getArchitectureFromPayloadsAndOrbits(childPayloads, childOrbits);

        return new Solution[]{child};
    }

    private ArrayList<Integer> getCandidateSatellitesForPayload(ArrayList<ArrayList<String>> currentPayloads, String payload) {
        ArrayList<Integer> candidateSatelliteIndices = new ArrayList<>();
        for (int i = 0; i < currentPayloads.size(); i++) {
            ArrayList<String> currentSatellitePayloads = currentPayloads.get(i);
            if (!currentSatellitePayloads.contains(payload)) {
                candidateSatelliteIndices.add(i);
            }
        }
        return candidateSatelliteIndices;
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
