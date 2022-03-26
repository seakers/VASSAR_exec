package seakers.vassarexecheur.search.operators.partitioning;

import jess.Fact;
import jess.JessException;
import jess.ValueVector;
import org.moeaframework.core.PRNG;
import org.moeaframework.core.Solution;
import org.moeaframework.core.Variation;
import jess.Rete;
import org.moeaframework.core.variable.BinaryVariable;
import seakers.architecture.util.IntegerVariable;
import seakers.vassarheur.BaseParams;
import seakers.vassarheur.QueryBuilder;
import seakers.vassarexecheur.search.problems.partitioning.PartitioningArchitecture;
import seakers.vassarheur.problems.PartitioningAndAssigning.ClimateCentricPartitioningParams;

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

public class RepairPackingEfficiencyPartitioning implements Variation {
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

    public RepairPackingEfficiencyPartitioning(double threshold, int xInstruments, BaseParams params) {
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
        PartitioningArchitecture parent = (PartitioningArchitecture) solutions[0];
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
        PartitioningArchitecture child = getArchitectureFromPayloadsAndOrbits(childPayloads, childOrbits);

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

    private PartitioningArchitecture getArchitectureFromPayloadsAndOrbits (ArrayList<ArrayList<String>> currentPayloads, ArrayList<String> currentOrbits) {
        // ORDER ORBITS AND INSTRUMENTS BASED ON CLIMATECENTRICPARAMS
        ArrayList<String> instrumentList = new ArrayList<>(Arrays.asList(params.getInstrumentList()));
        ArrayList<String> orbitList = new ArrayList<>(Arrays.asList(params.getOrbitList()));

        PartitioningArchitecture arch = new PartitioningArchitecture(instrumentList.size(), orbitList.size(), 2);

        for (int i = 0; i < currentOrbits.size(); i++) {
            ArrayList<String> currentOrbitPayloads = currentPayloads.get(i);

            for (int j = 0; j < currentOrbitPayloads.size(); j++) {
                int payloadIndex = instrumentList.indexOf(currentOrbitPayloads.get(j));

                IntegerVariable instrVar = new IntegerVariable(i, 0, instrumentList.size()-1);
                arch.setVariable(payloadIndex, instrVar);
            }

            int orbitIndex = orbitList.indexOf(currentOrbits.get(i));

            IntegerVariable orbitVar = new IntegerVariable(orbitIndex, 0, instrumentList.size()-1);
            arch.setVariable(instrumentList.size()+i, orbitVar);
        }

        IntegerVariable noOrbitVar = new IntegerVariable(-1, -1, -1);
        for (int k = 0; k < instrumentList.size()-currentOrbits.size(); k++) {
            arch.setVariable(instrumentList.size()+currentOrbits.size()+k, noOrbitVar);
        }

        return arch;
    }
}
