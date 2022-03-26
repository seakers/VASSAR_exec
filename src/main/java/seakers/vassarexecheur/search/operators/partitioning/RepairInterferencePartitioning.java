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
import seakers.vassarexecheur.search.problems.partitioning.PartitioningArchitecture;
import seakers.vassarheur.BaseParams;
import seakers.vassarheur.QueryBuilder;
import seakers.vassarexecheur.search.problems.assigning.AssigningArchitecture;
import seakers.vassarheur.Resource;
import seakers.vassarheur.ResourcePool;
import seakers.vassarheur.architecture.AbstractArchitecture;
import seakers.vassarheur.problems.PartitioningAndAssigning.Architecture;
import seakers.vassarheur.problems.PartitioningAndAssigning.ArchitectureEvaluator;
import seakers.vassarheur.problems.PartitioningAndAssigning.ClimateCentricPartitioningParams;
import seakers.vassarheur.utils.MatlabFunctions;

import java.lang.reflect.Array;
import java.util.*;

/**
 * Moves an instrument from one satellite to another if its interfering pair is present in the same satellite.
 * Interfering instrument pairs are predetermined for the set of candidate instruments.
 *
 * @author roshansuresh
 */

public class RepairInterferencePartitioning implements Variation {

    /**
     * Number of instrument moves
     */
    private int numberOfChanges;

    /**
     * Resource pool to obtain matlab functions
     */
    private final ResourcePool resourcePool;

    private final ArchitectureEvaluator evaluator;

    //private final ParallelPRNG pprng;

    private final BaseParams params;

    public RepairInterferencePartitioning(int numChanges, ResourcePool resourcePool, ArchitectureEvaluator evaluator, BaseParams params) {
        this.numberOfChanges = numChanges;
        this.resourcePool = resourcePool;
        this.evaluator = evaluator;
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

        int numPartitioningVariables = params.getNumInstr();
        int numAssignmentVariables = params.getNumInstr();

        int[] instrumentPartitioning = new int[numPartitioningVariables];
        int[] orbitAssignment = new int[numAssignmentVariables];

        for (int i = 0; i < numPartitioningVariables; i++) {
            instrumentPartitioning[i] = ((IntegerVariable) parent.getVariable(i)).getValue();
        }

        for (int i = 0; i < numAssignmentVariables; i++) {
            orbitAssignment[i] = ((IntegerVariable) parent.getVariable(numPartitioningVariables + i)).getValue();
        }

        // Check constraint
        double constraint = 1.0;
        if (!isFeasible(instrumentPartitioning, orbitAssignment)) {
            constraint = 0.0;
        }
        parent.setConstraint(0, constraint);

        AbstractArchitecture arch_abs = new Architecture(instrumentPartitioning, orbitAssignment, 1, params);

        Resource res = resourcePool.getResource();
        MatlabFunctions m = res.getM();
        Rete rete = res.getRete();
        QueryBuilder queryBuilder = res.getQueryBuilder();

        HashMap<String, String[]> interferingInstrumentsMap = getInstrumentInterferenceNameMap(params);

        evaluator.assertMissions(params, rete, arch_abs, m);
        ArrayList<Fact> satellites = queryBuilder.makeQuery("MANIFEST::Satellite");

        //// METHOD 1 - Get payloads and orbits from Satellite facts
        //ArrayList<ArrayList<String>> payloads = new ArrayList<>();
        //ArrayList<String> orbits = new ArrayList<>();
        //try {
        //payloads = getSatellitePayloadsFromSatelliteFacts(satellites);
        //orbits = getSatelliteOrbitsFromSatelliteFacts(satellites);
        //} catch (JessException e) {
        //e.printStackTrace();
        //}

        //// METHOD 2 - Get payloads and satellites from inherent methods
        ArrayList<ArrayList<String>> payloads = parent.getSatellitePayloads();
        ArrayList<String> orbits = parent.getSatelliteOrbits();

        ArrayList<ArrayList<Integer>> possibleInstrumentMoves = getValidMoves(rete, satellites, interferingInstrumentsMap);

        // Make choices of instrument move randomly
        int numberOfMoves = 0;
        while ((numberOfMoves < numberOfChanges) && (possibleInstrumentMoves.size() > 0)) {
            int moveChoiceIndex = PRNG.nextInt(possibleInstrumentMoves.size());
            ArrayList<Integer> moveChoice = possibleInstrumentMoves.get(moveChoiceIndex);

            // Update payloads
            ArrayList<String> removalOrbitPayload = payloads.get(moveChoice.get(1));
            String instrumentToMove = removalOrbitPayload.get(moveChoice.get(0));
            removalOrbitPayload.remove(moveChoice.get(0));

            ArrayList<String> additionOrbitPayload = payloads.get(moveChoice.get(2));
            additionOrbitPayload.add(instrumentToMove);

            payloads.set(moveChoice.get(1), removalOrbitPayload);
            payloads.set(moveChoice.get(2), additionOrbitPayload);

            possibleInstrumentMoves.remove(moveChoiceIndex);
            numberOfMoves += 1;
        }
        PartitioningArchitecture child = getArchitectureFromPayloadsAndOrbits(payloads, orbits);

        return new Solution[]{child};
    }

    private ArrayList<ArrayList<Integer>> getValidMoves (Rete r, ArrayList<Fact> allSatellites, HashMap<String, String[]> interferenceMap) {
        // Valid move -> ArrayList {instrument_index (in satellite_to_remove_from), index_of_satellite_to_remove_from, index_of_satellite_to_add_to}

        ArrayList<ArrayList<Integer>> validMoves = new ArrayList<>();
        ArrayList<String> interferenceMapKeys = new ArrayList<String>(interferenceMap.keySet());

        for (int i = 0; i < allSatellites.size(); i++) {
            Fact currentSatellite = allSatellites.get(i);
            try {
                ValueVector satelliteInstruments = currentSatellite.getSlotValue("instruments").listValue(r.getGlobalContext());

                for (int j = 0; j < satelliteInstruments.size(); j++) {
                    String currentInstrument  = satelliteInstruments.get(j).stringValue(r.getGlobalContext());

                    if (interferenceMapKeys.contains(currentInstrument)) { // If an instrument in the interference Map keyset is present in the satellite
                        String[] interferenceMapValues = interferenceMap.get(currentInstrument);
                        ArrayList<Integer> presentValueIntegers = checkInstrumentsinSatellite(r, satelliteInstruments, interferenceMapValues);
                        if (presentValueIntegers.size() > 0) {
                            presentValueIntegers.add(j);
                            ArrayList<Integer> possibleSatelliteMoves = getExclusiveSatelliteMoves(allSatellites.size(), i);
                            for (int k = 0; k < presentValueIntegers.size(); k++) {
                                for (int m = 0; m < possibleSatelliteMoves.size(); m++) {
                                    ArrayList<Integer> removalChoice = new ArrayList<>();
                                    removalChoice.add(presentValueIntegers.get(k));
                                    removalChoice.add(i);
                                    removalChoice.add(possibleSatelliteMoves.get(m));
                                    validMoves.addAll(Collections.singleton(removalChoice));
                                }
                            }
                        }
                    }
                }
            } catch (JessException e) {
                e.printStackTrace();
            }
        }
        return validMoves;
    }

    private ArrayList<Integer> checkInstrumentsinSatellite(Rete r, ValueVector instrumentsInSatellite, String[] checkInstruments) throws JessException {
        ArrayList<Integer> presentSatelliteIndices = new ArrayList<>(); // indices in instrumentsInSatellite
        for (int i = 0; i < checkInstruments.length; i++) {
            for (int j = 0; j < instrumentsInSatellite.size(); j++) {
                if (checkInstruments[i].equalsIgnoreCase(instrumentsInSatellite.get(j).stringValue(r.getGlobalContext()))) {
                    presentSatelliteIndices.add(j);
                    break;
                }
            }
        }
        return presentSatelliteIndices;
    }

    private ArrayList<Integer> getExclusiveSatelliteMoves (int numberOfSatellites, int excludedSatellite) {
        ArrayList<Integer> possibleMoves = new ArrayList<>();
        for (int i = 0; i < numberOfSatellites; i++) {
            if (i != excludedSatellite) {
                possibleMoves.add(i);
            }
        }
        return possibleMoves;
    }

    private boolean isFeasible(int[] instrumentPartitioning, int[] orbitAssignment){

        // Check if the number of satellites matches the number of orbit assignments
        Set<Integer> satIndices = new HashSet<>();
        for (Integer instrumentPartition: instrumentPartitioning) {
            satIndices.add(instrumentPartition);
        }

        for (int i = 0; i < orbitAssignment.length; i++) {
            if (orbitAssignment[i] < 0) {
                if (satIndices.size() != i) {
                    return false;
                }
            }
        }

        // Check if the index of the new satellite is +1 of the largest index
        int max = 0;
        for (Integer instrumentPartition: instrumentPartitioning) {
            if (instrumentPartition > max) {
                if (instrumentPartition == max + 1) {
                    max = instrumentPartition;
                }
                else {
                    return false;
                }
            }
        }

        return true;
    }

    /**
     * Creates instrument interference map used to compute the instrument interference violation heuristic (only formulated for the
     * Climate Centric problem for now)
     * @param params
     * @return Instrument interference hashmap
     */
    protected HashMap<String, String[]> getInstrumentInterferenceNameMap(BaseParams params) {
        HashMap<String, String[]> interferenceNameMap = new HashMap<>();
        if (params.getProblemName().equalsIgnoreCase("ClimateCentric")) {
            interferenceNameMap.put("ACE_LID", new String[]{"ACE_CPR", "DESD_SAR", "CLAR_ERB", "GACM_SWIR"});
            interferenceNameMap.put("ACE_CPR", new String[]{"ACE_LID", "DESD_SAR", "CNES_KaRIN", "CLAR_ERB", "ACE_POL", "ACE_ORCA", "GACM_SWIR"});
            interferenceNameMap.put("DESD_SAR", new String[]{"ACE_LID", "ACE_CPR"});
            interferenceNameMap.put("CLAR_ERB", new String[]{"ACE_LID", "ACE_CPR"});
            interferenceNameMap.put("CNES_KaRIN", new String[]{"ACE_CPR"});
            interferenceNameMap.put("ACE_POL", new String[]{"ACE_CPR"});
            interferenceNameMap.put("ACE_ORCA", new String[]{"ACE_CPR"});
            interferenceNameMap.put("GACM_SWIR", new String[]{"ACE_LID", "ACE_CPR"});
        }
        else {
            System.out.println("Interference Map fpr current problem not formulated");
        }
        return interferenceNameMap;
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
