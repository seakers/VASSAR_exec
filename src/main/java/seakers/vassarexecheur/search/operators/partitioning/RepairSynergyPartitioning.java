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
 * Moves an instrument from one orbit to another if instruments of a synergistic pair are present in
 * different orbits
 *
 * @author roshansuresh
 */

public class RepairSynergyPartitioning implements Variation {

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

    public RepairSynergyPartitioning(int numChanges, ResourcePool resourcePool, ArchitectureEvaluator evaluator, BaseParams params) {
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

        HashMap<String, String[]> instrumentSynergyMap = getInstrumentSynergyNameMap(params);

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

        ArrayList<ArrayList<Integer>> possibleInstrumentMoves = getValidMoves(rete, satellites, instrumentSynergyMap);

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

    private ArrayList<ArrayList<Integer>> getValidMoves (Rete r, ArrayList<Fact> allSatellites, HashMap<String, String[]> synergyMap) {
        // Valid move -> ArrayList {instrument_index (in satellite_to_remove_from), index_of_satellite_to_remove_from, index_of_satellite_to_add_to}

        ArrayList<ArrayList<Integer>> validMoves = new ArrayList<>();
        ArrayList<String> interferenceMapKeys = new ArrayList<String>(synergyMap.keySet());

        for (int i = 0; i < allSatellites.size(); i++) {
            Fact currentSatellite = allSatellites.get(i);
            try {
                ValueVector satelliteInstruments = currentSatellite.getSlotValue("instruments").listValue(r.getGlobalContext());

                for (int j = 0; j < satelliteInstruments.size(); j++) {
                    String currentInstrument  = satelliteInstruments.get(j).stringValue(r.getGlobalContext());

                    if (synergyMap.containsKey(currentInstrument)) {
                        String[] synergyInstruments = synergyMap.get(currentInstrument);
                        for (int k = 0; k < allSatellites.size(); k++) {
                            if (k == i) {
                                continue;
                            }
                            ValueVector otherSatelliteInstruments = allSatellites.get(k).getSlotValue("instruments").listValue(r.getGlobalContext());
                            ArrayList<Integer> presentValueIntegers = checkInstrumentsinSatellite(r, satelliteInstruments, synergyInstruments);
                            if (presentValueIntegers.size() > 0) {
                                for (int m = 0; m < presentValueIntegers.size(); m++) {
                                    ArrayList<Integer> shiftingChoice = new ArrayList<>();
                                    shiftingChoice.add(presentValueIntegers.get(m));
                                    shiftingChoice.add(k);
                                    shiftingChoice.add(i);
                                    validMoves.addAll(Collections.singleton(shiftingChoice));
                                }
                                ArrayList<Integer> shiftingChoice = new ArrayList<>();
                                shiftingChoice.add(j);
                                shiftingChoice.add(i);
                                shiftingChoice.add(k);
                                validMoves.addAll(Collections.singleton(shiftingChoice));
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
     * Creates instrument synergy map used to compute the instrument synergy violation heuristic (only formulated for the
     * Climate Centric problem for now) (added by roshansuresh)
     * @param params
     * @return Instrument synergy hashmap
     */
    protected HashMap<String, String[]> getInstrumentSynergyNameMap(BaseParams params) {
        HashMap<String, String[]> synergyNameMap = new HashMap<>();
        if (params.getProblemName().equalsIgnoreCase("ClimateCentric")) {
            synergyNameMap.put("ACE_ORCA", new String[]{"DESD_LID", "GACM_VIS", "ACE_POL", "HYSP_TIR", "ACE_LID"});
            synergyNameMap.put("DESD_LID", new String[]{"ACE_ORCA", "ACE_LID", "ACE_POL"});
            synergyNameMap.put("GACM_VIS", new String[]{"ACE_ORCA", "ACE_LID"});
            synergyNameMap.put("HYSP_TIR", new String[]{"ACE_ORCA", "POSTEPS_IRS"});
            synergyNameMap.put("ACE_POL", new String[]{"ACE_ORCA", "DESD_LID"});
            synergyNameMap.put("ACE_LID", new String[]{"ACE_ORCA", "CNES_KaRIN", "DESD_LID", "GACM_VIS"});
            synergyNameMap.put("POSTEPS_IRS", new String[]{"HYSP_TIR"});
            synergyNameMap.put("CNES_KaRIN", new String[]{"ACE_LID"});
        }
        else {
            System.out.println("Synergy Map for current problem not formulated");
        }
        return synergyNameMap;
    }

    private ArrayList<ArrayList<String>> getSatellitePayloadsFromSatelliteFacts (Rete r, ArrayList<Fact> allSatellites) throws JessException {
        ArrayList<ArrayList<String>> satellitePayloads = new ArrayList<>();
        for (int i = 0; i < allSatellites.size(); i++) {
            ArrayList<String> currentSatellitePayload = new ArrayList<>();
            ValueVector instrumentsString = allSatellites.get(i).getSlotValue("instruments").listValue(r.getGlobalContext());
            for (int j = 0; j < instrumentsString.size(); j++) {
                currentSatellitePayload.add(instrumentsString.get(j).stringValue(r.getGlobalContext()));
            }
            satellitePayloads.add(currentSatellitePayload);
        }
        return satellitePayloads;
    }

    private ArrayList<String> getSatelliteOrbitsFromSatelliteFacts (Rete r, ArrayList<Fact> allSatellites) throws JessException {
        ArrayList<String> satelliteOrbits = new ArrayList<>();

        for (int i = 0; i < allSatellites.size(); i++) {
            satelliteOrbits.add(allSatellites.get(i).getSlotValue("orbit-string").stringValue(r.getGlobalContext()));
        }

        return satelliteOrbits;
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
