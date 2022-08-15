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
import seakers.vassarheur.Resource;
import seakers.vassarheur.ResourcePool;
import seakers.vassarheur.architecture.AbstractArchitecture;
import seakers.vassarheur.problems.Assigning.Architecture;
import seakers.vassarheur.problems.Assigning.ArchitectureEvaluator;
import seakers.vassarheur.problems.Assigning.ClimateCentricAssigningParams;
import seakers.vassarheur.utils.MatlabFunctions;

import java.lang.reflect.Array;
import java.util.*;

/**
 * Moves an instrument from one orbit to another if instruments of a synergistic pair are present in
 * different orbits
 *
 * @author roshansuresh
 */

public class RepairSynergy implements Variation {
    /**
     * Rete object
     */
    private final Rete rete;

    /**
     * Query Builder object to extract slot data
     */
    private final QueryBuilder queryBuilder;

    /**
     * Number of instrument moves
     */
    private int numberOfChanges;

    /**
     * Resource pool to obtain matlab functions
     */
    private final ResourcePool resourcePool;

    private final ArchitectureEvaluator evaluator;

    private final ParallelPRNG pprng;

    private final BaseParams params;

    public RepairSynergy(int numChanges, ResourcePool resourcePool, ArchitectureEvaluator evaluator, Rete r, QueryBuilder qb, BaseParams params) {
        this.numberOfChanges = numChanges;
        this.resourcePool = resourcePool;
        this.evaluator = evaluator;
        this.rete = r;
        this.queryBuilder = qb;
        this.params = params;
        this.pprng = new ParallelPRNG();
    }

    @Override
    public int getArity() {
        return 1 ;
    }

    @Override
    public Solution[] evolve(Solution[] sols) {
        AssigningArchitecture parent = (AssigningArchitecture) sols[0];
        StringBuilder bitStringBuilder = new StringBuilder(parent.getNumberOfVariables());
        for (int i = 1; i < parent.getNumberOfVariables(); ++i) {
            bitStringBuilder.append(parent.getVariable(i).toString());
        }

        AbstractArchitecture arch_abs;
        arch_abs = new Architecture(bitStringBuilder.toString(), ((ClimateCentricAssigningParams) params).getNumSatellites()[0], (ClimateCentricAssigningParams) params);

        Resource res = resourcePool.getResource();
        MatlabFunctions m = res.getM();

        HashMap<String, String[]> instrumentSynergyMap = getInstrumentSynergyNameMap(params);

        evaluator.assertMissions(params, rete, arch_abs, m);
        ArrayList<Fact> satellites = queryBuilder.makeQuery("MANIFEST::Satellite");

        ArrayList<ArrayList<String>> payloads = new ArrayList<>();
        ArrayList<String> orbits = new ArrayList<>();
        try {
            payloads = getSatellitePayloadsFromSatelliteFacts(satellites);
            orbits = getSatelliteOrbitsFromSatelliteFacts(satellites);
        } catch (JessException e) {
            e.printStackTrace();
        }

        ArrayList<ArrayList<Integer>> possibleInstrumentMoves = getValidMoves(satellites, instrumentSynergyMap);

        // Make choices of instrument move randomly
        int numberOfMoves = 0;
        while ((numberOfMoves < numberOfChanges) && (possibleInstrumentMoves.size() > 0)) {
            int moveChoiceIndex = pprng.nextInt(possibleInstrumentMoves.size());
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

        AssigningArchitecture child = getArchitectureFromPayloadsAndOrbits(payloads, orbits);

        return new Solution[]{child};
    }

    private ArrayList<ArrayList<Integer>> getValidMoves (ArrayList<Fact> allSatellites, HashMap<String, String[]> synergyMap) {
        // Valid move -> ArrayList {instrument_index (in satellite_to_remove_from), index_of_satellite_to_remove_from, index_of_satellite_to_add_to}

        ArrayList<ArrayList<Integer>> validMoves = new ArrayList<>();
        ArrayList<String> interferenceMapKeys = new ArrayList<String>(synergyMap.keySet());

        for (int i = 0; i < allSatellites.size(); i++) {
            Fact currentSatellite = allSatellites.get(i);
            try {
                ValueVector satelliteInstruments = currentSatellite.getSlotValue("instruments").listValue(rete.getGlobalContext());

                for (int j = 0; j < satelliteInstruments.size(); j++) {
                    String currentInstrument  = satelliteInstruments.get(j).stringValue(rete.getGlobalContext());

                    if (synergyMap.containsKey(currentInstrument)) {
                        String[] synergyInstruments = synergyMap.get(currentInstrument);
                        for (int k = 0; k < allSatellites.size(); k++) {
                            if (k == i) {
                                continue;
                            }
                            ValueVector otherSatelliteInstruments = allSatellites.get(k).getSlotValue("instruments").listValue(rete.getGlobalContext());
                            ArrayList<Integer> presentValueIntegers = checkInstrumentsinSatellite(satelliteInstruments, synergyInstruments);
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

    private AssigningArchitecture getArchitectureFromPayloadsAndOrbits(ArrayList<ArrayList<String>> currentPayloads, ArrayList<String> currentOrbits) {

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

    private ArrayList<ArrayList<String>> getSatellitePayloadsFromSatelliteFacts (ArrayList<Fact> allSatellites) throws JessException {
        ArrayList<ArrayList<String>> satellitePayloads = new ArrayList<>();
        for (int i = 0; i < allSatellites.size(); i++) {
            ArrayList<String> currentSatellitePayload = new ArrayList<>();
            ValueVector instrumentsString = allSatellites.get(i).getSlotValue("instruments").listValue(rete.getGlobalContext());
            for (int j = 0; j < instrumentsString.size(); j++) {
                currentSatellitePayload.add(instrumentsString.get(j).stringValue(rete.getGlobalContext()));
            }
            satellitePayloads.add(currentSatellitePayload);
        }
        return satellitePayloads;
    }

    private ArrayList<String> getSatelliteOrbitsFromSatelliteFacts (ArrayList<Fact> allSatellites) throws JessException {
        ArrayList<String> satelliteOrbits = new ArrayList<>();

        for (int i = 0; i < allSatellites.size(); i++) {
            satelliteOrbits.add(allSatellites.get(i).getSlotValue("orbit-string").stringValue(rete.getGlobalContext()));
        }

        return satelliteOrbits;
    }

    private ArrayList<Integer> checkInstrumentsinSatellite(ValueVector instrumentsInSatellite, String[] checkInstruments) throws JessException {
        ArrayList<Integer> presentSatelliteIndices = new ArrayList<>(); // indices in instrumentsInSatellite
        for (int i = 0; i < checkInstruments.length; i++) {
            for (int j = 0; j < instrumentsInSatellite.size(); j++) {
                if (checkInstruments[i].equalsIgnoreCase(instrumentsInSatellite.get(j).stringValue(rete.getGlobalContext()))) {
                    presentSatelliteIndices.add(j);
                    break;
                }
            }
        }
        return presentSatelliteIndices;
    }
}
