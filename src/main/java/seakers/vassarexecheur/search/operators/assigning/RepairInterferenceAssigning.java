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
 * Moves an instrument from one satellite to another if its interfering pair is present in the same satellite.
 * Interfering instrument pairs are predetermined for the set of candidate instruments.
 *
 * @author roshansuresh
 */

public class RepairInterferenceAssigning implements Variation{

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

    public RepairInterferenceAssigning(int numChanges, ResourcePool resourcePool, ArchitectureEvaluator evaluator, BaseParams params) {
        this.numberOfChanges = numChanges;
        this.resourcePool = resourcePool;
        this.evaluator = evaluator;
        this.params = params;
        //this.pprng = new ParallelPRNG();
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

        AssigningArchitecture child = getArchitectureFromPayloadsAndOrbits(payloads, orbits);

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

    private ArrayList<Integer> getExclusiveSatelliteMoves (int numberOfSatellites, int excludedSatellite) {
        ArrayList<Integer> possibleMoves = new ArrayList<>();
        for (int i = 0; i < numberOfSatellites; i++) {
            if (i != excludedSatellite) {
                possibleMoves.add(i);
            }
        }
        return possibleMoves;
    }

    @Override
    public int getArity() {
        return 1;
    }


}
