package seakers.vassarexecheur.search.operators.assigning;

import jess.Fact;
import jess.JessException;
import jess.Rete;
import jess.ValueVector;
import org.moeaframework.core.PRNG;
import org.moeaframework.core.Solution;
import org.moeaframework.core.Variation;
import org.moeaframework.core.variable.BinaryVariable;
import org.moeaframework.core.variable.RealVariable;
import seakers.vassarexecheur.search.problems.assigning.AssigningArchitecture;
import seakers.vassarexecheur.search.problems.assigning.AssigningProblem;
import seakers.vassarheur.BaseParams;
import seakers.vassarheur.QueryBuilder;
import seakers.vassarheur.Resource;
import seakers.vassarheur.ResourcePool;
import seakers.vassarheur.architecture.AbstractArchitecture;
import seakers.vassarheur.problems.Assigning.ArchitectureEvaluator;
import seakers.vassarheur.utils.MatlabFunctions;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.stream.Collectors;

/**
 * Adds a synergistic instrument to one or more orbits to another if they are absent in those orbits
 *
 * @author roshansuresh
 */

public class RepairSynergyAdditionAssigning implements Variation {

    /**
     * Number of instrument additions
     */
    private int numberOfChanges;

    /**
     * Resource pool to obtain matlab functions
     */
    private final ResourcePool resourcePool;

    private final ArchitectureEvaluator evaluator;

    //private final ParallelPRNG pprng;

    private final BaseParams params;

    /**
     * Assigning Problem used to evaluate architectures
     */
    private AssigningProblem problem;

    private final HashMap<String, String[]> synergyMap;

    public RepairSynergyAdditionAssigning(int numChanges, ResourcePool resourcePool, ArchitectureEvaluator evaluator, BaseParams params, AssigningProblem problem, HashMap<String, String[]> synergyMap) {
        this.numberOfChanges = numChanges;
        this.resourcePool = resourcePool;
        this.evaluator = evaluator;
        this.params = params;
        this.problem = problem;
        this.synergyMap = synergyMap;
    }

    @Override
    public int getArity() {
        return 1;
    }

    @Override
    public Solution[] evolve(Solution[] solutions) {
        AssigningArchitecture parent = (AssigningArchitecture) solutions[0];
        //System.out.println("parent");
        //System.out.println(parent.getBitString());
        AbstractArchitecture arch_abs = problem.getAbstractArchitecture(parent);

        Resource res = resourcePool.getResource();
        MatlabFunctions m = res.getM();
        Rete rete = res.getRete();
        QueryBuilder queryBuilder = res.getQueryBuilder();

        //HashMap<String, String[]> instrumentSynergyMap = getInstrumentSynergyNameMap();

        evaluator.assertMissions(params, rete, arch_abs, m);
        ArrayList<Fact> satellites = queryBuilder.makeQuery("MANIFEST::Satellite");

        ArrayList<ArrayList<String>> payloads = new ArrayList<>();
        ArrayList<String> orbits = new ArrayList<>();
        try {
            payloads = getSatellitePayloadsFromSatelliteFacts(rete, satellites);
            orbits = getSatelliteOrbitsFromSatelliteFacts(rete, satellites);
        } catch (JessException e) {
            e.printStackTrace();
        }

        ArrayList<ArrayList<Integer>> possibleAdditionSatellites = getValidSatelliteSynergies(rete, satellites, synergyMap);

        // Make choices of instrument addition randomly
        int numberOfMoves = 0;
        while ((numberOfMoves < numberOfChanges) && (possibleAdditionSatellites.size() > 0)) {
            int satChoiceIndex = PRNG.nextInt(possibleAdditionSatellites.size());
            ArrayList<Integer> satChoice = possibleAdditionSatellites.get(satChoiceIndex);

            try {
                ValueVector satelliteInstruments = satellites.get(satChoice.get(1)).getSlotValue("instruments").listValue(rete.getGlobalContext());
                String synergyInstrument  = satelliteInstruments.get(satChoice.get(0)).stringValue(rete.getGlobalContext());
                String[] synergyInstruments = synergyMap.get(synergyInstrument);

                int synergisticInstrumentChoiceIndex = PRNG.nextInt(synergyInstruments.length);
                payloads.get(satChoice.get(1)).add(synergyInstruments[synergisticInstrumentChoiceIndex]);
                numberOfMoves++;
                possibleAdditionSatellites.remove(satChoice);
            } catch (JessException e) {
                e.printStackTrace();
            }
        }
        this.resourcePool.freeResource(res);
        AssigningArchitecture child = getArchitectureFromPayloadsAndOrbits(payloads, orbits);

        return new Solution[]{child};
    }

    private ArrayList<ArrayList<Integer>> getValidSatelliteSynergies (Rete r, ArrayList<Fact> allSatellites, HashMap<String, String[]> synergyMap) {
        // Valid move -> ArrayList {instrument_index (in satellite_to_add_to), index_of_satellite_to_add_to}

        ArrayList<ArrayList<Integer>> validSatelliteInstruments = new ArrayList<>();
        //ArrayList<String> synergyMapKeys = new ArrayList<String>(synergyMap.keySet());

        for (int i = 0; i < allSatellites.size(); i++) {
            Fact currentSatellite = allSatellites.get(i);
            try {
                ValueVector satelliteInstruments = currentSatellite.getSlotValue("instruments").listValue(r.getGlobalContext());

                for (int j = 0; j < satelliteInstruments.size(); j++) {
                    String currentInstrument  = satelliteInstruments.get(j).stringValue(r.getGlobalContext());

                    if (synergyMap.containsKey(currentInstrument)) {
                        String[] synergyInstruments = synergyMap.get(currentInstrument);

                        ArrayList<Integer> synergyInstrumentIndices = checkInstrumentsInSatellite(r, satelliteInstruments, synergyInstruments);
                        if (synergyInstrumentIndices.isEmpty()) {
                            ArrayList<Integer> shiftingChoice = new ArrayList<>();
                            shiftingChoice.add(j);
                            shiftingChoice.add(i);
                            validSatelliteInstruments.addAll(Collections.singleton(shiftingChoice));
                        }
                    }
                }
            } catch (JessException e) {
                e.printStackTrace();
            }
        }
        return validSatelliteInstruments;
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

    private int checkInstrumentInSatellite(Rete r, ValueVector instrumentsInSatellite, String checkInstrument) throws JessException {
        int presentSatelliteIndex = -1; // index in instrumentsInSatellite
        for (int j = 0; j < instrumentsInSatellite.size(); j++) {
            if (checkInstrument.equalsIgnoreCase(instrumentsInSatellite.get(j).stringValue(r.getGlobalContext()))) {
                presentSatelliteIndex = j;
                break;
            }
        }
        return presentSatelliteIndex;
    }

    private ArrayList<Integer> checkInstrumentsInSatellite(Rete r, ValueVector instrumentsInSatellite, String[] checkInstruments) throws JessException {
        ArrayList<Integer> presentSatelliteIndices = new ArrayList<>(); // indices in instrumentsInSatellite
        for (int i = 0; i < checkInstruments.length; i++) {
            int instrumentIndex = checkInstrumentInSatellite(r, instrumentsInSatellite, checkInstruments[i]);
            if (instrumentIndex != -1) {
                presentSatelliteIndices.add(instrumentIndex);
            }
        }
        return presentSatelliteIndices;
    }

    private AssigningArchitecture getArchitectureFromPayloadsAndOrbits(ArrayList<ArrayList<String>> currentPayloads, ArrayList<String> currentOrbits) {

        ArrayList<String> instrumentList = new ArrayList<>(Arrays.asList(params.getInstrumentList()));
        ArrayList<String> orbitList = new ArrayList<>(Arrays.asList(params.getOrbitList()));

        AssigningArchitecture arch = new AssigningArchitecture(new int[]{1}, params.getNumInstr(), params.getNumOrbits(), 2);

        RealVariable var0 = new RealVariable(0.0, 0.0, 0.0);
        arch.setVariable(0, var0);
        for (int i = 0; i < params.getNumOrbits(); i++) {
            if (!currentOrbits.contains(orbitList.get(i))) {
                for (int j = 0; j < params.getNumInstr(); j++) {
                    BinaryVariable var = new BinaryVariable(1);
                    var.set(0, false);
                    int decisionIndex = params.getNumInstr()*i + j + 1;
                    arch.setVariable(decisionIndex, var);
                }
            } else {
                int orbitIndex = currentOrbits.indexOf(orbitList.get(i));
                ArrayList<String> currentInstruments = currentPayloads.get(orbitIndex);
                for (int j = 0; j < params.getNumInstr(); j++) {
                    BinaryVariable var = new BinaryVariable(1);
                    int decisionIndex = params.getNumInstr()*i + j + 1;
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
