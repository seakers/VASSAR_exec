package seakers.vassarexecheur.search.operators.assigning;

import jess.*;
import org.moeaframework.core.PRNG;
import org.moeaframework.core.Solution;
import org.moeaframework.core.Variation;
import org.moeaframework.core.variable.BinaryVariable;
import org.moeaframework.core.variable.RealVariable;
import seakers.vassarexecheur.search.problems.assigning.AssigningProblem;
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
import java.util.stream.Collectors;

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

    /**
     * Toggles movement of removed instruments into other spacecraft
     */
    private boolean moveInstruments;

    /**
     * Assigning Problem used to evaluate architectures
     */
    private AssigningProblem problem;

    private final HashMap<String, String[]> interferenceMap;

    public RepairInterferenceAssigning(int numChanges, ResourcePool resourcePool, ArchitectureEvaluator evaluator, BaseParams params, AssigningProblem problem, HashMap<String, String[]> interferenceMap) {
        this.numberOfChanges = numChanges;
        this.resourcePool = resourcePool;
        this.evaluator = evaluator;
        this.params = params;
        this.moveInstruments = true;
        this.problem = problem;
        this.interferenceMap = interferenceMap;
        //this.pprng = new ParallelPRNG();
    }

    public RepairInterferenceAssigning(int numChanges, ResourcePool resourcePool, ArchitectureEvaluator evaluator, BaseParams params, AssigningProblem problem, HashMap<String, String[]> interferenceMap, boolean moveInstruments) {
        this.numberOfChanges = numChanges;
        this.resourcePool = resourcePool;
        this.evaluator = evaluator;
        this.params = params;
        this.moveInstruments = moveInstruments;
        this.problem = problem;
        this.interferenceMap = interferenceMap;
        //this.pprng = new ParallelPRNG();
    }

    @Override
    public Solution[] evolve(Solution[] sols) {
        AssigningArchitecture parent = (AssigningArchitecture) sols[0];
        //System.out.println("parent");
        //System.out.println(parent.getBitString());
        AbstractArchitecture arch_abs = problem.getAbstractArchitecture(parent);

        Resource res = resourcePool.getResource();
        MatlabFunctions m = res.getM();
        Rete rete = res.getRete();
        QueryBuilder queryBuilder = res.getQueryBuilder();

        //HashMap<String, String[]> interferingInstrumentsMap = getInstrumentInterferenceNameMap(params);

        evaluator.assertMissions(params, rete, arch_abs, m);
        ArrayList<Fact> satellites = queryBuilder.makeQuery("MANIFEST::Satellite");

        //// METHOD 1 - Get payloads and orbits from Satellite facts
        ArrayList<ArrayList<String>> payloads = new ArrayList<>();
        ArrayList<String> orbits = new ArrayList<>();
        try {
            payloads = getSatellitePayloadsFromSatelliteFacts(rete, satellites);
            orbits = getSatelliteOrbitsFromSatelliteFacts(rete, satellites);
        } catch (JessException e) {
            e.printStackTrace();
        }

        //// METHOD 2 - Get payloads and satellites from inherent methods
        //ArrayList<ArrayList<String>> payloads = parent.getSatellitePayloads();
        //ArrayList<String> orbits = parent.getSatelliteOrbits();

        ArrayList<ArrayList<Integer>> possibleInstrumentMoves = getValidMoves(rete, satellites, interferenceMap);
        ArrayList<ArrayList<Integer>> possibleUniqueInstrumentMoves = (ArrayList<ArrayList<Integer>>) possibleInstrumentMoves.stream().distinct().collect(Collectors.toList());

        // Make choices of instrument move randomly
        int numberOfMoves = 0;
        while ((numberOfMoves < numberOfChanges) && (possibleUniqueInstrumentMoves.size() > 0)) {
            //System.out.println("Move Count: " + numberOfMoves);
            int numberOfTries = 0;
            boolean feasibleMove = false;
            while ((!feasibleMove) && (numberOfTries < 1) && (possibleUniqueInstrumentMoves.size() > 0)) {
                int moveChoiceIndex = PRNG.nextInt(possibleUniqueInstrumentMoves.size());
                ArrayList<Integer> moveChoice = possibleUniqueInstrumentMoves.get(moveChoiceIndex);
                //System.out.println("Selected Move: " + moveChoice);

                // Update payloads
                ArrayList<String> removalOrbitPayload = payloads.get(moveChoice.get(1));
                String instrumentToMove = removalOrbitPayload.get(moveChoice.get(0));
                removalOrbitPayload.remove(instrumentToMove);

                payloads.set(moveChoice.get(1), removalOrbitPayload);

                ArrayList<String> additionOrbitPayload;
                if (moveInstruments) {
                    additionOrbitPayload = payloads.get(moveChoice.get(2));
                    additionOrbitPayload.add(instrumentToMove);

                    payloads.set(moveChoice.get(2), additionOrbitPayload);
                }

                // Check if instrument movement improves heuristic satisfaction
                AssigningArchitecture candidateChild = getArchitectureFromPayloadsAndOrbits(payloads, orbits);
                AbstractArchitecture candidateChild_abs = problem.getAbstractArchitecture(candidateChild);

                try {
                    rete.reset();
                } catch (JessException e) {
                    e.printStackTrace();
                }
                evaluator.assertMissions(params, rete, candidateChild_abs, m);
                ArrayList<ArrayList<Double>> satelliteHeuristics;
                ArrayList<Double> archHeuristics = null;
                try {
                    evaluator.evaluateHeuristicParameters(rete, candidateChild_abs, queryBuilder, m);
                    satelliteHeuristics = evaluator.computeHeuristics(rete, candidateChild_abs, queryBuilder, params);
                    archHeuristics = evaluator.computeHeuristicsArchitecture(satelliteHeuristics);
                } catch (JessException e) {
                    e.printStackTrace();
                }

                // archHeuristics -> [dutyCycleViolation, instrumentOrbitAssignmentViolation, interferenceViolation, packingEfficiencyViolation, spacecraftMassViolation, synergyViolation]
                assert archHeuristics != null;
                double childInterference = archHeuristics.get(2);

                if (childInterference < (double) parent.getAttribute("InterInstrViolation")) {
                    feasibleMove = true; // instrument move is feasible

                    possibleUniqueInstrumentMoves.remove(moveChoiceIndex);
                    numberOfMoves += 1;
                } else  {
                    // Replace moved instrument and try again
                    removalOrbitPayload.add(instrumentToMove);
                    payloads.set(moveChoice.get(1), removalOrbitPayload);

                    if (moveInstruments) {
                        additionOrbitPayload = payloads.get(moveChoice.get(2));
                        additionOrbitPayload.remove(instrumentToMove);
                        payloads.set(moveChoice.get(2), additionOrbitPayload);
                    }

                    possibleUniqueInstrumentMoves.remove(moveChoiceIndex);
                    numberOfTries += 1;
                }
            }
            numberOfMoves += 1;
        }
        this.resourcePool.freeResource(res);
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
                                int presentValueInteger = presentValueIntegers.get(k);

                                for (int m = 0; m < possibleSatelliteMoves.size(); m++) {
                                    int possibleMoveSatellite = possibleSatelliteMoves.get(m);
                                    ValueVector otherSatelliteInstruments = allSatellites.get(possibleMoveSatellite).getSlotValue("instruments").listValue(r.getGlobalContext());
                                    ArrayList<String> otherSatelliteInstrumentsArray = getStringArrayFromValueVector(r, otherSatelliteInstruments);
                                    String instrument = satelliteInstruments.get(presentValueInteger).stringValue(r.getGlobalContext());
                                    if (otherSatelliteInstrumentsArray.contains(instrument)) { // Check if current instrument is present in the other candidate satellite
                                        continue;
                                    } else {
                                        ArrayList<Integer> removalChoice = new ArrayList<>();
                                        removalChoice.add(presentValueInteger);
                                        removalChoice.add(i);
                                        removalChoice.add(possibleMoveSatellite);
                                        if (!validMoves.contains(removalChoice)) {
                                            validMoves.addAll(Collections.singleton(removalChoice));
                                        }
                                    }
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

    private ArrayList<String> getStringArrayFromValueVector (Rete rete, ValueVector vv) throws JessException {
        ArrayList<String> stringArray = new ArrayList<>();
        for (int i = 0; i < vv.size(); i++) {
            stringArray.add(vv.get(i).stringValue(rete.getGlobalContext()));
        }
        return  stringArray;
    }

    @Override
    public int getArity() {
        return 1;
    }


}
