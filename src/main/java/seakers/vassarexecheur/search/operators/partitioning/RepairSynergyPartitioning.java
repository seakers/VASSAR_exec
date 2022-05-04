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
import seakers.vassarexecheur.search.problems.partitioning.PartitioningProblem;
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
import java.util.stream.Collectors;

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

    /**
     * Assigning Problem used to evaluate architectures
     */
    private PartitioningProblem problem;

    private final HashMap<String, String[]> synergyMap;

    public RepairSynergyPartitioning(int numChanges, ResourcePool resourcePool, ArchitectureEvaluator evaluator, BaseParams params, PartitioningProblem problem, HashMap<String, String[]> synergyMap) {
        this.numberOfChanges = numChanges;
        this.resourcePool = resourcePool;
        this.evaluator = evaluator;
        this.params = params;
        this.problem = problem;
        this.synergyMap = synergyMap;
        //this.pprng = new ParallelPRNG();
    }

    @Override
    public int getArity() {
        return 1;
    }

    @Override
    public Solution[] evolve(Solution[] solutions) {
        PartitioningArchitecture parent = (PartitioningArchitecture) solutions[0];

        AbstractArchitecture arch_abs = problem.getAbstractArchitecture(parent);

        Resource res = resourcePool.getResource();
        MatlabFunctions m = res.getM();
        Rete rete = res.getRete();
        QueryBuilder queryBuilder = res.getQueryBuilder();

        //HashMap<String, String[]> instrumentSynergyMap = getInstrumentSynergyNameMap(params);

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

        ArrayList<ArrayList<Integer>> possibleInstrumentMoves = getValidMoves(rete, satellites, synergyMap);
        ArrayList<ArrayList<Integer>> possibleUniqueInstrumentMoves = (ArrayList<ArrayList<Integer>>) possibleInstrumentMoves.stream().distinct().collect(Collectors.toList());

        // Make choices of instrument move randomly
        int numberOfMoves = 0;
        while ((numberOfMoves < numberOfChanges) && (possibleUniqueInstrumentMoves.size() > 0)) {
            int numberOfTries = 0;
            boolean feasibleMove = false;
            while ((!feasibleMove) && (numberOfTries < 5) && (possibleUniqueInstrumentMoves.size() > 0)) {
                int moveChoiceIndex = PRNG.nextInt(possibleUniqueInstrumentMoves.size());
                ArrayList<Integer> moveChoice = possibleUniqueInstrumentMoves.get(moveChoiceIndex);

                // Update payloads
                ArrayList<String> removalOrbitPayload = payloads.get(moveChoice.get(1));
                String instrumentToMove = removalOrbitPayload.get(moveChoice.get(0));
                removalOrbitPayload.remove(instrumentToMove);

                payloads.set(moveChoice.get(1), removalOrbitPayload);

                ArrayList<String> additionOrbitPayload = payloads.get(moveChoice.get(2));
                additionOrbitPayload.add(instrumentToMove);

                payloads.set(moveChoice.get(2), additionOrbitPayload);

                // Check if instrument movement improves heuristic satisfaction
                PartitioningArchitecture candidateChild = getArchitectureFromPayloadsAndOrbits(payloads, orbits);
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
                double childSynergy = archHeuristics.get(5);

                if (childSynergy < (double) parent.getAttribute("SynergyViolation")) {
                    feasibleMove = true; // instrument move is feasible

                    possibleUniqueInstrumentMoves.remove(moveChoiceIndex);
                    numberOfMoves += 1;
                } else  {
                    // Replace moved instrument and try again
                    removalOrbitPayload.add(instrumentToMove);
                    payloads.set(moveChoice.get(1), removalOrbitPayload);

                    additionOrbitPayload.remove(instrumentToMove);
                    payloads.set(moveChoice.get(2), additionOrbitPayload);

                    possibleUniqueInstrumentMoves.remove(moveChoiceIndex);
                    numberOfTries += 1;
                }
            }
            numberOfMoves += 1;
        }
        this.resourcePool.freeResource(res);
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
                            ArrayList<Integer> presentValueIntegers = checkInstrumentsinSatellite(r, otherSatelliteInstruments, synergyInstruments);
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
        ArrayList<String> instrumentList = new ArrayList<>(Arrays.asList(params.getInstrumentList()));
        ArrayList<String> orbitList = new ArrayList<>(Arrays.asList(params.getOrbitList()));

        PartitioningArchitecture arch = new PartitioningArchitecture(instrumentList.size(), orbitList.size(), 2, params);

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
