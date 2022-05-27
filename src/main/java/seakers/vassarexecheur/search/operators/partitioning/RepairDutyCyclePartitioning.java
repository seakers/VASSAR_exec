package seakers.vassarexecheur.search.operators.partitioning;

import jess.Fact;
import jess.JessException;
import jess.ValueVector;
import org.moeaframework.core.PRNG;
import org.moeaframework.core.Solution;
import org.moeaframework.core.Variation;
import jess.Rete;
import seakers.architecture.util.IntegerVariable;
import seakers.vassarexecheur.search.problems.partitioning.PartitioningArchitecture;
import seakers.vassarexecheur.search.problems.partitioning.PartitioningProblem;
import seakers.vassarheur.BaseParams;
import seakers.vassarheur.QueryBuilder;
import seakers.vassarheur.Resource;
import seakers.vassarheur.ResourcePool;
import seakers.vassarheur.architecture.AbstractArchitecture;
import seakers.vassarheur.problems.PartitioningAndAssigning.ArchitectureEvaluator;
import seakers.vassarheur.utils.MatlabFunctions;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Checks that the data rate duty cycle of the spacecraft fall within the
 * acceptable bounds. If not one or more random instruments are moved to another
 * satellite to try to alleviate the situation. User can define whether to change
 * one or multiple spacecraft. Class modified from the original RepairDutyCycle.java
 * in the EOSS repo authored by nozomihitomi
 *
 * @author roshansuresh
 */

public class RepairDutyCyclePartitioning implements Variation {

    /**
     * The duty cycle that a spacecraft must be at or higher
     */
    private final double threshold;

    /**
     * The number of instruments to remove from each satellite that does not
     * meet the threshold
     */
    private final int xInstruments;

    //private final ParallelPRNG pprng;

    private final BaseParams params;

    /**
     * Partitioning Problem used to evaluate architectures
     */
    private PartitioningProblem problem;

    /**
     * Resource pool to obtain matlab functions
     */
    private final ResourcePool resourcePool;

    private final ArchitectureEvaluator evaluator;

    public RepairDutyCyclePartitioning(double threshold, int xInstruments, BaseParams params, PartitioningProblem problem, ResourcePool resourcePool, ArchitectureEvaluator evaluator) {
        this.xInstruments = xInstruments;
        this.threshold = threshold;
        this.params = params;
        this.problem = problem;
        this.resourcePool = resourcePool;
        this.evaluator = evaluator;
        //this.pprng = new ParallelPRNG();
    }

    @Override
    public int getArity() {
        return 1;
    }

    @Override
    public Solution[] evolve(Solution[] solutions) {
        PartitioningArchitecture parent = (PartitioningArchitecture) solutions[0];
        if (parent.getSatellitePayloads() != null) {
            AbstractArchitecture arch_abs = problem.getAbstractArchitecture(parent);

            Resource res = resourcePool.getResource();
            MatlabFunctions m = res.getM();
            Rete rete = res.getRete();
            QueryBuilder queryBuilder = res.getQueryBuilder();

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

            //ArrayList<ArrayList<String>> payloads = parent.getSatellitePayloads();
            //ArrayList<String> orbits = parent.getSatelliteOrbits();

            ArrayList<String> allOrbits = addEmptyOrbits(orbits, params.getOrbitList());
            ArrayList<ArrayList<String>> allPayloads = addEmptyPayloads(payloads, orbits, params.getOrbitList());

            ArrayList<ArrayList<String>> childPayloads = new ArrayList<>(allPayloads);
            ArrayList<String> childOrbits = new ArrayList<>(allOrbits);

            ArrayList<ArrayList<Double>> operatorParameters = parent.getOperatorParameters(); //{duty cycle, wet mass, packing efficiency} for each satellite
            ArrayList<ArrayList<Double>> allOperatorParameters = orderOperatorParameters(operatorParameters, orbits, params.getOrbitList());
            ArrayList<Integer> candidateSatellites = new ArrayList<>();

            for (int i = 0; i < allPayloads.size(); i++) {
                if ((allOperatorParameters.get(i).get(0)) < threshold) {
                    candidateSatellites.add(i);
                }
            }

            if ((candidateSatellites.size() > 0) && (candidateSatellites.size() >= xInstruments)) {
                int moveCount = 0;
                int satelliteIndex;
                while ((moveCount < xInstruments) && (candidateSatellites.size() > 0)) {
                    int maxFeasibleTries = 3;
                    boolean feasibleMove = false;
                    int numberOfFeasibleTries = 0;

                    satelliteIndex = PRNG.nextInt(candidateSatellites.size());
                    int candidateSatelliteIndex = candidateSatellites.get(satelliteIndex);
                    if (childPayloads.get(candidateSatelliteIndex).isEmpty()) {
                        candidateSatellites.remove(satelliteIndex);
                        //System.out.println("Duty Cycle continue 1");
                        continue;
                    }
                    else {
                        // Remove payload from current satellite
                        int payloadIndex = PRNG.nextInt(childPayloads.get(candidateSatelliteIndex).size());
                        String payload = childPayloads.get(candidateSatelliteIndex).get(payloadIndex);
                        ArrayList<Integer> candidatePayloadSatellites = getCandidateSatellitesForPayload(childPayloads, payload);
                        if (candidatePayloadSatellites.size() == 0) {
                            //System.out.println("Duty Cycle continue 2");
                            continue;
                        }
                        //for (int m = 0; m < childPayloads.get(satelliteIndex).size(); m++) {
                        //if (m != payloadIndex) {
                        //currentPayloads.add(childPayloads.get(satelliteIndex).get(m));
                        //}
                        //}
                        childPayloads.get(candidateSatelliteIndex).remove(payload);
                        String originalPayload = payload;

                        // Add payload to different satellite
                        int candidatePayloadSatelliteIndex = PRNG.nextInt(candidatePayloadSatellites.size());
                        int candidatePayloadSatellite = candidatePayloadSatellites.get(candidatePayloadSatelliteIndex);

                        boolean breakNest = false;

                        while ((!feasibleMove) && (numberOfFeasibleTries < maxFeasibleTries)) {
                            if (breakNest) {
                                moveCount--;
                                break;
                            }

                            //ArrayList<String> candidateSatellitePayload = childPayloads.get(candidatePayloadSatelliteIndex);
                            //candidateSatellitePayload.add(payload);
                            childPayloads.get(candidatePayloadSatellite).add(payload);
                            //childPayloads.set(satelliteIndex, currentPayloads);

                            PartitioningArchitecture candidateChild = getArchitectureFromPayloadsAndOrbits(childPayloads, childOrbits);
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
                            double childDutyCycle = archHeuristics.get(0);

                            if (childDutyCycle < (double) parent.getAttribute("DCViolation")) {
                                feasibleMove = true;
                            } else {
                                // Remove the added instrument to try again
                                childPayloads.get(candidatePayloadSatellite).remove(payload);
                                candidatePayloadSatellites.remove(candidatePayloadSatelliteIndex);

                                if (candidatePayloadSatellites.size() == 0) { // If the candidate instrument cannot be added to other satellites, choose a different instrument from the same satellite and start again
                                    childPayloads.get(candidateSatelliteIndex).add(originalPayload);

                                    ArrayList<String> candidateSatellitePayloads = new ArrayList<>();
                                    candidateSatellitePayloads.addAll(childPayloads.get(candidateSatelliteIndex));
                                    candidateSatellitePayloads.remove(payload);

                                    if (candidateSatellitePayloads.size() == 0) { // Try again with a different instrument after replacing original instrument in its satellite
                                        // move counter is not reset to prevent infinite looping over all instruments in all satellites  (limits to moveCounter)
                                        breakNest = true;
                                        continue;
                                    }

                                    int newPayloadIndex = PRNG.nextInt(candidateSatellitePayloads.size());
                                    payload = candidateSatellitePayloads.get(newPayloadIndex);
                                    candidatePayloadSatellites = getCandidateSatellitesForPayload(childPayloads, payload);
                                    if (candidatePayloadSatellites.size() == 0) {
                                        breakNest = true;
                                        continue;
                                    }
                                    originalPayload = payload;
                                    childPayloads.get(candidateSatelliteIndex).remove(payload);
                                    //candidatePayloadSatelliteIndex = PRNG.nextInt(candidatePayloadSatellites.size());
                                    //candidatePayloadSatellite = candidatePayloadSatellites.get(candidatePayloadSatelliteIndex);
                                }

                                candidatePayloadSatelliteIndex = PRNG.nextInt(candidatePayloadSatellites.size());
                                candidatePayloadSatellite = candidatePayloadSatellites.get(candidatePayloadSatelliteIndex);

                                numberOfFeasibleTries++;
                                if (numberOfFeasibleTries == maxFeasibleTries) {
                                    childPayloads.get(candidateSatelliteIndex).add(originalPayload);
                                }
                            }
                        }
                    }
                    candidateSatellites.remove(satelliteIndex);
                    moveCount++;
                }
            }
            this.resourcePool.freeResource(res);
            PartitioningArchitecture child = getArchitectureFromPayloadsAndOrbits(childPayloads, childOrbits);

            return new Solution[]{child};
        } else {
            return new Solution[]{parent};
        }
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
        ArrayList<String> instrumentList = new ArrayList<>(Arrays.asList(params.getInstrumentList()));
        ArrayList<String> orbitList = new ArrayList<>(Arrays.asList(params.getOrbitList()));

        PartitioningArchitecture arch = new PartitioningArchitecture(instrumentList.size(), orbitList.size(), 2, params); // ADD HEURISTIC OBJECTIVES AND CONSTRAINTS
        int partitionIndex = 0;

        for (int i = 0; i < currentOrbits.size(); i++) {
            ArrayList<String> currentOrbitPayloads = currentPayloads.get(i);

            for (int j = 0; j < currentOrbitPayloads.size(); j++) {
                int payloadIndex = instrumentList.indexOf(currentOrbitPayloads.get(j));

                IntegerVariable instrVar = new IntegerVariable(partitionIndex, 0, instrumentList.size()-1);
                arch.setVariable(payloadIndex, instrVar);
            }

            if (currentOrbitPayloads.size() != 0) {
                int orbitIndex = orbitList.indexOf(currentOrbits.get(i));

                IntegerVariable orbitVar = new IntegerVariable(orbitIndex, -1, orbitList.size()-1);
                arch.setVariable(instrumentList.size()+partitionIndex, orbitVar);

                partitionIndex += 1;
            }
        }

        //// NOT REQUIRED SINCE PARTITIONING ARCHITECTURE IS INITIALIZED WITH ALL ASSIGNING VARIABLES AS -1
        //IntegerVariable noOrbitVar = new IntegerVariable(-1, -1, orbitList.size()-1);
        //for (int k = 0; k < instrumentList.size()-currentOrbits.size(); k++) {
            //arch.setVariable(instrumentList.size()+currentOrbits.size()+k, noOrbitVar);
        //}

        return arch;
    }

    private ArrayList<String> addEmptyOrbits (ArrayList<String> archOrbits, String[] orbitsList) {
        ArrayList<String> allOrbits = new ArrayList<>();
        allOrbits.addAll(archOrbits);
        allOrbits.addAll(Arrays.asList(orbitsList));
        return allOrbits;
    }

    private ArrayList<ArrayList<String>> addEmptyPayloads (ArrayList<ArrayList<String>> archPayloads, ArrayList<String> archOrbits, String[] orbitsList) {
        ArrayList<ArrayList<String>> archAllPayloads = new ArrayList<>();
        for (int i = 0; i < archOrbits.size(); i++) {
            archAllPayloads.add(archPayloads.get(i));
        }
        for (int j = 0; j < orbitsList.length; j++) {
            archAllPayloads.add(new ArrayList<>());
        }
        return archAllPayloads;
    }

    private ArrayList<ArrayList<Double>> orderOperatorParameters (ArrayList<ArrayList<Double>> archOperatorParameters, ArrayList<String> archOrbits, String[] orbitsList) {
        ArrayList<ArrayList<Double>> allOperatorParameters = new ArrayList<>(archOperatorParameters);
        for (int i = 0; i < orbitsList.length; i++) {
            ArrayList<Double> emptyOrbitParameters = new ArrayList<>();
            emptyOrbitParameters.add(1.0); // duty cycle
            emptyOrbitParameters.add(0.0); // wet mass in kg
            emptyOrbitParameters.add(1.0); // packing efficiency
            allOperatorParameters.add(emptyOrbitParameters);
        }
        return allOperatorParameters;
    }
}
