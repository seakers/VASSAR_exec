package seakers.vassarexecheur.search.operators.assigning;

import jess.Fact;
import jess.JessException;
import jess.ValueVector;
import org.moeaframework.core.PRNG;
import org.moeaframework.core.Solution;
import org.moeaframework.core.Variation;
import jess.Rete;
import org.moeaframework.core.variable.BinaryVariable;
import org.moeaframework.core.variable.RealVariable;
import seakers.vassarexecheur.search.problems.assigning.AssigningProblem;
import seakers.vassarheur.BaseParams;
import seakers.vassarheur.QueryBuilder;
import seakers.vassarexecheur.search.problems.assigning.AssigningArchitecture;
import seakers.vassarheur.Resource;
import seakers.vassarheur.ResourcePool;
import seakers.vassarheur.architecture.AbstractArchitecture;
import seakers.vassarheur.problems.Assigning.ArchitectureEvaluator;
import seakers.vassarheur.problems.Assigning.ClimateCentricAssigningParams;
import seakers.vassarheur.utils.MatlabFunctions;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Checks to see if the spacecraft wet mass is less than the threshold.
 * If not, one or more instruments are moved from one or more satellites to alleviate the issue.
 *
 * @author roshansuresh
 */

public class RepairMassAssigning implements Variation {

    /**
     * The wet mass that a spacecraft must be at or higher
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
     * Toggles movement of removed instruments into other spacecraft
     */
    private boolean moveInstruments;

    /**
     * Assigning Problem used to evaluate architectures
     */
    private AssigningProblem problem;

    /**
     * Resource pool to obtain matlab functions
     */
    private final ResourcePool resourcePool;

    private final ArchitectureEvaluator evaluator;

    public RepairMassAssigning(double threshold, int xInstruments, BaseParams params, AssigningProblem problem, ResourcePool resourcePool, ArchitectureEvaluator evaluator) {
        this.xInstruments = xInstruments;
        this.threshold = threshold;
        this.params = params;
        this.moveInstruments = true;
        this.problem = problem;
        this.resourcePool = resourcePool;
        this.evaluator = evaluator;
        //this.pprng = new ParallelPRNG();
    }

    public RepairMassAssigning(double threshold, int xInstruments, BaseParams params, boolean moveInstruments, AssigningProblem problem, ResourcePool resourcePool, ArchitectureEvaluator evaluator) {
        this.xInstruments = xInstruments;
        this.threshold = threshold;
        this.params = params;
        this.moveInstruments = moveInstruments;
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
        AssigningArchitecture parent = (AssigningArchitecture) solutions[0];
        ArrayList<ArrayList<String>> payloads = parent.getSatellitePayloads();
        ArrayList<String> orbits = parent.getSatelliteOrbits();

        ArrayList<ArrayList<String>> childPayloads = new ArrayList<>(payloads);
        ArrayList<String> childOrbits = new ArrayList<>(orbits);

        ArrayList<ArrayList<Double>> operatorParameters = parent.getOperatorParameters(); //{duty cycle, wet mass, packing efficiency} for each satellite
        ArrayList<Integer> candidateSatellites = new ArrayList<>();

        Resource res = resourcePool.getResource();
        MatlabFunctions m = res.getM();
        Rete rete = res.getRete();
        QueryBuilder queryBuilder = res.getQueryBuilder();

        for (int i = 0; i < payloads.size(); i++) {
            if ((operatorParameters.get(i).get(1)) > threshold) {
                candidateSatellites.add(i);
            }
        }

        if ((candidateSatellites.size() > 0) && (candidateSatellites.size() >= xInstruments)) {
            int moveCount = 0;
            int satelliteIndex;
            while ((moveCount < xInstruments) && (candidateSatellites.size() > 0)) {
                if (moveInstruments) {
                    int maxFeasibleTries = 3;
                    boolean feasibleMove = false;
                    int numberOfFeasibleTries = 0;

                    satelliteIndex = PRNG.nextInt(candidateSatellites.size());
                    int candidateSatelliteIndex = candidateSatellites.get(satelliteIndex);
                    //ArrayList<String> currentPayloads = new ArrayList<>();
                    if (childPayloads.get(candidateSatelliteIndex).isEmpty()) {
                        //childOrbits.remove(satelliteIndex);
                        //childPayloads.remove(satelliteIndex);
                        moveCount--;
                        continue;
                    }
                    else {
                        // Remove payload from current satellite
                        int payloadIndex = PRNG.nextInt(childPayloads.get(candidateSatelliteIndex).size());
                        String payload = childPayloads.get(candidateSatelliteIndex).get(payloadIndex);
                        ArrayList<Integer> candidatePayloadSatellites = getCandidateSatellitesForPayload(childPayloads, payload);
                        if (candidatePayloadSatellites.size() == 0) {
                            moveCount--;
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

                        while ((!feasibleMove) && (numberOfFeasibleTries < maxFeasibleTries)) {

                            //ArrayList<String> candidateSatellitePayload = childPayloads.get(candidatePayloadSatelliteIndex);
                            //candidateSatellitePayload.add(payload);
                            childPayloads.get(candidatePayloadSatellite).add(payload);
                            //childPayloads.set(satelliteIndex, currentPayloads);

                            AssigningArchitecture candidateChild = getArchitectureFromPayloadsAndOrbits(childPayloads, childOrbits);
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
                            double childMass = archHeuristics.get(4);

                            if (childMass < (double) parent.getAttribute("SpMassViolation")) {
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
                                        break;
                                    }

                                    int newPayloadIndex = PRNG.nextInt(candidateSatellitePayloads.size());
                                    payload = candidateSatellitePayloads.get(newPayloadIndex);
                                    candidatePayloadSatellites = getCandidateSatellitesForPayload(childPayloads, payload);
                                    if (candidatePayloadSatellites.size() == 0) {
                                        moveCount--;
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
                } else {
                    satelliteIndex = PRNG.nextInt(candidateSatellites.size());
                    int candidateSatelliteIndex = candidateSatellites.get(satelliteIndex);
                    //ArrayList<String> currentPayloads = new ArrayList<>();
                    if (childPayloads.get(candidateSatelliteIndex).isEmpty()) {
                        //childOrbits.remove(satelliteIndex);
                        //childPayloads.remove(satelliteIndex);
                        moveCount--;
                        continue;
                    }

                    // Remove payload from current satellite
                    int payloadIndex = PRNG.nextInt(childPayloads.get(candidateSatelliteIndex).size());
                    String payload = childPayloads.get(candidateSatelliteIndex).get(payloadIndex);
                    ArrayList<Integer> candidatePayloadSatellites = getCandidateSatellitesForPayload(childPayloads, payload);
                    if (candidatePayloadSatellites.size() == 0) {
                        moveCount--;
                        continue;
                    }
                    //for (int m = 0; m < childPayloads.get(satelliteIndex).size(); m++) {
                    //if (m != payloadIndex) {
                    //currentPayloads.add(childPayloads.get(satelliteIndex).get(m));
                    //}
                    //}
                    childPayloads.get(candidateSatelliteIndex).remove(payload);
                }
                candidateSatellites.remove(satelliteIndex);
                moveCount++;
            }
        }
        this.resourcePool.freeResource(res);
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
