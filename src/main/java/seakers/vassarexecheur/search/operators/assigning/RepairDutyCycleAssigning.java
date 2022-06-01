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
 * Checks that the data rate duty cycle of the spacecraft fall within the
 * acceptable bounds. If not one or more random instruments are moved to different
 * satellites to try to alleviate the situation. User can define whether to change
 * one or multiple spacecraft. Class modified from the original RepairDutyCycle.java in the EOSS
 * repo authored by nozomihitomi
 *
 * @author roshansuresh
 */

public class RepairDutyCycleAssigning implements Variation {

    /**
     * The duty cycle that a spacecraft must be at or higher
     */
    private final double threshold;

    /**
     * The number of instruments to move from each satellite that does not
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

    public RepairDutyCycleAssigning(double threshold, int xInstruments, BaseParams params, AssigningProblem problem, ResourcePool resourcePool, ArchitectureEvaluator evaluator) {
        this.xInstruments = xInstruments;
        this.threshold = threshold;
        this.params = params;
        this.moveInstruments = true;
        this.problem = problem;
        this.resourcePool = resourcePool;
        this.evaluator = evaluator;
        //this.pprng = new ParallelPRNG();
    }

    public RepairDutyCycleAssigning(double threshold, int xInstruments, BaseParams params, boolean moveInstruments, AssigningProblem problem, ResourcePool resourcePool, ArchitectureEvaluator evaluator) {
        this.xInstruments = xInstruments;
        this.threshold = threshold;
        this.params = params;
        this.moveInstruments = moveInstruments;
        this.problem = problem;
        this.resourcePool = resourcePool;
        this.evaluator = evaluator;
        //this.pprng = new ParallelPRNG();
    }

    /**
     * Moves x number of instruments from the payload of y number of satellite
     * that does not meet the data rate duty cycle threshold
     *
     * @param solutions
     * @return
     */
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
            //int continueCounter = 0;
            while ((moveCount < xInstruments) && (candidateSatellites.size() > 0)) {
                if (moveInstruments) {
                    //System.out.println("Move Count: " + moveCount);
                    int maxFeasibleTries = 3;
                    boolean feasibleMove = false;
                    int numberOfFeasibleTries = 0;

                    if (candidateSatellites.isEmpty()) {
                        break;
                    }

                    satelliteIndex = PRNG.nextInt(candidateSatellites.size());
                    int candidateSatelliteIndex = candidateSatellites.get(satelliteIndex);
                    //ArrayList<String> currentPayloads = new ArrayList<>();
                    if (childPayloads.get(candidateSatelliteIndex).isEmpty()) {
                        //childOrbits.remove(satelliteIndex);
                        //childPayloads.remove(satelliteIndex);
                        candidateSatellites.remove(satelliteIndex);
                        continue;
                    }
                    else {
                        // Remove payload from current satellite
                        int payloadIndex = PRNG.nextInt(childPayloads.get(candidateSatelliteIndex).size());
                        String payload = childPayloads.get(candidateSatelliteIndex).get(payloadIndex);
                        ArrayList<Integer> candidatePayloadSatellites = getCandidateSatellitesForPayload(childPayloads, payload);
                        if (candidatePayloadSatellites.size() == 0) {
                            //System.out.println("Candidate Payload cannot be moved to other satellites");
                            continue;
                        }
                        //for (int m = 0; m < childPayloads.get(satelliteIndex).size(); m++) {
                        //if (m != payloadIndex) {
                        //currentPayloads.add(childPayloads.get(satelliteIndex).get(m));
                        //}
                        //}
                        childPayloads.get(candidateSatelliteIndex).remove(payload);
                        String originalPayload = payload;
                        //System.out.println("Move choice");
                        //System.out.println("Removal satellite: " + candidateSatelliteIndex);
                        //System.out.println("Payload: " + payload);

                        // Add payload to different satellite
                        int candidatePayloadSatelliteIndex = PRNG.nextInt(candidatePayloadSatellites.size());
                        int candidatePayloadSatellite = candidatePayloadSatellites.get(candidatePayloadSatelliteIndex);

                        boolean breakNest = false;

                        while ((!feasibleMove) && (numberOfFeasibleTries < maxFeasibleTries)) {
                            //System.out.println("Number of feasible tries: " + numberOfFeasibleTries);
                            if (breakNest) {
                                moveCount--;
                                break;
                            }

                            //ArrayList<String> candidateSatellitePayload = childPayloads.get(candidatePayloadSatelliteIndex);
                            //candidateSatellitePayload.add(payload);
                            //System.out.println("Addition satellite: " + candidatePayloadSatellite);
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
                            double childDutyCycle = archHeuristics.get(0);

                            if (childDutyCycle < (double) parent.getAttribute("DCViolation")) {
                                feasibleMove = true;
                                //System.out.println("Feasible Move");
                            } else {
                                //System.out.println("Infeasible Move");
                                // Remove the added instrument to try again
                                childPayloads.get(candidatePayloadSatellite).remove(payload);
                                candidatePayloadSatellites.remove(candidatePayloadSatelliteIndex);

                                if (candidatePayloadSatellites.size() == 0) { // If the candidate instrument cannot be added to other satellites, choose a different instrument from the same satellite and start again
                                    //System.out.println("New Payload choice");
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
                                    //System.out.println("New Payload choice: " + payload);
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
                                //System.out.println("New Satellite choice: " + candidatePayloadSatellite);

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
                        continue;
                    }

                    // Remove payload from current satellite
                    int payloadIndex = PRNG.nextInt(childPayloads.get(candidateSatelliteIndex).size());
                    String payload = childPayloads.get(candidateSatelliteIndex).get(payloadIndex);
                    ArrayList<Integer> candidatePayloadSatellites = getCandidateSatellitesForPayload(childPayloads, payload);
                    if (candidatePayloadSatellites.size() == 0) {
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

    private ArrayList<String> addEmptyOrbits(ArrayList<String> archOrbits, String[] orbitsList) {
        ArrayList<String> allOrbits = new ArrayList<>();
        for (int i = 0; i < orbitsList.length; i++) {
            if (archOrbits.contains(orbitsList[i])) {
                int orbitIndex = archOrbits.indexOf(orbitsList[i]);
                allOrbits.add(archOrbits.get(orbitIndex));
            } else {
                allOrbits.add(orbitsList[i]);
            }
        }
        return allOrbits;
    }

    private ArrayList<ArrayList<String>> addEmptyPayloads(ArrayList<ArrayList<String>> archPayloads, ArrayList<String> archOrbits, String[] orbitsList) {
        ArrayList<ArrayList<String>> allPayloads = new ArrayList<>();
        for (int i = 0; i < orbitsList.length; i++) {
            if (archOrbits.contains(orbitsList[i])) {
                int orbitIndex = archOrbits.indexOf(orbitsList[i]);
                allPayloads.add(archPayloads.get(orbitIndex));
            } else {
                allPayloads.add(new ArrayList<>());
            }
        }
        return allPayloads;
    }

    private ArrayList<ArrayList<Double>> orderOperatorParameters(ArrayList<ArrayList<Double>> archOperatorParameters, ArrayList<String> archOrbits, String[] orbitsList) {
        ArrayList<ArrayList<Double>> allOperatorParameters = new ArrayList<>();
        for (int i = 0; i < orbitsList.length; i++) {
            if (archOrbits.contains(orbitsList[i])) {
                int orbitIndex = archOrbits.indexOf(orbitsList[i]);
                allOperatorParameters.add(archOperatorParameters.get(orbitIndex));
            } else {
                ArrayList<Double> emptyOrbitParameters = new ArrayList<>();
                emptyOrbitParameters.add(1.0); // duty cycle
                emptyOrbitParameters.add(0.0); // wet mass in kg
                emptyOrbitParameters.add(1.0); // packing efficiency
                allOperatorParameters.add(emptyOrbitParameters);
            }
        }
        return allOperatorParameters;
    }

    @Override
    public int getArity() {
        return 1;
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

}
