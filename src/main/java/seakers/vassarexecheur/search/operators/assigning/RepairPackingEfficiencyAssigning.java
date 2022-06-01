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
import seakers.vassarheur.problems.Assigning.ArchitectureEvaluator;
import seakers.vassarheur.problems.Assigning.ClimateCentricAssigningParams;
import seakers.vassarheur.utils.MatlabFunctions;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

/**
 * Checks to see if packing efficiency of the architecture falls within the acceptable bounds.
 * If not, one or more instruments are moved from one or more orbits to try to alleviate the situation
 *
 * @author roshansuresh
 */

public class RepairPackingEfficiencyAssigning implements Variation {
    /**
     * The packing efficiency that a spacecraft must be at or higher
     */
    private final double threshold;

    /**
     * The number of instruments to remove from each satellite that does not
     * meet the threshold
     */
    private final int xInstruments; // NOT USED IN CURRENT VERSION

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

    public RepairPackingEfficiencyAssigning(double threshold, int xInstruments, BaseParams params, AssigningProblem problem, ResourcePool resourcePool, ArchitectureEvaluator evaluator) {
        this.xInstruments = xInstruments;
        this.threshold = threshold;
        this.params = params;
        this.moveInstruments = true;
        this.problem = problem;
        this.resourcePool = resourcePool;
        this.evaluator = evaluator;
        //this.pprng = new ParallelPRNG();
    }

    public RepairPackingEfficiencyAssigning(double threshold, int xInstruments, BaseParams params, boolean moveInstruments, AssigningProblem problem, ResourcePool resourcePool, ArchitectureEvaluator evaluator) {
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
        //System.out.println("parent");
        //System.out.println(parent.getBitString());

        AbstractArchitecture arch_abs = problem.getAbstractArchitecture(parent);

        Resource res = resourcePool.getResource();
        MatlabFunctions m = res.getM();
        Rete rete = res.getRete();
        QueryBuilder queryBuilder = res.getQueryBuilder();

        evaluator.assertMissions(params, rete, arch_abs, m);

        try {
            evaluator.evaluateHeuristicParameters(rete, arch_abs, queryBuilder, m);
        } catch (JessException e) {
            e.printStackTrace();
        }

        ArrayList<Fact> satellites = queryBuilder.makeQuery("MANIFEST::Satellite");
        ArrayList<Fact> instrumentFacts = queryBuilder.makeQuery("DATABASE::Instrument");

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

        ArrayList<ArrayList<String>> childPayloads = new ArrayList<>(payloads);
        ArrayList<String> childOrbits = new ArrayList<>(orbits);

        ArrayList<ArrayList<Double>> operatorParameters = parent.getOperatorParameters(); //{duty cycle, wet mass, packing efficiency} for each satellite

        // Store satellites which violate heuristic
        //for (int i = 0; i < payloads.size(); i++) {
            //if ((operatorParameters.get(i).get(2)) < threshold) {
                //int satelliteFactIndex = 0;
                //for (int j = 0; j < satellites.size(); j++) {
                    //String satelliteOrbit = null;
                    //try {
                        //satelliteOrbit = satellites.get(j).getSlotValue("orbit-string").stringValue(rete.getGlobalContext());
                    //} catch (JessException e) {
                        //e.printStackTrace();
                    //}
                    //if (satelliteOrbit.equalsIgnoreCase(orbits.get(i))) {
                        //satelliteFactIndex = j;
                        //break;
                    //}
                //}
                //candidateSatellites.add(satelliteFactIndex);
                //satellitesFactToList.put(satelliteFactIndex, i);
            //}
        //}

        // Determine worst heuristic violating satellite
        int worstSatelliteListIndex = getWorstSatellite(operatorParameters, payloads);
        //System.out.println("Worst Satellite: " + worstSatelliteListIndex);
        if (worstSatelliteListIndex != -1) {
            int worstSatelliteFactIndex = 0;
            for (int j = 0; j < satellites.size(); j++) {
                String satelliteOrbit = null;
                try {
                    satelliteOrbit = getOrbitStringFromSatellite(satellites.get(j), rete);
                } catch (JessException e) {
                    e.printStackTrace();
                }
                if (satelliteOrbit.equalsIgnoreCase(orbits.get(worstSatelliteListIndex))) {
                    worstSatelliteFactIndex = j;
                    break;
                }
            }

            int satelliteIndex;
            boolean noChange = false;
            if (moveInstruments) {
                //satelliteIndex = PRNG.nextInt(candidateSatellites.size());
                //int candidateSatelliteIndex = candidateSatellites.get(satelliteIndex);
                //int candidateSatelliteListIndex = satellitesFactToList.get(candidateSatelliteIndex);
                //if (childPayloads.get(candidateSatelliteListIndex).isEmpty()) {
                //while (!noChange) {
                //candidateSatellites.remove(candidateSatelliteIndex); // Remove current satellite from candidates list
                //if (candidateSatellites.size() == 0) {
                //noChange = true;
                //}
                //if (!noChange) {
                //satelliteIndex = PRNG.nextInt(candidateSatellites.size()); // find a different satellite
                //candidateSatelliteIndex = candidateSatellites.get(satelliteIndex);
                //candidateSatelliteListIndex = satellitesFactToList.get(candidateSatelliteIndex);
                //if (childPayloads.get(candidateSatelliteListIndex).size() > 0) {
                //break;
                //}
                //}
                //}
                //}

                try {
                    ValueVector satelliteInstruments = satellites.get(worstSatelliteFactIndex).getSlotValue("instruments").listValue(rete.getGlobalContext());
                    for (int i = 0; i < satelliteInstruments.size(); i++) {
                        String instrument = satelliteInstruments.get(i).stringValue(rete.getGlobalContext());
                        boolean noFeasibleSatelliteToMoveTo = true;
                        for (int j = 0; j < satellites.size(); j++) {
                            if (j == worstSatelliteFactIndex) {
                                continue;
                            }
                            double satelliteLaunchMass = satellites.get(j).getSlotValue("satellite-launch-mass").floatValue(rete.getGlobalContext());
                            Fact instrumentFact = getInstrumentFact(instrumentFacts, instrument, rete);
                            double instrumentMass = instrumentFact.getSlotValue("mass#").floatValue(rete.getGlobalContext());
                            String launchVehicle = satellites.get(j).getSlotValue("launch-vehicle").stringValue(rete.getGlobalContext());
                            double orbitAltitude = satellites.get(j).getSlotValue("orbit-altitude#").floatValue(rete.getGlobalContext());
                            String orbitType = satellites.get(j).getSlotValue("orbit-type").stringValue(rete.getGlobalContext());
                            String orbitInclination = satellites.get(j).getSlotValue("orbit-inclination").stringValue(rete.getGlobalContext());
                            double launchVehicleCapacity = getLaunchVehicleCapacity(res, rete, launchVehicle, orbitType, orbitInclination, orbitAltitude);
                            int numberOfLaunches = satellites.get(j).getSlotValue("num-launches").intValue(rete.getGlobalContext());
                            int moveSatelliteListIndex = getCorrespondingListIndex(j, satellites, orbits, rete);
                            if ((satelliteLaunchMass + (3*instrumentMass) < (numberOfLaunches*launchVehicleCapacity)) && (!payloads.get(moveSatelliteListIndex).contains(instrument))) {
                                //System.out.println("Feasibly moving instrument " + instrument + " from satellite " + worstSatelliteListIndex + " to satellite " + moveSatelliteListIndex);
                                childPayloads.get(worstSatelliteListIndex).remove(instrument);
                                childPayloads.get(moveSatelliteListIndex).add(instrument);
                                noFeasibleSatelliteToMoveTo = false;
                                break;
                            }
                        }
                        if (noFeasibleSatelliteToMoveTo) {
                            ArrayList<Integer> candidatePayloadSatellites = getCandidateSatellitesForPayload(childPayloads, instrument);
                            if (candidatePayloadSatellites.size() != 0) {
                                int candidatePayloadSatelliteIndex = PRNG.nextInt(candidatePayloadSatellites.size());
                                int candidatePayloadSatellite = candidatePayloadSatellites.get(candidatePayloadSatelliteIndex);
                                //System.out.println("Infeasibly moving instrument " + instrument + " from satellite " + worstSatelliteListIndex + " to satellite " + candidatePayloadSatellite);
                                payloads.get(worstSatelliteListIndex).remove(instrument);
                                payloads.get(candidatePayloadSatellite).add(instrument);
                            }
                        }
                    }
                } catch (JessException e) {
                    e.printStackTrace();
                }
            } else {
                ValueVector satelliteInstruments = null;
                try {
                    satelliteInstruments = satellites.get(worstSatelliteFactIndex).getSlotValue("instruments").listValue(rete.getGlobalContext());
                    for (int i = 0; i < satelliteInstruments.size(); i++) {
                        String instrument = satelliteInstruments.get(i).stringValue(rete.getGlobalContext());
                        childPayloads.get(worstSatelliteListIndex).remove(instrument);
                    }
                } catch (JessException e) {
                    e.printStackTrace();
                }
            }
        }

        this.resourcePool.freeResource(res);
        AssigningArchitecture child = getArchitectureFromPayloadsAndOrbits(childPayloads, childOrbits);

        return new Solution[]{child};
    }

    private String getOrbitStringFromSatellite(Fact satellite, Rete r) throws JessException {
        String satelliteOrbitType = satellite.getSlotValue("orbit-type").stringValue(r.getGlobalContext());
        int satelliteOrbitAltitude = (int) satellite.getSlotValue("orbit-altitude#").floatValue(r.getGlobalContext());
        String satelliteOrbitInclination = satellite.getSlotValue("orbit-inclination").stringValue(r.getGlobalContext());
        String satelliteOrbitRAAN = satellite.getSlotValue("orbit-RAAN").stringValue(r.getGlobalContext());
        return satelliteOrbitType + "-" + satelliteOrbitAltitude + "-" + satelliteOrbitInclination + "-" + satelliteOrbitRAAN;
    }

    private int getWorstSatellite(ArrayList<ArrayList<Double>> operatorParams, ArrayList<ArrayList<String>> payloadStrings) {
        int worstSatelliteIndex = -1;
        double worstPackingEfficiency = 1;
        for (int i = 0; i < operatorParams.size(); i++) {
            if ((operatorParams.get(i).get(2) < worstPackingEfficiency) && (payloadStrings.get(i).size() != 0) && (operatorParams.get(i).get(2) < threshold)) {
                worstPackingEfficiency = operatorParams.get(i).get(2);
                worstSatelliteIndex = i;
            }
        }
        return worstSatelliteIndex;
    }

    private int getCorrespondingListIndex(int satelliteFactIndex, ArrayList<Fact> satelliteFacts, ArrayList<String> orbits, Rete r) throws JessException {
        int satelliteListIndex = 0;
        String satelliteOrbit = satelliteFacts.get(satelliteFactIndex).getSlotValue("orbit-string").stringValue(r.getGlobalContext());
        for (int i = 0; i < orbits.size(); i++) {
            if (orbits.get(i).equalsIgnoreCase(satelliteOrbit)) {
                satelliteListIndex = i;
                break;
            }
        }
        return satelliteListIndex;
    }

    /**
     * Extracts the Instrument fact corresponding to the instrument name as part of the payload in a satellite
     * @param instrumentFacts
     * @param instrumentName
     * @return
     */
    private Fact getInstrumentFact(ArrayList<Fact> instrumentFacts, String instrumentName, Rete r) throws JessException {
        Fact instrumentFact = null;
        for (int i = 0; i < instrumentFacts.size(); i++) {
            Fact currentInstrumentFact = instrumentFacts.get(i);
            if (instrumentName.equalsIgnoreCase(currentInstrumentFact.getSlotValue("Name").stringValue(r.getGlobalContext()))) {
                instrumentFact = currentInstrumentFact;
                break;
            }
        }
        return instrumentFact;
    }

    private double getLaunchVehicleCapacity(Resource resource, Rete r, String launchVehicle, String orbitType, String orbitInclination, double orbitAltitude) throws JessException {
        //ValueVector vv = new ValueVector();
        String payloadCoeffsKey = orbitType + "-" + orbitInclination;
        //vv.add(launchVehicle);
        //vv.add(payloadCoeffsKey);
        ArrayList<String> inputs = new ArrayList<>();
        inputs.add(launchVehicle);
        inputs.add(payloadCoeffsKey);
        //Value launchVehicleCoefficients = resource.getM().getLaunchVehiclePerformanceCoeffs((Funcall) vv, r.getGlobalContext());
        Value launchVehicleCoefficients = resource.getM().getLaunchVehiclePerformanceCoeffs2(inputs, r.getGlobalContext());
        return launchVehicleCoefficients.listValue(r.getGlobalContext()).get(0).floatValue(r.getGlobalContext())*1 +
                launchVehicleCoefficients.listValue(r.getGlobalContext()).get(1).floatValue(r.getGlobalContext())*orbitAltitude +
                launchVehicleCoefficients.listValue(r.getGlobalContext()).get(2).floatValue(r.getGlobalContext())*Math.pow(orbitAltitude, 2);
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
