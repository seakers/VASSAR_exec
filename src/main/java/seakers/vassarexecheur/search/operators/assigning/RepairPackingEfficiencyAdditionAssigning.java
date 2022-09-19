package seakers.vassarexecheur.search.operators.assigning;

import jess.*;
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
import java.util.HashMap;

/**
 * Checks to see if packing efficiency of the architecture falls within the acceptable bounds.
 * If not, one or more instruments are added to one or more orbits to try to alleviate the situation
 *
 * @author roshansuresh
 */

public class RepairPackingEfficiencyAdditionAssigning implements Variation {

    /**
     * The packing efficiency that a spacecraft must be at or higher
     */
    private final double threshold;

    /**
     * The number of instruments to add to each satellite that does not
     * meet the threshold
     */
    private final int xInstruments;

    /**
     * The number of satellites to which to add instruments
     */
    private final int ySatellites;

    //private final ParallelPRNG pprng;

    private final BaseParams params;

    /**
     * Assigning Problem used to evaluate architectures
     */
    private AssigningProblem problem;

    /**
     * Resource pool to obtain matlab functions
     */
    private final ResourcePool resourcePool;

    private final ArchitectureEvaluator evaluator;

    public RepairPackingEfficiencyAdditionAssigning(double threshold, int xInstruments, int ySatellites, BaseParams params, AssigningProblem problem, ResourcePool resourcePool, ArchitectureEvaluator evaluator) {
        this.xInstruments = xInstruments;
        this.threshold = threshold;
        this.params = params;
        this.ySatellites = ySatellites;
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

        ArrayList<Integer> candidateSatellites = new ArrayList<>();
        HashMap<Integer, Integer> satellitesFactToList = new HashMap<>();

        // Store satellites which violate heuristic
        for (int i = 0; i < payloads.size(); i++) {
            if ((operatorParameters.get(i).get(2)) < threshold) {
                int satelliteFactIndex = 0;
                for (int j = 0; j < satellites.size(); j++) {
                    String satelliteOrbit = null;
                    try {
                        satelliteOrbit = satellites.get(j).getSlotValue("orbit-string").stringValue(rete.getGlobalContext());
                    } catch (JessException e) {
                        e.printStackTrace();
                    }
                    assert satelliteOrbit != null;
                    if (satelliteOrbit.equalsIgnoreCase(orbits.get(i))) {
                        satelliteFactIndex = j;
                        break;
                    }
                }
                candidateSatellites.add(satelliteFactIndex);
                satellitesFactToList.put(satelliteFactIndex, i);
            }
        }

        if (!candidateSatellites.isEmpty()) {
            int numSatellitesModified = 0;
            int satelliteIndex;
            while (numSatellitesModified < ySatellites) {
                satelliteIndex = PRNG.nextInt(candidateSatellites.size());
                int candidateSatelliteIndex = candidateSatellites.get(satelliteIndex);
                int candidateSatelliteListIndex = satellitesFactToList.get(candidateSatelliteIndex);
                try {
                    ValueVector satelliteInstruments = satellites.get(candidateSatelliteListIndex).getSlotValue("instruments").listValue(rete.getGlobalContext());
                    ArrayList<String> instrumentsNotInSatellite = getInstrumentsNotPresent(satelliteInstruments, params.getInstrumentList(), rete);
                    int numInstrumentsAdded = 0;
                    double totalAddedInstrumentsMass = 0.0;
                    while (numInstrumentsAdded < xInstruments) {
                        ArrayList<String> candidateAdditionInstruments = getCandidateInstrumentsToAdd(instrumentsNotInSatellite, totalAddedInstrumentsMass, satellites, candidateSatelliteListIndex, instrumentFacts, res, rete);
                        if (!candidateAdditionInstruments.isEmpty()) {
                            int instrumentToAddIndex = PRNG.nextInt(candidateAdditionInstruments.size());
                            childPayloads.get(candidateSatelliteListIndex).add(candidateAdditionInstruments.get(instrumentToAddIndex));
                            Fact instrumentFact = getInstrumentFact(instrumentFacts, candidateAdditionInstruments.get(instrumentToAddIndex), rete);
                            totalAddedInstrumentsMass += instrumentFact.getSlotValue("mass#").floatValue(rete.getGlobalContext());
                            numInstrumentsAdded += 1;
                        } else {
                            candidateSatellites.remove(candidateSatelliteListIndex);
                            break;
                        }
                    }
                    if (numInstrumentsAdded > 0) {
                        numSatellitesModified += 1;
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

    private ArrayList<String> getCandidateInstrumentsToAdd (ArrayList<String> absentInstruments, double totalAddedInstrumentMass, ArrayList<Fact> allSatellites, Integer satelliteIndex, ArrayList<Fact> instrumentFacts, Resource res, Rete r) throws JessException {
        ArrayList<String> candidateInstruments = new ArrayList<>();
        for (int i = 0; i < absentInstruments.size(); i++) {
            double satelliteLaunchMass = allSatellites.get(satelliteIndex).getSlotValue("satellite-launch-mass").floatValue(r.getGlobalContext());
            Fact instrumentFact = getInstrumentFact(instrumentFacts, absentInstruments.get(i), r);
            double instrumentMass = instrumentFact.getSlotValue("mass#").floatValue(r.getGlobalContext());
            String launchVehicle = allSatellites.get(satelliteIndex).getSlotValue("launch-vehicle").stringValue(r.getGlobalContext());
            double orbitAltitude = allSatellites.get(satelliteIndex).getSlotValue("orbit-altitude#").floatValue(r.getGlobalContext());
            String orbitType = allSatellites.get(satelliteIndex).getSlotValue("orbit-type").stringValue(r.getGlobalContext());
            String orbitInclination = allSatellites.get(satelliteIndex).getSlotValue("orbit-inclination").stringValue(r.getGlobalContext());
            double launchVehicleCapacity = getLaunchVehicleCapacity(res, r, launchVehicle, orbitType, orbitInclination, orbitAltitude);
            int numberOfLaunches = allSatellites.get(satelliteIndex).getSlotValue("num-launches").intValue(r.getGlobalContext());
            if (satelliteLaunchMass + totalAddedInstrumentMass + (3*instrumentMass) < (numberOfLaunches*launchVehicleCapacity)) {
                candidateInstruments.add(absentInstruments.get(i));
            }
        }
        return candidateInstruments;
    }

    private ArrayList<String> getInstrumentsNotPresent (ValueVector satelliteInstruments, String[] instrumentsList, Rete r) throws JessException {
        ArrayList<String> absentInstruments = new ArrayList<>();
        for (int i = 0; i < instrumentsList.length; i++) {
            boolean instrumentPresent = false;
            for (int j = 0; j < satelliteInstruments.size(); j++) {
                if (instrumentsList[i].equalsIgnoreCase(satelliteInstruments.get(j).stringValue(r.getGlobalContext()))) {
                    instrumentPresent = true;
                    break;
                }
            }
            if (!instrumentPresent) {
                absentInstruments.add(instrumentsList[i]);
            }
        }
        return absentInstruments;
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
