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

import java.util.ArrayList;
import java.util.Arrays;

/**
 * Moves incompatible instruments from one satellite to another. Incompatibility involves a chemistry concept
 * instrument in a PM orbit, passive illumination instrument in a DD orbit or a slant geometry instrument in an orbit
 * with altitude 400km or less. The number of moves is specified by the user.
 * Class modified from the original RepairInstrumentOrbit.java class in the EOSS-AIAA repo by nozomihitomi
 *
 * @author roshansuresh
 */

public class RepairInstrumentOrbit implements Variation {

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

    public RepairInstrumentOrbit(int numChanges, ResourcePool resourcePool, ArchitectureEvaluator evaluator, Rete r, QueryBuilder qb, BaseParams params) {
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
        return 1;
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

        ArrayList<String> instrumentList = new ArrayList<>(Arrays.asList(params.getInstrumentList()));
        //ArrayList<String> orbitList = new ArrayList<>(Arrays.asList(params.getOrbitList()));

        ArrayList<Fact> instrumentFacts = queryBuilder.makeQuery("CAPABILITIES::Manifested-instrument");

        // Find the valid instrument moves that can be made
        ArrayList<ArrayList<Integer>> candidateInstrumentMoves = getValidMoves(satellites, instrumentFacts, instrumentList);

        // Make the required moves at random
        int numberOfMoves = 0;
        while ((numberOfMoves < numberOfChanges) && (candidateInstrumentMoves.size() != 0)) {
            ArrayList<Integer> selectedMove = candidateInstrumentMoves.get(pprng.nextInt(candidateInstrumentMoves.size()));
            String selectedRemovalOrbit = orbits.get(selectedMove.get(1));

            // Remove instrument from selected satellite
            String selectedInstrument = instrumentList.get(selectedMove.get(0));
            ArrayList<String> removalOrbitPayload = payloads.get(selectedMove.get(1));
            removalOrbitPayload.remove(selectedInstrument);

            // Update payloads
            payloads.set(selectedMove.get(1), removalOrbitPayload);

            // Add instrument to selected satellite
            ArrayList<String> additionOrbitPayload = payloads.get(selectedMove.get(2));
            additionOrbitPayload.add(selectedInstrument);

            // Update payloads
            payloads.set(selectedMove.get(2), additionOrbitPayload);

            numberOfMoves += 1;
        }

        AssigningArchitecture child = getArchitectureFromPayloadsAndOrbits(payloads, orbits);

        return new Solution[]{child};
    }

    private ArrayList<ArrayList<Integer>> getValidMoves (ArrayList<Fact> allSatellites, ArrayList<Fact> allInstruments, ArrayList<String> listOfInstruments) {
        // Valid move -> ArrayList {instrument_index (in instrumentsList), index_of_satellite_to_remove_from, index_of_satellite_to_add_to}

        ArrayList<ArrayList<Integer>> validMoves = new ArrayList<>();

        for (int i = 0; i < allSatellites.size(); i++) {
            Fact currentSatellite = allSatellites.get(i);
            try {
                double orbitAltitude = currentSatellite.getSlotValue("orbit-altitude#").floatValue(rete.getGlobalContext());
                String orbitRAAN = currentSatellite.getSlotValue("orbit-RAAN").stringValue(rete.getGlobalContext());
                ValueVector satelliteInstruments = currentSatellite.getSlotValue("instruments").listValue(rete.getGlobalContext());

                for (int j = 0; j < satelliteInstruments.size(); j++) {
                    String currentInstrument  = satelliteInstruments.get(j).stringValue(rete.getGlobalContext());
                    Fact instrumentFact = getInstrumentFact(allInstruments, currentInstrument, rete);
                    String instrumentConcept = instrumentFact.getSlotValue("Concept").stringValue(rete.getGlobalContext());
                    String instrumentIllumination = instrumentFact.getSlotValue("Illumination").stringValue(rete.getGlobalContext());
                    String instrumentGeometry = instrumentFact.getSlotValue("Geometry").stringValue(rete.getGlobalContext());

                    if (chemistryInstrumentInPMOrbit(instrumentConcept, orbitRAAN)) { // If a chemistry concept instrument is present in a PM orbit
                        for (int k = 0; k < allSatellites.size(); k++) {
                            if (k == i) {
                                continue;
                            }
                            Fact otherSatellite = allSatellites.get(k);
                            String otherOrbitRAAN = otherSatellite.getSlotValue("orbit-RAAN").stringValue(rete.getGlobalContext());
                            ArrayList<String> otherSatellitePayload = new ArrayList<>();
                            ValueVector otherSatelliteInstruments = otherSatellite.getSlotValue("instruments").listValue(rete.getGlobalContext());
                            for (int m = 0; m < otherSatelliteInstruments.size(); m++) {
                                otherSatellitePayload.add(otherSatelliteInstruments.get(j).stringValue(rete.getGlobalContext()));
                            }
                            if (!chemistryInstrumentInPMOrbit(instrumentConcept, otherOrbitRAAN) && !otherSatellitePayload.contains(currentInstrument)) {
                                ArrayList<Integer> validMove = new ArrayList<>();
                                validMove.add(listOfInstruments.indexOf(currentInstrument));
                                validMove.add(i);
                                validMove.add(k);
                                validMoves.add(validMove);
                            }
                        }
                    }

                    if (passiveInstrumentInDDOribt(instrumentIllumination, orbitRAAN)) { // If a passive illumination instrument is present in DD orbit
                        for (int k = 0; k < allSatellites.size(); k++) {
                            if (k == i) {
                                continue;
                            }
                            Fact otherSatellite = allSatellites.get(k);
                            String otherOrbitRAAN = otherSatellite.getSlotValue("orbit-RAAN").stringValue(rete.getGlobalContext());
                            ArrayList<String> otherSatellitePayload = new ArrayList<>();
                            ValueVector otherSatelliteInstruments = otherSatellite.getSlotValue("instruments").listValue(rete.getGlobalContext());
                            for (int m = 0; m < otherSatelliteInstruments.size(); m++) {
                                otherSatellitePayload.add(otherSatelliteInstruments.get(j).stringValue(rete.getGlobalContext()));
                            }
                            if (!passiveInstrumentInDDOribt(instrumentIllumination, otherOrbitRAAN) && !otherSatellitePayload.contains(currentInstrument)) {
                                ArrayList<Integer> validMove = new ArrayList<>();
                                validMove.add(listOfInstruments.indexOf(currentInstrument));
                                validMove.add(i);
                                validMove.add(k);
                                validMoves.add(validMove);
                            }
                        }
                    }

                    if (slantInstrumentInLowAltitudeOrbit(instrumentGeometry, orbitAltitude)) { // If a slant geometry instrument is present in a low altitude orbit
                        for (int k = 0; k < allSatellites.size(); k++) {
                            if (k == i) {
                                continue;
                            }
                            Fact otherSatellite = allSatellites.get(k);
                            double otherOrbitAltitude = otherSatellite.getSlotValue("orbit-altitude#").floatValue(rete.getGlobalContext());
                            ArrayList<String> otherSatellitePayload = new ArrayList<>();
                            ValueVector otherSatelliteInstruments = otherSatellite.getSlotValue("instruments").listValue(rete.getGlobalContext());
                            for (int m = 0; m < otherSatelliteInstruments.size(); m++) {
                                otherSatellitePayload.add(otherSatelliteInstruments.get(j).stringValue(rete.getGlobalContext()));
                            }
                            if (!slantInstrumentInLowAltitudeOrbit(instrumentGeometry, otherOrbitAltitude) && !otherSatellitePayload.contains(currentInstrument)) {
                                ArrayList<Integer> validMove = new ArrayList<>();
                                validMove.add(listOfInstruments.indexOf(currentInstrument));
                                validMove.add(i);
                                validMove.add(k);
                                validMoves.add(validMove);
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

    private boolean chemistryInstrumentInPMOrbit (String instrConcept, String orbRAAN) {
        return (orbRAAN.equals("PM") && instrConcept.contains("chemistry"));
    }

    private boolean passiveInstrumentInDDOribt (String instrIllumination, String orbRAAN) {
        return (orbRAAN.equals("DD") && instrIllumination.equals("Passive"));
    }

    private boolean slantInstrumentInLowAltitudeOrbit (String instrGeometry, double orbAlt) {
        return (instrGeometry.equals("slant") && (orbAlt <= 400.0));
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
}
