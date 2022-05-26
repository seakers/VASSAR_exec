package seakers.vassarexecheur.search.operators.assigning;

import jess.Fact;
import jess.JessException;
import jess.ValueVector;
import jxl.Workbook;
import jxl.read.biff.BiffException;
import org.moeaframework.core.PRNG;
import org.moeaframework.core.Solution;
import org.moeaframework.core.Variation;
import jess.Rete;
import org.moeaframework.core.variable.BinaryVariable;
import org.moeaframework.core.variable.RealVariable;
import seakers.vassarexecheur.search.problems.assigning.AssigningProblem;
import seakers.vassarheur.*;
import seakers.vassarexecheur.search.problems.assigning.AssigningArchitecture;
import seakers.vassarheur.architecture.AbstractArchitecture;
import seakers.vassarheur.problems.Assigning.Architecture;
import seakers.vassarheur.problems.Assigning.ArchitectureEvaluator;
import seakers.vassarheur.problems.Assigning.ClimateCentricAssigningParams;
import seakers.vassarheur.utils.MatlabFunctions;

import java.io.File;
import java.io.IOException;
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

public class RepairInstrumentOrbitAssigning implements Variation {

    /**
     * Number of instrument moves
     */
    private int numberOfChanges;

    /**
     * Resource pool to obtain matlab functions
     */
    private final ResourcePool resourcePool;

    private final ArchitectureEvaluator evaluator;

    //private final PRNG pprng;

    private final BaseParams params;

    /**
     * Assigning Problem used to evaluate architectures
     */
    private AssigningProblem problem;

    /**
     * Toggles movement of removed instruments into other spacecraft
     */
    private boolean moveInstruments;

    public RepairInstrumentOrbitAssigning(int numChanges, ResourcePool resourcePool, ArchitectureEvaluator evaluator, BaseParams params, AssigningProblem problem) {
        this.numberOfChanges = numChanges;
        this.resourcePool = resourcePool;
        this.evaluator = evaluator;
        this.params = params;
        this.problem = problem;
        this.moveInstruments = true;
        //this.pprng = new PRNG();
    }

    public RepairInstrumentOrbitAssigning(int numChanges, ResourcePool resourcePool, ArchitectureEvaluator evaluator, BaseParams params, AssigningProblem problem, boolean moveInstruments) {
        this.numberOfChanges = numChanges;
        this.resourcePool = resourcePool;
        this.evaluator = evaluator;
        this.params = params;
        this.problem = problem;
        this.moveInstruments = moveInstruments;
        //this.pprng = new PRNG();
    }

    @Override
    public int getArity() {
        return 1;
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

        ArrayList<String> instrumentList = new ArrayList<>(Arrays.asList(params.getInstrumentList()));
        //ArrayList<String> orbitList = new ArrayList<>(Arrays.asList(params.getOrbitList()));

        ArrayList<Fact> instrumentFacts = queryBuilder.makeQuery("DATABASE::Instrument");

        // Find the valid instrument moves that can be made
        ArrayList<ArrayList<Integer>> candidateInstrumentMoves = getValidMoves(rete, satellites, instrumentFacts, instrumentList);

        // Make the required moves at random
        int numberOfMoves = 0;
        while ((numberOfMoves < numberOfChanges) && (candidateInstrumentMoves.size() != 0)) {
            //System.out.println("Move Count: " + numberOfMoves);
            int selectedMoveIndex = PRNG.nextInt(candidateInstrumentMoves.size());
            ArrayList<Integer> selectedMove = candidateInstrumentMoves.get(selectedMoveIndex);
            //System.out.println("Selected Move: " + selectedMove);
            String selectedRemovalOrbit = orbits.get(selectedMove.get(1));

            // Remove instrument from selected satellite
            String selectedInstrument = instrumentList.get(selectedMove.get(0));
            ArrayList<String> removalOrbitPayload = payloads.get(selectedMove.get(1));
            removalOrbitPayload.remove(selectedInstrument);

            // Update payloads
            payloads.set(selectedMove.get(1), removalOrbitPayload);

            if (moveInstruments) {
                // Add instrument to selected satellite
                ArrayList<String> additionOrbitPayload = payloads.get(selectedMove.get(2));
                additionOrbitPayload.add(selectedInstrument);

                // Update payloads
                payloads.set(selectedMove.get(2), additionOrbitPayload);
            }

            numberOfMoves += 1;
            candidateInstrumentMoves.remove(selectedMoveIndex);
        }
        this.resourcePool.freeResource(res);
        AssigningArchitecture child = getArchitectureFromPayloadsAndOrbits(payloads, orbits);

        return new Solution[]{child};
    }

    private ArrayList<ArrayList<Integer>> getValidMoves (Rete r, ArrayList<Fact> allSatellites, ArrayList<Fact> allInstruments, ArrayList<String> listOfInstruments) {
        // Valid move -> ArrayList {instrument_index (in instrumentsList), index_of_satellite_to_remove_from, index_of_satellite_to_add_to}

        ArrayList<ArrayList<Integer>> validMoves = new ArrayList<>();

        for (int i = 0; i < allSatellites.size(); i++) {
            Fact currentSatellite = allSatellites.get(i);
            try {
                double orbitAltitude = currentSatellite.getSlotValue("orbit-altitude#").floatValue(r.getGlobalContext());
                String orbitRAAN = currentSatellite.getSlotValue("orbit-RAAN").stringValue(r.getGlobalContext());
                ValueVector satelliteInstruments = currentSatellite.getSlotValue("instruments").listValue(r.getGlobalContext());

                for (int j = 0; j < satelliteInstruments.size(); j++) {
                    String currentInstrument  = satelliteInstruments.get(j).stringValue(r.getGlobalContext());
                    Fact instrumentFact = getInstrumentFact(allInstruments, currentInstrument, r);
                    String instrumentConcept = instrumentFact.getSlotValue("Concept").stringValue(r.getGlobalContext());
                    String instrumentIllumination = instrumentFact.getSlotValue("Illumination").stringValue(r.getGlobalContext());
                    String instrumentGeometry = instrumentFact.getSlotValue("Geometry").stringValue(r.getGlobalContext());

                    if (chemistryInstrumentInPMOrbit(instrumentConcept, orbitRAAN)) { // If a chemistry concept instrument is present in a PM orbit
                        for (int k = 0; k < allSatellites.size(); k++) {
                            if (k == i) {
                                continue;
                            }
                            Fact otherSatellite = allSatellites.get(k);
                            String otherOrbitRAAN = otherSatellite.getSlotValue("orbit-RAAN").stringValue(r.getGlobalContext());
                            ArrayList<String> otherSatellitePayload = new ArrayList<>();
                            ValueVector otherSatelliteInstruments = otherSatellite.getSlotValue("instruments").listValue(r.getGlobalContext());
                            for (int m = 0; m < otherSatelliteInstruments.size(); m++) {
                                otherSatellitePayload.add(otherSatelliteInstruments.get(m).stringValue(r.getGlobalContext()));
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
                            String otherOrbitRAAN = otherSatellite.getSlotValue("orbit-RAAN").stringValue(r.getGlobalContext());
                            ArrayList<String> otherSatellitePayload = new ArrayList<>();
                            ValueVector otherSatelliteInstruments = otherSatellite.getSlotValue("instruments").listValue(r.getGlobalContext());
                            for (int m = 0; m < otherSatelliteInstruments.size(); m++) {
                                otherSatellitePayload.add(otherSatelliteInstruments.get(m).stringValue(r.getGlobalContext()));
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
                            double otherOrbitAltitude = otherSatellite.getSlotValue("orbit-altitude#").floatValue(r.getGlobalContext());
                            ArrayList<String> otherSatellitePayload = new ArrayList<>();
                            ValueVector otherSatelliteInstruments = otherSatellite.getSlotValue("instruments").listValue(r.getGlobalContext());
                            for (int m = 0; m < otherSatelliteInstruments.size(); m++) {
                                otherSatellitePayload.add(otherSatelliteInstruments.get(m).stringValue(r.getGlobalContext()));
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
