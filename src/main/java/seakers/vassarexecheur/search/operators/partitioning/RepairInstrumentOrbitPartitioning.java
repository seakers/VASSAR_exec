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
import seakers.vassarexecheur.search.problems.partitioning.PartitioningArchitecture;
import seakers.vassarexecheur.search.problems.partitioning.PartitioningProblem;
import seakers.vassarheur.BaseParams;
import seakers.vassarheur.QueryBuilder;
import seakers.vassarexecheur.search.problems.assigning.AssigningArchitecture;
import seakers.vassarheur.Resource;
import seakers.vassarheur.ResourcePool;
import seakers.vassarheur.architecture.AbstractArchitecture;
import seakers.vassarheur.problems.PartitioningAndAssigning.ArchitectureEvaluator;
import seakers.vassarheur.problems.PartitioningAndAssigning.ClimateCentricPartitioningParams;
import seakers.vassarheur.problems.PartitioningAndAssigning.Architecture;
import seakers.vassarheur.utils.MatlabFunctions;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

/**
 * Moves incompatible instruments from one satellite to another. Incompatibility involves a chemistry concept
 * instrument in a PM orbit, passive illumination instrument in a DD orbit or a slant geometry instrument in an orbit
 * with altitude 400km or less. The number of moves is specified by the user.
 * Class modified from the original RepairInstrumentOrbit.java class in the EOSS-AIAA repo by nozomihitomi
 *
 * @author roshansuresh
 */

public class RepairInstrumentOrbitPartitioning implements Variation {

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

    /**
     * Assigning Problem used to evaluate architectures
     */
    private PartitioningProblem problem;

    private final BaseParams params;

    public RepairInstrumentOrbitPartitioning(int numChanges, ResourcePool resourcePool, ArchitectureEvaluator evaluator, BaseParams params, PartitioningProblem problem) {
        this.numberOfChanges = numChanges;
        this.resourcePool = resourcePool;
        this.evaluator = evaluator;
        this.params = params;
        this.problem = problem;
        //this.pprng = new PRNG();
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
                int selectedMoveIndex = PRNG.nextInt(candidateInstrumentMoves.size());
                ArrayList<Integer> selectedMove = candidateInstrumentMoves.get(selectedMoveIndex);
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
                candidateInstrumentMoves.remove(selectedMoveIndex);
            }
            this.resourcePool.freeResource(res);
            PartitioningArchitecture child = getArchitectureFromPayloadsAndOrbits(payloads, orbits);

            return new Solution[]{child};
        } else {
            return new Solution[]{parent};
        }
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
}
