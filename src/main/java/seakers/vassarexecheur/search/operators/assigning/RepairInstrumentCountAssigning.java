package seakers.vassarexecheur.search.operators.assigning;

import jess.Fact;
import jess.JessException;
import jess.Rete;
import jess.ValueVector;
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

public class RepairInstrumentCountAssigning implements Variation {
    /**
     * The number of instruments to remove from each satellite
     */
    private final int xInstruments;

    /**
     * The number of satellites to remove instruments from
     */
    private final int ySatellites;

    private final double instrCountThreshold;

    /**
     * Assigning Problem used to evaluate architectures
     */
    private AssigningProblem problem;

    /**
     * Resource pool to obtain matlab functions
     */
    private final ResourcePool resourcePool;

    private final ArchitectureEvaluator evaluator;

    private final BaseParams params;

    public RepairInstrumentCountAssigning(int xInstruments, int ySatellites, double instrCountThreshold, AssigningProblem problem, ResourcePool resourcePool, ArchitectureEvaluator evaluator, BaseParams params) {
        this.xInstruments = xInstruments;
        this.ySatellites = ySatellites;
        this.instrCountThreshold = instrCountThreshold;
        this.problem = problem;
        this.resourcePool = resourcePool;
        this.evaluator = evaluator;
        this.params = params;
    }

    @Override
    public int getArity() {
        return 1;
    }

    @Override
    public Solution[] evolve(Solution[] solutions) {
        AssigningArchitecture parent = (AssigningArchitecture) solutions[0];
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

        ArrayList<String> allOrbits = addEmptyOrbits(orbits, params.getOrbitList());
        ArrayList<ArrayList<String>> allPayloads = addEmptyPayloads(payloads, orbits, params.getOrbitList());

        ArrayList<ArrayList<String>> childPayloads = new ArrayList<>(allPayloads);
        ArrayList<String> childOrbits = new ArrayList<>(allOrbits);

        if (problem.getNumberOfInstruments(parent) > instrCountThreshold) {
            for (int i = 0; i < ySatellites; i++) {
                int randomSatIndex = PRNG.nextInt(childOrbits.size());
                if (childPayloads.get(randomSatIndex).isEmpty()) {
                    i--;
                    continue;
                } else {
                    for (int j = 0; j < xInstruments; j++) {
                        int randomInstrumentIndex = PRNG.nextInt(childPayloads.get(randomSatIndex).size());
                        childPayloads.get(randomSatIndex).remove(randomInstrumentIndex);
                        if (childPayloads.get(randomSatIndex).isEmpty()) {
                            break;
                        }
                    }
                }
                if (getNonEmptySatellites(childPayloads) < (ySatellites - i)) {
                    break;
                }
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

    private int getNonEmptySatellites(ArrayList<ArrayList<String>> completePayloads) {
        int numNonEmptySatellites = 0;
        for (int i = 0; i < completePayloads.size(); i++) {
            if (!completePayloads.get(i).isEmpty()) {
                numNonEmptySatellites++;
            }
        }
        return numNonEmptySatellites;
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
