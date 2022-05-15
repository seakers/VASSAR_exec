package seakers.vassarexecheur.search.problems.assigning;

import org.apache.commons.lang3.SerializationUtils;
import org.moeaframework.core.Solution;
import org.moeaframework.core.variable.BinaryVariable;
import org.moeaframework.core.variable.RealVariable;
import seakers.architecture.pattern.ArchitecturalDecision;
import seakers.architecture.pattern.Assigning;
import seakers.architecture.pattern.Combining;
import seakers.architecture.Architecture;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Map;

public class AssigningArchitecture extends Architecture {

    public static final long serialVersionUID = 8776271523866899534L;
    /**
     * Tag used for the assigning decision
     */
    private static final String assignTag = "inst";
    /**
     * Tag used for the combining decision
     */
    private static final String combineTag = "nSat";
    /**
     * The available options of the number of satellites
     */
    private final int[] alternativesForNumberOfSatellites;
    private boolean alreadyEvaluated;
    private ArrayList<ArrayList<Double>> operatorParameters;
    private ArrayList<ArrayList<String>> satellitePayloads;
    private ArrayList<String> satelliteOrbits;

    public AssigningArchitecture(int[] alternativesForNumberOfSatellites, int numberOfInstruments, int numberOfOrbits, int numberOfObjectives) {
        super(numberOfObjectives, 0, createDecisions(alternativesForNumberOfSatellites, numberOfInstruments, numberOfOrbits));
        this.alternativesForNumberOfSatellites = alternativesForNumberOfSatellites;
        this.alreadyEvaluated = false;
        this.operatorParameters = null;
        this.satellitePayloads = null;
        this.satelliteOrbits = null;
    }

    public AssigningArchitecture(int[] alternativesForNumberOfSatellites, int numberOfInstruments, int numberOfOrbits, int numberOfObjectives, int numberOfConstraints) {
        super(numberOfObjectives, numberOfConstraints, createDecisions(alternativesForNumberOfSatellites, numberOfInstruments, numberOfOrbits));
        this.alternativesForNumberOfSatellites = alternativesForNumberOfSatellites;
        this.alreadyEvaluated = false;
        this.operatorParameters = null;
        this.satellitePayloads = null;
        this.satelliteOrbits = null;
    }

    private static ArrayList<ArchitecturalDecision> createDecisions(int[] altnertivesForNumberOfSatellites, int numberOfInstruments, int numberOfOrbits) {
        ArrayList<ArchitecturalDecision> out = new ArrayList<>();
        out.add(new Combining(new int[]{altnertivesForNumberOfSatellites.length}, combineTag));
        out.add(new Assigning(numberOfInstruments, numberOfOrbits, assignTag));
        return out;
    }

    /**
     * makes a copy solution from the input solution
     *
     * @param solution
     */
    private AssigningArchitecture(Solution solution) {
        super(solution);
        this.alternativesForNumberOfSatellites = ((AssigningArchitecture) solution).alternativesForNumberOfSatellites;
        this.operatorParameters = ((AssigningArchitecture) solution).operatorParameters;
        this.satellitePayloads = ((AssigningArchitecture) solution).satellitePayloads;
    }

    public void setAlreadyEvaluated(boolean alreadyEvaluated) {
        this.alreadyEvaluated = alreadyEvaluated;
    }

    public boolean getAlreadyEvaluated() {
        return this.alreadyEvaluated;
    }

    public void setOperatorParameters(ArrayList<ArrayList<Double>> operatorParameters) { this.operatorParameters = operatorParameters; }

    public ArrayList<ArrayList<Double>> getOperatorParameters() { return this.operatorParameters; }

    public void setSatellitePayloads(ArrayList<ArrayList<String>> satellitePayloads) { this.satellitePayloads = satellitePayloads; }

    public ArrayList<ArrayList<String>> getSatellitePayloads() { return this.satellitePayloads; }

    public void setSatelliteOrbits(ArrayList<String> satelliteOrbits) { this.satelliteOrbits = satelliteOrbits; }

    public ArrayList<String> getSatelliteOrbits() { return this.satelliteOrbits; }

    @Override
    public Solution copy() { return new AssigningArchitecture(this); }

    public Solution deepCopy() {
        AssigningArchitecture copy = (AssigningArchitecture) this.copy();

        // Copy Variables
        for (int i = 0; i < this.getNumberOfVariables(); i++) {
            copy.setVariable(i, this.getVariable(i));
        }

        // Copy Objectives and Constraints
        for (int i = 0; i < this.getNumberOfObjectives(); i++) {
            copy.setObjective(i, this.getObjective(i));
        }

        for (int i = 0; i < this.getNumberOfConstraints(); i++) {
            copy.setConstraint(i, this.getConstraint(i));
        }

        // Copy Attributes
        Iterator iterator = this.getAttributes().entrySet().iterator();
        while(iterator.hasNext()) {
            Map.Entry entry = (Map.Entry)iterator.next();
            copy.setAttribute((String)entry.getKey(), SerializationUtils.clone((Serializable)entry.getValue()));
        }
        return copy;
    }

    public void setVariablesfromString(String archString) {
        // Populate arch with bits from architectureString
        RealVariable var0 = new RealVariable(0.0, 0.0, 0.0);
        this.setVariable(0, var0);
        for (int i = 0; i < archString.length(); i++) {
            BinaryVariable var = new BinaryVariable(1);
            String bit = archString.substring(i,i+1);
            if (bit.equalsIgnoreCase("1")) {
                var.set(0, true);
            } else {
                var.set(0, false);
            }
            this.setVariable(i+1, var);
        }
    }

    public String getBitString() {
        String bitString = "";
        for (int i = 1; i < this.getNumberOfVariables(); ++i) {
            bitString += this.getVariable(i).toString();
        }

        return bitString;
    }

}
