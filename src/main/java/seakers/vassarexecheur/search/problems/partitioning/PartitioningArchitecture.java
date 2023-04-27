package seakers.vassarexecheur.search.problems.partitioning;

import org.apache.commons.lang3.SerializationUtils;
import org.moeaframework.core.Solution;
import seakers.architecture.Architecture;
import seakers.architecture.pattern.ArchitecturalDecision;
import seakers.architecture.pattern.Partitioning;
import seakers.architecture.util.IntegerVariable;
import seakers.engineerserver.search.problems.PartitioningAndAssigning.AssigningPatternCategorical;
import seakers.vassarexecheur.search.problems.assigning.AssigningArchitecture;
import seakers.vassarheur.BaseParams;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.Map;
import java.util.regex.Pattern;

public class PartitioningArchitecture extends Architecture {
    /**
     * Tag used for the assigning decision
     */
    public static final String assignTag = "assigning";
    /**
     * Tag used for the partitioning decision
     */
    public static final String partitionTag = "partitioning";
    private static final long serialVersionUID = 8776271523867135862L;
    private int numberOfInstruments;
    private int numberOfOrbits;
    private BaseParams params;
    private boolean alreadyEvaluated;
    private ArrayList<ArrayList<Double>> operatorParameters;
    private ArrayList<ArrayList<String>> satellitePayloads;
    private ArrayList<String> satelliteOrbits;

    public PartitioningArchitecture(int numberOfInstruments, int numberOfOrbits, int numberOfObjectives, BaseParams params) {
        super(numberOfObjectives, 1, createDecisions(numberOfInstruments, numberOfOrbits));
        this.numberOfInstruments = numberOfInstruments;
        this.numberOfOrbits = numberOfOrbits;
        this.params = params;
        this.alreadyEvaluated = false;
        this.operatorParameters = null;
        this.satellitePayloads = null;
        this.satelliteOrbits = null;
    }

    public PartitioningArchitecture(int numberOfInstruments, int numberOfOrbits, int numberOfObjectives, int numberOfConstraints, BaseParams params) {
        super(numberOfObjectives, 1+numberOfConstraints, createDecisions(numberOfInstruments, numberOfOrbits));
        this.numberOfInstruments = numberOfInstruments;
        this.numberOfOrbits = numberOfOrbits;
        this.params = params;
        this.alreadyEvaluated = false;
        this.operatorParameters = null;
        this.satellitePayloads = null;
        this.satelliteOrbits = null;
    }

    private static ArrayList<ArchitecturalDecision> createDecisions(int numberOfInstruments, int numberOfOrbits) {
        ArrayList<ArchitecturalDecision> out = new ArrayList<>();
        out.add(new Partitioning(numberOfInstruments, partitionTag));
        out.add(new AssigningPatternCategorical(numberOfInstruments, numberOfOrbits, assignTag));
        return out;
    }

    private PartitioningArchitecture(Solution solution) { super(solution); }

    public void setAlreadyEvaluated(boolean alreadyEvaluated) {this.alreadyEvaluated = alreadyEvaluated; }

    public boolean getAlreadyEvaluated() {return this.alreadyEvaluated; }

    public void setOperatorParameters(ArrayList<ArrayList<Double>> operatorParameters) { this.operatorParameters = operatorParameters; }

    public ArrayList<ArrayList<Double>> getOperatorParameters() { return this.operatorParameters; }

    public void setSatellitePayloads(ArrayList<ArrayList<String>> satellitePayloads) { this.satellitePayloads = satellitePayloads; }

    public ArrayList<ArrayList<String>> getSatellitePayloads() { return this.satellitePayloads; }

    public void setSatelliteOrbits(ArrayList<String> satelliteOrbits) { this.satelliteOrbits = satelliteOrbits; }

    public ArrayList<String> getSatelliteOrbits() { return this.satelliteOrbits; }

    @Override
    public Solution copy() { return new PartitioningArchitecture(this); }

    public Solution deepCopy() {
        PartitioningArchitecture copy = (PartitioningArchitecture) this.copy();

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

    public void setVariablesFromPartitionArrays(int[] instrumentPartitions, int[] orbitAssignments) {
        // Populate arch with assigning and partitioning values
        for (int i = 0; i < 2*instrumentPartitions.length; i++) {
            if (i < instrumentPartitions.length) {
                IntegerVariable var = new IntegerVariable(instrumentPartitions[i], 0, instrumentPartitions.length);
                this.setVariable(i, var);
            }
            if (i >= instrumentPartitions.length) {
                IntegerVariable var = new IntegerVariable(orbitAssignments[i - instrumentPartitions.length], -1, numberOfOrbits);
                this.setVariable(i, var);
            }
        }
    }

    public int[] getInstrumentPartitionsFromString(String archString) {
        String[] partitionStrings = archString.split(Pattern.quote("|"),2);
        ArrayList<Integer> instrumentPartitioningArrayList = new ArrayList<>();
        String instrumentPartitionString = partitionStrings[0];
        for (int i = 0; i < instrumentPartitionString.length(); i++) {
            String val = instrumentPartitionString.substring(i,i+1);
            if (!val.equalsIgnoreCase(" ")) {
                instrumentPartitioningArrayList.add(Integer.parseInt(val));
            }
        }
        return instrumentPartitioningArrayList.stream().mapToInt(i->i).toArray();
    }

    public int[] getOrbitAssignmentsFromString(String archString) {
        String[] partitionStrings = archString.split(Pattern.quote("|"),2);
        int[] orbitAssignment = new int[numberOfInstruments];
        Arrays.fill(orbitAssignment, -1);
        String orbitAssignmentString = partitionStrings[1];
        int index = 0;
        for (int i = 0; i < orbitAssignmentString.length(); i++) {
            String val = orbitAssignmentString.substring(i,i+1);
            if (!val.equalsIgnoreCase(" ")) {
                if (val.equalsIgnoreCase("-")) {
                    break;
                } else {
                    orbitAssignment[index] = Integer.parseInt(val);
                    index += 1;
                }
            }
        }
        return orbitAssignment;
    }

    public String getString() {

        int[] instrumentPartitioning = new int[numberOfInstruments];
        int[] orbitAssignment = new int[numberOfInstruments];

        for (int i = 0; i < numberOfInstruments; i++) {
            instrumentPartitioning[i] = ((IntegerVariable)this.getVariable(i)).getValue();
        }

        for (int i = 0; i < numberOfInstruments; i++) {
            orbitAssignment[i] = ((IntegerVariable) this.getVariable(numberOfInstruments + i)).getValue();
        }

        seakers.vassarheur.problems.PartitioningAndAssigning.Architecture arch_abs= new seakers.vassarheur.problems.PartitioningAndAssigning.Architecture(instrumentPartitioning, orbitAssignment, 1, params);

        return arch_abs.toString("");
    }
}
