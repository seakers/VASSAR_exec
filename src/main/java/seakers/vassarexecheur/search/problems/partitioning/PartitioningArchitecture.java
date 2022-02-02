package seakers.vassarexecheur.search.problems.partitioning;

import org.moeaframework.core.Solution;
import seakers.architecture.Architecture;
import seakers.architecture.pattern.ArchitecturalDecision;
import seakers.architecture.pattern.Partitioning;
import seakers.engineerserver.search.problems.PartitioningAndAssigning.AssigningPatternCategorical;

import java.util.ArrayList;

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
    private boolean alreadyEvaluated;
    private ArrayList<ArrayList<Double>> operatorParameters;
    private ArrayList<ArrayList<String>> satellitePayloads;
    private ArrayList<String> satelliteOrbits;

    public PartitioningArchitecture(int numberOfInstruments, int numberOfOrbits, int numberOfObjectives) {
        super(numberOfObjectives, 1, createDecisions(numberOfInstruments, numberOfOrbits));
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
}
