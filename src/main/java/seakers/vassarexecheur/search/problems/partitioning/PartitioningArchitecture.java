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

    public PartitioningArchitecture(int numberOfInstruments, int numberOfOrbits, int numberOfObjectives) {
        super(numberOfObjectives, 1, createDecisions(numberOfInstruments, numberOfOrbits));
        this.alreadyEvaluated = false;
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

    @Override
    public Solution copy() { return new PartitioningArchitecture(this); }
}
