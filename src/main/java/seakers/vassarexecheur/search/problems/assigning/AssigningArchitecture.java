package seakers.vassarexecheur.search.problems.assigning;

import org.moeaframework.core.Solution;
import seakers.architecture.pattern.ArchitecturalDecision;
import seakers.architecture.pattern.Assigning;
import seakers.architecture.pattern.Combining;
import seakers.architecture.Architecture;

import java.util.ArrayList;

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

    public AssigningArchitecture(int[] alternativesForNumberOfSatellites, int numberOfInstruments, int numberOfOrbits, int numberOfObjectives) {
        super(numberOfObjectives, 0, createDecisions(alternativesForNumberOfSatellites, numberOfInstruments, numberOfOrbits));
        this.alternativesForNumberOfSatellites = alternativesForNumberOfSatellites;
        this.alreadyEvaluated = false;
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
    }

    public void setAlreadyEvaluated(boolean alreadyEvaluated) {
        this.alreadyEvaluated = alreadyEvaluated;
    }

    public boolean getAlreadyEvaluated() {
        return this.alreadyEvaluated;
    }

    @Override
    public Solution copy() { return new AssigningArchitecture(this); }
}
