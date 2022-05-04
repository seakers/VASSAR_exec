package seakers.vassarexecheur.search.operators.partitioning;

import org.moeaframework.core.Solution;
import org.moeaframework.core.Variation;
import org.moeaframework.core.PRNG;
import seakers.architecture.Architecture;
import seakers.architecture.util.IntegerVariable;
import seakers.vassarheur.BaseParams;
import seakers.vassarexecheur.search.problems.partitioning.PartitioningArchitecture;

import java.util.Arrays;
import java.util.ArrayList;

public class PartitioningCrossover implements Variation {

    protected double probability;
    protected BaseParams params;

    public PartitioningCrossover(double probability, BaseParams params){
        this.probability = probability;
        this.params = params;
    }

    @Override
    public int getArity() {
        return 2;
    }

    @Override
    public Solution[] evolve(Solution[] parents) {
        Solution[] out = new Solution[2];

        Architecture arch1 = (PartitioningArchitecture) parents[0];
        Architecture arch2 = (PartitioningArchitecture) parents[1];
        int[] intVars1 = getIntVariables(arch1);
        int[] intVars2 = getIntVariables(arch2);
        int[] partitioning1 = Arrays.copyOfRange(intVars1, 0, params.getNumInstr());
        int[] partitioning2 = Arrays.copyOfRange(intVars2, 0, params.getNumInstr());
        int[] assigning1 = Arrays.copyOfRange(intVars1, params.getNumInstr(), 2 * params.getNumInstr());
        int[] assigning2 = Arrays.copyOfRange(intVars2, params.getNumInstr(), 2 * params.getNumInstr());

        Architecture newArch1 = new PartitioningArchitecture(params.getNumInstr(), params.getNumOrbits(), 2, params);
        Architecture newArch2 = new PartitioningArchitecture(params.getNumInstr(), params.getNumOrbits(), 2, params);
        int[] newPartitioning1 = new int[partitioning1.length];
        int[] newPartitioning2 = new int[partitioning2.length];
        int[] newAssigning1 = new int[assigning1.length];
        int[] newAssigning2 = new int[assigning2.length];

        int split = PRNG.nextInt(partitioning1.length); // Single-point crossover

        int[] orbitAssigned1 = new int[partitioning1.length];
        for(int i = 0; i < partitioning1.length; i++){
            int sat = partitioning1[i];
            int orb = assigning1[sat];
            orbitAssigned1[i] = orb;
        }

        int[] orbitAssigned2 = new int[partitioning2.length];
        for(int i = 0; i < partitioning2.length; i++){
            int sat = partitioning2[i];
            int orb = assigning2[sat];
            orbitAssigned2[i] = orb;
        }

        ArrayList<int[]> orbitAssignedSwapped = swapSubarrays(orbitAssigned1, orbitAssigned2, split);

        int[] orbitAssignmentInfo1 = orbitAssignedSwapped.get(0);
        int[] orbitAssignmentInfo2 = orbitAssignedSwapped.get(1);

        ArrayList<int[]> temp1 = extractPartitioningAndAssigning(orbitAssignmentInfo1);
        ArrayList<int[]> temp2 = extractPartitioningAndAssigning(orbitAssignmentInfo2);

        newPartitioning1 = temp1.get(0);
        newAssigning1 = temp1.get(1);
        newPartitioning2 = temp2.get(0);
        newAssigning2 = temp2.get(1);

        int[] newIntVars1 = new int[partitioning1.length + assigning1.length];
        int[] newIntVars2 = new int[partitioning1.length + assigning1.length];
        for(int i = 0; i < newPartitioning1.length;i ++){
            newIntVars1[i] = newPartitioning1[i];
            newIntVars2[i] = newPartitioning2[i];
        }
        for(int i = 0; i < newAssigning1.length;i ++){
            newIntVars1[i + newPartitioning1.length] = newAssigning1[i];
            newIntVars2[i + newPartitioning1.length] = newAssigning2[i];
        }
        setIntVariables(newArch1, newIntVars1);
        setIntVariables(newArch2, newIntVars2);

        out[0] = newArch1;
        out[1] = newArch2;
        return out;
    }

    private ArrayList<int[]> extractPartitioningAndAssigning(int[] assignedOrbit){
        int[] partitioning = new int[assignedOrbit.length];
        int[] assigning = new int[assignedOrbit.length];

        for(int i = 0; i < assignedOrbit.length; i++){
            partitioning[i] = -1;
            assigning[i] = -1;
        }

        int satIndex = 0;
        ArrayList<Integer> orbits = new ArrayList<>();

        int currentOrb = assignedOrbit[0];
        orbits.add(currentOrb);
        for(int i = 0; i < assignedOrbit.length; i++){
            int orb = assignedOrbit[i];

            if (currentOrb == orb) {
                partitioning[i] = satIndex;
            }
            else {
                currentOrb = orb;
                orbits.add(currentOrb);
                satIndex++;
                partitioning[i] = satIndex;
            }
        }

        for (int j = 0; j < orbits.size(); j++) {
            assigning[j] = orbits.get(j);
        }

        ArrayList<int[]> out = new ArrayList<>();
        out.add(partitioning);
        out.add(assigning);
        return out;
    }

    private ArrayList<int[]> swapSubarrays(int[] arr1, int[] arr2, int split){
        int[] out1 = new int[arr1.length];
        int[] out2 = new int[arr2.length];
        for(int i = 0; i < arr1.length; i++){
            if(i < split){
                out1[i] = arr1[i];
                out2[i] = arr2[i];
            }else{
                out2[i] = arr1[i];
                out1[i] = arr2[i];
            }
        }
        ArrayList<int[]> out = new ArrayList<>();
        out.add(out1);
        out.add(out2);
        return out;
    }

    private int[] getIntVariables(Architecture arch){
        int[] out = new int[arch.getNumberOfVariables()];
        for(int i = 0; i < out.length; i++){
            out[i] = ((IntegerVariable) arch.getVariable(i)).getValue();
        }
        return out;
    }

    private void setIntVariables(Architecture arch, int[] values){
        int[] out = new int[arch.getNumberOfVariables()];
        for(int i = 0; i < out.length; i++){
            ((IntegerVariable) arch.getVariable(i)).setValue(values[i]);
        }
    }
}
