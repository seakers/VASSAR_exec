package seakers.vassartest.search.problems.Assigning;


import org.moeaframework.core.Solution;
import org.moeaframework.core.variable.BinaryVariable;
import seakers.architecture.Architecture;
import seakers.architecture.pattern.ArchitecturalDecision;
import seakers.architecture.pattern.Assigning;
import seakers.architecture.pattern.Combining;

import java.util.*;


public class GigaArchitecture extends Solution {

    public String bit_string;
    public boolean already_evaluated;

    // Random Generators
    public Random rand = new Random();
    public static Random rand_stat = new Random();

    public GigaArchitecture(String bit_string){
        super(1, 2, 0);
        this.bit_string = bit_string;
        this.already_evaluated = false;

        BinaryVariable var = new BinaryVariable(0);
        this.setVariable(0, var);
    }

    protected GigaArchitecture(Solution solution) {
        super(solution);
        if (!(solution instanceof GigaArchitecture)) {
            throw new ClassCastException("Solution is not an instance of Architecture");
        } else {
            GigaArchitecture arch = (GigaArchitecture)solution;
            this.bit_string = arch.bit_string;
            this.already_evaluated = arch.already_evaluated;
        }
    }



    @Override
    public String toString(){
        return this.bit_string;
    }

    @Override
    public int hashCode(){
        return this.bit_string.hashCode();
    }

    @Override
    public boolean equals(Object obj){
        if (obj == null) {
            return false;
        } else if (this.getClass() != obj.getClass()) {
            return false;
        } else {
            GigaArchitecture other = (GigaArchitecture)obj;
            return this.bit_string.equalsIgnoreCase(other.bit_string);
        }
    }
}
