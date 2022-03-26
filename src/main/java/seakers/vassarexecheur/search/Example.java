package seakers.vassarexecheur.search;

import java.util.ArrayList;

public class Example {

    public static void main(String[] args) {
        ArrayList<ArrayList<String>> arrays1 = new ArrayList<>();

        int numRows = 3;
        int numCols = 5;

        for (int i = 0; i < numRows; i++) {
            ArrayList<String> array = new ArrayList<>();
            for (int j = 0; j < numCols; j++) {
                array.add(Integer.toString((i+1)*(j+1)));
            }
            arrays1.add(array);
        }
        System.out.println("Full Array");
        System.out.println(arrays1);

        //System.out.println("Truncated Array 1");
        //arrays1.remove(0);
        //System.out.println(arrays1);

        System.out.println("Truncated Array 1");
        arrays1.get(0).remove("3");
        System.out.println(arrays1);

        //System.out.println("Truncated Array 2");
        //ArrayList<String> arrayTrunc = arrays1.get(0);
        //arrayTrunc.remove(1);
        //arrays1.set(0,arrayTrunc);
        //System.out.println(arrays1);

    }
}
