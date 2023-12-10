package seakers.vassartest.search.problems.Assigning;

import java.util.*;
import java.util.stream.Collectors;

public class DesignBuilder {

    public static Random rand = new Random();

    public static String randomPartitionBitString(int numInstruments, int numOrbs, int maxInsts){
        Random rand = new Random();

        int numInsts = rand.nextInt(maxInsts) + 1;
        List<Integer> insts = Collections.nCopies(numInsts, -1);
//        System.out.println("--> NUM INSTS: " + numInsts);

        List<List<Integer>> partitions = new ArrayList<>();
        if (numInsts == 1) {
            partitions.add(new ArrayList<>(Arrays.asList(-1)));
        } else {
            int numPartitions = rand.nextInt(numInsts) + 1;
            for (int i = 0; i < numPartitions; i++) {
                partitions.add(new ArrayList<>());
            }
            for (int inst : insts) {
                int randomPartitionIdx = rand.nextInt(numPartitions);
                partitions.get(randomPartitionIdx).add(inst);
            }
            partitions = partitions.stream().filter(p -> !p.isEmpty()).collect(Collectors.toList());
        }
//        System.out.println("--> PARTITIONS: " + partitions.size() + " " + partitions);

        for (List<Integer> partition : partitions) {
            int[] instruments = rand.ints(0, numInstruments).distinct().limit(partition.size()).toArray();
            Arrays.sort(instruments);
            partition.clear();
            for (int inst : instruments) {
                partition.add(inst);
            }
        }
//        System.out.println("--> PARTITIONS FILLED: " + partitions);



        int[] orbits = rand.ints(0, numOrbs).distinct().limit(partitions.size()).toArray();
        Arrays.sort(orbits);
//        for(int i: orbits){
//            System.out.println("--> ORBIT: " + i);
//        }

        int[] design = new int[numOrbs * numInstruments];
        for (int orbNum = 0; orbNum < numOrbs; orbNum++) {
            if (contains(orbits, orbNum)) {
                int partitionIdx = indexOf(orbits, orbNum);
                List<Integer> partition = partitions.get(partitionIdx);
                for (int instNum = 0; instNum < numInstruments; instNum++) {
                    design[orbNum * numInstruments + instNum] = partition.contains(instNum) ? 1 : 0;
                }
            } else {
                for (int instNum = 0; instNum < numInstruments; instNum++) {
                    design[orbNum * numInstruments + instNum] = 0;
                }
            }
        }
        String design_str = "";
        for(int i: design){
            design_str += Integer.toString(i);
        }
//        System.out.println("--> 1 BITS: " + counter + " " + design_str.length());
        return design_str;
    }

    public static boolean contains(int[] arr, int target) {
        for (int val : arr) {
            if (val == target) {
                return true;
            }
        }
        return false;
    }

    public static int indexOf(int[] arr, int target) {
        for (int i = 0; i < arr.length; i++) {
            if (arr[i] == target) {
                return i;
            }
        }
        return -1;
    }











    public static String randomBitString(int num_insts, int num_orbs, int max_bits){
        int length = num_insts * num_orbs;
        int[] bitArray = new int[length];
        for (int i = 0; i < max_bits; i++) {
            int index = rand.nextInt(length);
            bitArray[index] = rand.nextInt(2);
        }

        StringBuilder bitString = new StringBuilder();
        int oneCount = 0;
        for (int bit : bitArray) {
            bitString.append(bit);
            if (bit == 1) {
                oneCount++;
            }
        }

        if (oneCount == 0) {
            int index = rand.nextInt(length);
            bitArray[index] = 1;
            bitString.setCharAt(index, '1');
        }

        return bitString.toString();
    }


    public static String randomBitStringFixed(int num_insts, int num_orbs, int num_bits){
        int length = num_insts * num_orbs;
        char[] bitString = new char[length];
        Arrays.fill(bitString, '0');
        Random rand = new Random();
        HashSet<Integer> flippedIndices = new HashSet<>();
        while (flippedIndices.size() < num_bits) {
            int randomIndex = GigaArchitecture.rand_stat.nextInt(length);
            if (!flippedIndices.contains(randomIndex)) {
                bitString[randomIndex] = '1';
                flippedIndices.add(randomIndex);
            }
        }
        return new String(bitString);
    }

    public static String randomBitStringOrbs(int num_insts, int num_orbs, int used_orbs){
        String bit_string = "";

        ArrayList<Integer> orb_idxs = new ArrayList<>();
        for(int x = 0; x < used_orbs; x++){
            orb_idxs.add(GigaArchitecture.rand_stat.nextInt(num_orbs));
        }

        for(int x = 0; x < num_orbs; x++){
            int inst_idx = GigaArchitecture.rand_stat.nextInt(num_insts);
            for(int y = 0; y < num_insts; y++){
                if(y == inst_idx && orb_idxs.contains(x)){
                    bit_string += "1";
                }
                else{
                    bit_string += "0";
                }
            }
        }
        return bit_string;
    }

    public static String specificIdxBitString(int num_insts, int num_orbs, ArrayList<Integer> indices){
        String bit_string = "";
        int count = 0;
        for(int x = 0; x < num_orbs; x++) {
            for (int y = 0; y < num_insts; y++) {
                if (indices.contains(count)) {
                    bit_string += "1";
                } else {
                    bit_string += "0";
                }
                count++;
            }
        }
        return bit_string;
    }


    // 74, 140, 163 for np orbit
    public static String specificOrbInstBitString(int num_insts, int num_orbs, int orb_idx, int inst_idx){
        String bit_string = "";
        for(int x = 0; x < num_orbs; x++) {
            for (int y = 0; y < num_insts; y++) {
                if (x == orb_idx && y == inst_idx) {
                    bit_string += "1";
                } else {
                    bit_string += "0";
                }
            }
        }
        return bit_string;
    }



}
