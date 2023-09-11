package seakers.vassartest.search.problems.Assigning;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Random;

public class DesignBuilder {

    public static Random rand = new Random();


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
