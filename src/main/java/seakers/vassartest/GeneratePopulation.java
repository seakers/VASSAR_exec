package seakers.vassartest;

import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import seakers.vassar.problems.Assigning.GigaAssigningParams;
import seakers.vassartest.search.problems.Assigning.DesignBuilder;

import java.io.*;
import java.nio.file.*;
import java.util.*;

public class GeneratePopulation {

    public static void main (String[] args) {
        System.out.println("--> RUNNING GeneratePopulation");
//        generate();
//        stats();
//        getPopulationChunkPruned(0, 10000);
//        generateEvalPop();
//        countEvalPopulationJsonDupes();
//        deleteFilesWithScienceNegativeOne();
        splitFiles();
    }

    public static void deleteFilesWithScienceNegativeOne() {
        String source_dir = "/home/ec2-user/designs2_100k/designs_all";
        try{
            // Initialize Gson's JsonParser
            JsonParser parser = new JsonParser();

            // Getting the list of JSON files in the directory
            DirectoryStream<Path> directoryStream = Files.newDirectoryStream(Paths.get(source_dir), "*.json");
            int del_count = 0;
            for (Path path : directoryStream) {
                try (FileReader reader = new FileReader(path.toFile())) {
                    // Parse JSON file
                    JsonObject jsonObject = (JsonObject) parser.parse(reader);

                    // Check if the "science" key exists and its value is -1
                    if (jsonObject.has("science") && jsonObject.get("science").getAsInt() == -1) {
                        // Delete the file
                        Files.delete(path);
                        System.out.println("--> DELETED FILE");
                        del_count++;
                    }
                }
            }
            System.out.println("--> DELETED: " + del_count);
        }
        catch (Exception ex){
            ex.printStackTrace();
        }

    }

    public static void splitFiles(){
        String source_dir = "/home/ec2-user/designs2_100k/designs_all";
        String train_dir = "/home/ec2-user/designs2_100k/designs_train";
        String val_dir = "/home/ec2-user/designs2_100k/designs_val";
        // Getting the list of JSON files in the source directory
        try{
            DirectoryStream<Path> directoryStream = Files.newDirectoryStream(Paths.get(source_dir), "*.json");
            List<Path> fileList = new ArrayList<>();
            for (Path path : directoryStream) {
                fileList.add(path);
            }

            // Shuffling the file list
            Collections.shuffle(fileList);

            // Calculating the split indices
            int totalFiles = fileList.size();
            int eightPercentIndex = (int) (totalFiles * 0.1);

            // Getting the subsets of files for each destination directory
            List<Path> eightPercentFiles = fileList.subList(0, eightPercentIndex);
            List<Path> ninetyTwoPercentFiles = fileList.subList(eightPercentIndex, totalFiles);

            // Copying files to the 8% directory
            for (Path path : eightPercentFiles) {
                Files.copy(path, Paths.get(val_dir).resolve(path.getFileName()), StandardCopyOption.REPLACE_EXISTING);
            }

            // Copying files to the 92% directory
            for (Path path : ninetyTwoPercentFiles) {
                Files.copy(path, Paths.get(train_dir).resolve(path.getFileName()), StandardCopyOption.REPLACE_EXISTING);
            }

        }
        catch (Exception ex){
            ex.printStackTrace();
        }
    }



    public static void generate(){
        String filePath = "/home/ec2-user/designs2_100k/population.txt";
        int pop_size = 100000;

        // Problem Parameters
        int orekit_threads = 1;
        String resourcesPath = "/home/ec2-user/vassar/giga/VASSAR_resources";
        GigaAssigningParams params = new GigaAssigningParams(resourcesPath, "FUZZY-CASES", "test", "normal", orekit_threads);
        int num_insts = params.getNumInstr();
        int num_orbs = params.getNumOrbits();

        // Create Population
        HashSet<String> uniqueBitStrings = new HashSet<>();
        while (uniqueBitStrings.size() < pop_size) {
            uniqueBitStrings.add(DesignBuilder.randomBitString(num_insts, num_orbs, 40));
        }

        // Write Population
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(filePath))) {
            for (String bitString : uniqueBitStrings) {
                writer.write(bitString);
                writer.newLine();
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    public static void generateEvalPop(){
        String filePath = "/home/ec2-user/designs2_100k/population_eval.txt";
        HashSet<String> eval_designs = getEvalPopulationStringsJson();
        System.out.println("--> EVALED DESIGNS: " + eval_designs.size());

        // Write Population
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(filePath))) {
            for (String bitString : eval_designs) {
                writer.write(bitString);
                writer.newLine();
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }



    public static void stats(){
        List<String> designs = getPopulation();
        Map<Integer, Integer> statistics = computePositiveBitsStatistics(designs);
        for(Map.Entry<Integer, Integer> entry : statistics.entrySet()) {
            System.out.println("Number of 1s: " + entry.getKey() + " - Count: " + entry.getValue());
        }
    }

    public static List<String> getPopulation(){
        String filePath = "/home/ec2-user/designs2_100k/population.txt";
        List<String> designs = new ArrayList<>();
        try{
            designs = Files.readAllLines(Paths.get(filePath));
        }
        catch (Exception ex){
            ex.printStackTrace();
        }
        return designs;
    }

    public static HashSet<String> getEvalPopulation(){
        String filePath = "/home/ec2-user/designs2_100k/population_eval.txt";
        List<String> designs = new ArrayList<>();
        try{
            designs = Files.readAllLines(Paths.get(filePath));
        }
        catch (Exception ex){
            ex.printStackTrace();
        }
        HashSet<String> designs_hash = new HashSet<>(designs);
        return designs_hash;
    }




    public static List<String> getPopulationChunk(int begin, int end){
        String filePath = "/home/ec2-user/designs2_100k/population.txt";
        List<String> designs = new ArrayList<>();
        try{
            designs = Files.readAllLines(Paths.get(filePath));
        }
        catch (Exception ex){
            ex.printStackTrace();
        }

        List<String> designs_chunk = new ArrayList<>();
        for(int x = 0; x < designs.size(); x++){
            if(x >= begin && x < end){
                designs_chunk.add(designs.get(x));
            }
        }
        return designs_chunk;
    }

    public static List<String> getPopulationChunkPruned(int begin, int end){
        String filePath = "/home/ec2-user/designs2_100k/population.txt";
        List<String> designs = new ArrayList<>();
        try{
            designs = Files.readAllLines(Paths.get(filePath));
        }
        catch (Exception ex){
            ex.printStackTrace();
        }

        // Prune Evaluated Designs
        System.out.println("--> BEFORE: " + designs.size());
        HashSet<String> eval_designs = getEvalPopulation();
        designs.removeIf(eval_designs::contains);
        System.out.println("--> AFTER: " + designs.size());

        List<String> designs_chunk = new ArrayList<>();
        for(int x = 0; x < designs.size(); x++){
            if(x >= begin && x < end){
                designs_chunk.add(designs.get(x));
            }
        }
        return designs_chunk;
    }



    public static HashSet<String> getEvalPopulationStringsJson(){
        String path_str = "/home/ec2-user/designs2_100k/designs_all";
        HashSet<String> all_designs = new HashSet<>();

        File directory = new File(path_str);
        File[] files = directory.listFiles((dir, name) -> name.endsWith(".json"));
        if (files != null) {
            JsonParser parser = new JsonParser();
            for (File file : files) {
                try (FileReader reader = new FileReader(file)) {
                    JsonObject jsonObject = parser.parse(reader).getAsJsonObject();
                    all_designs.add(jsonObject.get("design").getAsString());
                } catch (IOException e) {
                    System.err.println("Error reading file: " + file.getName());
                }
            }
        }
        return all_designs;
    }

    public static HashSet<String> countEvalPopulationJsonDupes(){
        String path_str = "/home/ec2-user/designs2_100k/designs_all";
        HashSet<String> all_designs = new HashSet<>();
        File directory = new File(path_str);
        File[] files = directory.listFiles((dir, name) -> name.endsWith(".json"));
        int dup_counter = 0;
        if (files != null) {
            JsonParser parser = new JsonParser();
            for (File file : files) {
                try (FileReader reader = new FileReader(file)) {
                    JsonObject jsonObject = parser.parse(reader).getAsJsonObject();
                    if(all_designs.contains(jsonObject.get("design").getAsString())){
                        dup_counter++;
                    }
                    all_designs.add(jsonObject.get("design").getAsString());
                } catch (IOException e) {
                    System.err.println("Error reading file: " + file.getName());
                }
            }
        }
        System.out.println("--> DUPLICATES: " + dup_counter);
        return all_designs;
    }



    public static Map<Integer, Integer> computePositiveBitsStatistics(List<String> bitStrings) {
        Map<Integer, Integer> statistics = new HashMap<>();
        for (String bitString : bitStrings) {
            int positiveBitsCount = 0;
            for (char bit : bitString.toCharArray()) {
                if (bit == '1') {
                    positiveBitsCount++;
                }
            }
            statistics.put(positiveBitsCount, statistics.getOrDefault(positiveBitsCount, 0) + 1);
        }
        return statistics;
    }



}
