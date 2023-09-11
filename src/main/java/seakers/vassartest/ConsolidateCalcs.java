package seakers.vassartest;

import java.io.IOException;
import java.nio.file.*;

public class ConsolidateCalcs {
    public static void main(String[] args) {
        String path_str = "/home/ec2-user/vassar/giga/VASSAR_resources/orekit/CoverageDatabase";
//        manyToOne(path_str);
        oneToMany(path_str);
    }


    public static void manyToOne(String path_str){
        Path destinationDir = Paths.get(path_str);
        try {
            // Create the destination directory if it doesn't exist
            Files.createDirectories(destinationDir);
            for (int i = 0; i <= 87; i++) {
                Path sourceDir = Paths.get(path_str + "_" + i);
                if (Files.exists(sourceDir)) {
                    try (DirectoryStream<Path> stream = Files.newDirectoryStream(sourceDir)) {
                        for (Path file : stream) {
                            if (Files.isRegularFile(file)) {
                                Files.copy(file, destinationDir.resolve(file.getFileName()), StandardCopyOption.REPLACE_EXISTING);
                            }
                        }
                    }
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }



    public static void oneToMany(String path_str){
        int num_copies = 87;
        Path sourceDir = Paths.get(path_str);
        try {
            if (Files.exists(sourceDir)) {
                try (DirectoryStream<Path> stream = Files.newDirectoryStream(sourceDir)) {
                    for (Path file : stream) {
                        if (Files.isRegularFile(file)) {
                            for (int i = 0; i <= num_copies; i++) {
                                Path destinationDir = Paths.get(path_str + "_" + i);
                                if(!Files.exists(destinationDir)){
                                    Files.createDirectories(destinationDir);
                                }
                                Files.copy(file, destinationDir.resolve(file.getFileName()), StandardCopyOption.REPLACE_EXISTING);
                            }
                        }
                    }
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }




}
