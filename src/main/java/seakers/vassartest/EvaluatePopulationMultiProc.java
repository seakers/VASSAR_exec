package seakers.vassartest;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public class EvaluatePopulationMultiProc {

    public static void main(String[] args) {

        int num_procs = 8;
        int proc_chunks = 5000;
        List<Process> processes = new ArrayList<>();
        try{
            for(int x = 0; x < num_procs; x++){
                System.out.println("--> STARTING PROC: " + x);
                int proc_num = x;
                int pop_lb = proc_chunks * x;
                int pop_ub = proc_chunks * (x + 1);

                List<String> command = new ArrayList<>();
                command.add("java");
                command.add("-cp"); // specify the classpath, you need to fill in the correct classpath
                command.add(System.getProperty("java.class.path"));
                command.add("seakers.vassartest.EvaluatePopulationProc");
                command.add(Integer.toString(proc_num));
                command.add(Integer.toString(pop_lb));
                command.add(Integer.toString(pop_ub));

                ProcessBuilder processBuilder = new ProcessBuilder(command);
                Process process = processBuilder.start();
                processes.add(process);

                try {
                    Thread.sleep(5000); // pauses for 2000 milliseconds, i.e., 2 seconds
                } catch (InterruptedException e) {
                    // Handle the exception
                    Thread.currentThread().interrupt();
                }
            }

            // --> Hook in output
            ExecutorService executorService = Executors.newFixedThreadPool(2);
            executorService.submit(() -> {
                try (BufferedReader reader = new BufferedReader(new InputStreamReader(processes.get(0).getInputStream()))) {
                    String line;
                    while ((line = reader.readLine()) != null) {
                        System.out.println(line);
                    }
                } catch (IOException e) {
                    e.printStackTrace();
                }
            });
            executorService.submit(() -> {
                try (BufferedReader errorReader = new BufferedReader(new InputStreamReader(processes.get(0).getErrorStream()))) {
                    String errorLine;
                    while ((errorLine = errorReader.readLine()) != null) {
                        System.err.println(errorLine);
                    }
                } catch (IOException e) {
                    e.printStackTrace();
                }
            });

            // --> Wait for processes to finish
            for (Process running_proc : processes) {
                running_proc.waitFor();
            }
            executorService.shutdown();
        }
        catch (Exception ex){
            ex.printStackTrace();
        }
    }
}
