package test;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

//Executor service is non blocking
public class ESNB {
    public int NUM_THREADS = 8;
    class abc implements Runnable {
        int id;

        public abc(int id) {
            this.id = id;
        }

        public void run() {
            System.out.printf("sda,");
            longRunningCall();
            System.out.printf("I'm done" + id);

        }

        public void longRunningCall()
        {
            try {
                TimeUnit.SECONDS.sleep(10);
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }

    }
    public void createESAndRun() {
        ExecutorService executorService = Executors.newFixedThreadPool(NUM_THREADS);
        for (int i = 0; i < 8; i++) {
            executorService.execute(new abc(i + 1) //this is still not passing a class in


            );
        }
        //we can either add them to a ist and run them all at once
        //or we can do the online running (for loop)
        //how is that different from just a thread thing though...


        executorService.shutdown();
    }

}

