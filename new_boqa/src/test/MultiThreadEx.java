package test;


import java.util.Arrays;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

//this will be a success if we can show threads doing different things (i.e. using their different args to start at
//different places in an array)
public class MultiThreadEx {


    //different threads will add up all the elements

    public static void populateArray(int[] arr) {
        int j = 1;
        for (int i = 0; i < arr.length; i++) {
            //every 1000 elts, multiply by a 10 factor
            arr[i] = i * j;
//            if (i%1000 == 0)
//            {
//                j*=2;
//            }
        }


    }


    //taking it slow and just explicitly declaring the runnable

    public static class MySummmer implements Runnable {
        int id;
        //if you are the last thread, consume until the end
        boolean last;
        int startInd;
        int endInd;
        int[] arr;

        public MySummmer(int x, int startInd, int endInd) {
            id = x;
            this.startInd = startInd;
            this.endInd = endInd;
        }

        public void anymethod() {
            //make the class, then call.run on it?
        }

        @Override
        public void run() {
            sumArray(arr);


        }

        public void sumArray(int[] arr) {
            int sum = 0;
            for (int i = startInd; i < endInd; i++) {
                try {
                    Thread.sleep(5);
                } catch (InterruptedException ie) {

                }
                sum += arr[i];

            }

            System.out.printf("sum from %d to %d is %d\n", this.startInd, this.endInd, sum);


        }

        //specify whether this is the last thread
        public void setLast(boolean lst) {
            this.last = lst;

        }

        public void setData(int[] data) {
            this.arr = data;

        }
    }

    public static void main(String[] args) {

        final int NUM_THREADS = 8;
        //MultiThreadEx me = new MultiThreadEx();

        // MySummmer m1 = new MySummmer()
        //Runnable r = new MySummmer(3,3,3);
        //MultiThreadEx asd= new MultiThreadEx();
        class abc implements Runnable {
            int id;

            public abc(int id) {
                this.id = id;
            }

            public void run() {
                System.out.printf("sda,");
            }

        }
        ExecutorService executorService = Executors.newFixedThreadPool(NUM_THREADS);
        for (int i = 0; i < 8; i++) {
            executorService.execute(new abc(i) //this is still not passing a class in


            );
        }
        //we can either add them to a ist and run them all at once
        //or we can do the online running (for loop)
        //how is that different from just a thread thing though...


        executorService.execute(new MySummmer(3, 3, 3));

        executorService.shutdown();


    }
}

