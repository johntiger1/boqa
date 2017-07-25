package test;


import java.util.Arrays;
import java.util.concurrent.TimeUnit;

//this will be a success if we can show threads doing different things (i.e. using their different args to start at
//different places in an array)
public class MultiThreadEx {



    //different threads will add up all the elements

    public static  void populateArray(int [] arr)
    {
        int j = 1;
        for (int i = 0; i< arr.length; i++)
        {
            //every 1000 elts, multiply by a 10 factor
            arr[i] = i*j;
//            if (i%1000 == 0)
//            {
//                j*=2;
//            }
        }


    }



    //taking it slow and just explicitly declaring the runnable

    public static class MySummmer implements Runnable
    {
        int id;
        //if you are the last thread, consume until the end
        boolean last;
        int startInd;
        int endInd;
        int [] arr;
        public MySummmer(int x, int startInd, int endInd)
        {
            id = x;
            this.startInd = startInd;
            this.endInd = endInd;
        }
        public void anymethod()
        {
            //make the class, then call.run on it?
        }

        @Override
        public void run() {
            sumArray(arr);


        }
        public void sumArray(int [] arr)
        {
            int sum = 0;
            for (int i = startInd; i< endInd; i++)
            {
                try{Thread.sleep(5);}
                catch (InterruptedException ie)
                {

                }
                sum += arr[i];

            }

            System.out.printf("sum from %d to %d is %d\n" ,this.startInd,this.endInd, sum);


        }
        //specify whether this is the last thread
        public void setLast(boolean lst)
        {
            this.last = lst;

        }

        public void setData(int []data)
        {
            this.arr = data;

        }
    }

    public static void main(String[]args)
    {
        final int SIZE = 1000;

        int []arr = new int[SIZE];
        populateArray(arr);
        MySummmer temp;
        long start = System.nanoTime();
        final int NUM_THREADS = 8;


        //System.out.println(Arrays.toString(arr));


        Thread  [] threads = new Thread[NUM_THREADS];

        int start_ind;
        int end;
        for (int i = 0; i < 8; i++){

            start_ind = i*(SIZE/NUM_THREADS);
            end = (i+1)*(SIZE/NUM_THREADS);

            temp = new MySummmer(i, start_ind ,end );

            temp.setData(arr);
            threads[i] = new Thread(temp); //|point to the same thing"
            threads[i].start();

        }int i;
        for(i = 0; i < threads.length; i++) {
            try{
                threads[i].join();}
            catch (InterruptedException ie){
                System.out.println("oop:");
            }
        }
        System.out.println("took" + (System.nanoTime()-start));
        temp = new MySummmer(10,0,arr.length);
        temp.setData(arr);
        temp.sumArray(arr);
        System.out.println("took" + (System.nanoTime()-start));
       //MultiThreadEx me = new MultiThreadEx();

       // MySummmer m1 = new MySummmer()
        //Runnable r = new MySummmer(3,3,3);
        //MultiThreadEx asd= new MultiThreadEx();

//        ExecutorService executorService = Executors.newFixedThreadPool(10);
//
//        executorService.execute(new Runnable() {
//            public void run() {
//                System.out.println("Asynchronous task");
//            }
//        });
//
//        executorService.shutdown();
    }

}

