package test;


import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

//this will be a success if we can show threads doing different things (i.e. using their different args to start at
//different places in an array)
public class multiThreadEx {



    //different threads will add up all the elements

    public void populateArray(int [] arr)
    {
        int j = 1;
        for (int i = 0; i< arr.length; i++)
        {
            //every 1000 elts, multiply by a 10 factor
            arr[i] = i*j;
            if (i%1000 == 0)
            {
                j*=10;
            }
        }


    }

    //taking it slow and just explicitly declaring the runnable

    public class mySummmer implements Runnable
    {
        int id;
        //if you are the last thread, consume until the end
        boolean last;
        int [] arr;
        public mySummmer(int x)
        {
            id = x;
        }
        public void anymethod()
        {
            //make the class, then call.run on it?
        }

        @Override
        public void run() {
            for (int i =id*; i < )


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

        ExecutorService executorService = Executors.newFixedThreadPool(10);

        executorService.execute(new Runnable() {
            public void run() {
                System.out.println("Asynchronous task");
            }
        });

        executorService.shutdown();
    }

}


