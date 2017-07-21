package calculation_john_pack;

import java.util.Arrays;

public class testIncrement {


    public static void main(String []args)
    {
        int size = 10000;
        int i=0;
        int j = 34;
        int [] array = new int[size];
        while (j < 10000)
        {
            array[i++] = j++;

        }

        System.out.println(Arrays.toString(array));

        int[]other = new int[size];

        j=34;
        i=0;
        while (j < 10000)
        {
            other[i] = j;
            i++;
            j++;

        }

        System.out.println(Arrays.toString(other));

//        for (int i =0; i< 10000; i++)
//        {
//
//
//        }
    }
}

