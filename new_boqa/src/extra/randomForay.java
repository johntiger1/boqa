package extra;

import calculation_john_pack.ReducedBoqa;
import ontologizer.types.ByteString;

import java.util.Random;

public class randomForay {

    public static void main(String[]args)
    {
        Random r = new Random(23);
        System.out.println("do some stuff");
        ReducedBoqa rb = new ReducedBoqa();
        System.out.println(12*124);

        for (int i = 0; i < 100; i++)
        System.out.println(r.nextInt(12));

        ByteString [] bs = new ByteString[11];
        bs[1] = new ByteString("aa");

    }
}

