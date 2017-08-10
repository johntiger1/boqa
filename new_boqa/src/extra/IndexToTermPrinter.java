package extra;

import ontologizer.go.Term;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.HashMap;
import java.util.Map;

public class IndexToTermPrinter {


    public IndexToTermPrinter(Object []vertices)
    {

    }

    public static void printMapping(HashMap<Term, Integer> vertices)
    {
        PrintWriter writer = null;
        try {
            writer = new PrintWriter("mapping.txt", "UTF-8");
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (UnsupportedEncodingException e) {
            e.printStackTrace();
        }
        Term t;

        for (Map.Entry e: vertices.entrySet())
        {
            writer.print(e.getValue());
            writer.print("\t");
            writer.print(e.getKey());
            writer.print("\n");
        }
        writer.close();

    }
}
