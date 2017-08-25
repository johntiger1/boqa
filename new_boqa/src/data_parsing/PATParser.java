package data_parsing;

import com.github.phenomics.ontolib.formats.hpo.HpoDiseaseAnnotation;
import com.github.phenomics.ontolib.io.base.TermAnnotationParserException;
import com.github.phenomics.ontolib.io.obo.hpo.HpoDiseaseAnnotationParser;

import java.io.File;
import java.io.IOException;

public class PATParser {
    public static void main(String[]args)
    {

        File inputFile = new File("phenotype_annotation.tab");
        
        try {
            HpoDiseaseAnnotationParser parser = new HpoDiseaseAnnotationParser(inputFile);
            while (parser.hasNext()) {
                HpoDiseaseAnnotation anno = parser.next();
                System.out.println(anno);

                //we will probably get some dictionary that allows us
                //to map between the String and a numeric value
//                anno.getFrequencyModifier();
                // work with anno
            }
        } catch (IOException e) {
        System.err.println("Problem reading from file.");
    } catch (TermAnnotationParserException e) {
        System.err.println("Problem parsing file.");
    }
    }

}
