package data_parsing.Adapter;

import com.github.phenomics.ontolib.ontology.data.TermId;
import ontologizer.go.TermID;

public class TermIDUpgrade {

    //helper class that interconverts between TermId and TermID and vice versa

    public static TermID IdToID(TermId termId)
    {
        return new TermID(termId.getIdWithPrefix());
    }
}
