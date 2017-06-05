package sonumina.boqa.calculation;

import ontologizer.go.TermID;

/**
 * Created by johnp on 2017-06-02.
 */
public interface PhenotypeSelector {


    TermID getBestPhenotype(Object o);
}
