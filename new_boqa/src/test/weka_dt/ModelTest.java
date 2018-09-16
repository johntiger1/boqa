package test.weka_dt;
import test.weka_dt.ModelClassifier;
import test.weka_dt.ModelGenerator;
import weka.classifiers.functions.MultilayerPerceptron;


import weka.core.Debug;
import weka.core.Instances;
import weka.core.pmml.jaxbbindings.DecisionTree;
import weka.filters.Filter;
import weka.filters.unsupervised.attribute.Normalize;

//i would likke to import: Decision trees
//import weka.classifiers.trees.Id3;

import weka.classifiers.trees.J48;

import weka.classifiers.*;
public class ModelTest {

        public static final String DATASETPATH = "C:\\Users\\johnp\\Desktop\\git_stuff\\boqa\\weka_resources\\data\\iris.2D.arff";
        public static final String MODElPATH = "C:\\Users\\johnp\\Desktop\\git_stuff\\boqa\\weka_resources\\model\\model.bin";

        public static void main(String[] args) throws Exception {

            ModelGenerator mg = new ModelGenerator();

            Instances dataset = mg.loadDataset(DATASETPATH);

            Filter filter = new Normalize();

            // divide dataset to train dataset 80% and test dataset 20%
            int trainSize = (int) Math.round(dataset.numInstances() * 0.8);
            int testSize = dataset.numInstances() - trainSize;

            dataset.randomize(new Debug.Random(1));// if you comment this line the accuracy of the model will be droped from 96.6% to 80%

            //Normalize dataset
            filter.setInputFormat(dataset);
            Instances datasetnor = Filter.useFilter(dataset, filter);

            Instances traindataset = new Instances(datasetnor, 0, trainSize);
            Instances testdataset = new Instances(datasetnor, trainSize, testSize);

            // build classifier with train dataset
//            MultilayerPerceptron ann = (MultilayerPerceptron ) mg.buildClassifier(traindataset);

            J48 ann = (J48) mg.buildClassifierDT(traindataset);

            // Evaluate classifier with test dataset
            String evalsummary = mg.evaluateModel(ann, traindataset, testdataset);
            System.out.println("Evaluation: " + evalsummary);

            //Save model
            mg.saveModel(ann, MODElPATH);

            //classifiy a single instance
            ModelClassifier cls = new ModelClassifier();
            String classname =cls.classifyDt(Filter.useFilter(cls.createInstance(1.6, 0.2, 0), filter), MODElPATH);
            System.out.println("\n The class name for the instance with petallength = 1.6 and petalwidth =0.2 is  " +classname);

        }


}
