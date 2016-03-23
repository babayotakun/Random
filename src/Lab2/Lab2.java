package Lab2;

import org.apache.commons.math3.distribution.HypergeometricDistribution;
import org.apache.commons.math3.stat.inference.KolmogorovSmirnovTest;

/**
 * Created by d.kalach on 23.03.2016.
 */
public class Lab2 {
    private static final int N = 30;
    private static final int M = 20;
    private static final int n = 20;
    private static final int SAMPLE_SIZE = 100;


    public static void main(String[] args) {
        double[] hyper = new double[SAMPLE_SIZE];
        double[] sample = new double[SAMPLE_SIZE];
        HypergeometricDistribution distribution = new HypergeometricDistribution(N, M, n);
        for (int i = 0; i < SAMPLE_SIZE; i++) {
            hyper[i] = HyperGeometric.staticNextInt(N, M, n);
            sample[i] = distribution.sample();
        }
        System.out.println("P-value for hyper geometric: " + new KolmogorovSmirnovTest().kolmogorovSmirnovTest(sample, hyper));
    }

}
