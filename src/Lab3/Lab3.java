package Lab3;

import java.util.Random;
import org.apache.commons.math3.distribution.WeibullDistribution;
import org.apache.commons.math3.stat.inference.KolmogorovSmirnovTest;

/**
 * Created by d.kalach on 23.03.2016.
 */
public class Lab3 {
    private static final int SAMPLE_SIZE = 100;
    private static final double SHAPE = 1;
    private static final double SCALE = 1;
    private static final Random random = new Random();

    public static void main(String[] args) {
        double[] sample = new double[SAMPLE_SIZE];
        for (int i = 0; i < SAMPLE_SIZE; i++) {
            sample[i] = nextDouble();
        }
        System.out.println("P-value for Weibull: " +
            new KolmogorovSmirnovTest().kolmogorovSmirnovTest(new WeibullDistribution(SHAPE, SCALE), sample));

    }

    public static double nextDouble() {
        double uniform = random.nextDouble();
        return Math.pow( - 1.0 / SCALE * Math.log(uniform), 1.0 / SHAPE);
    }
}
