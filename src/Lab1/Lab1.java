package Lab1;

import java.util.Random;
import org.apache.commons.math3.stat.inference.ChiSquareTest;
import org.apache.commons.math3.stat.inference.KolmogorovSmirnovTest;

/**
 * Created by d.kalach on 23.03.2016.
 */
public class Lab1 {
    private static final Random random = new Random();
    private static final long M = (long) Math.pow(2.0, 24.0);
    private static final long A = 25214903917L;
    private static final int K = 64;
    private static final int SAMPLE_SIZE = 100;
    private static int[] V;
    private static long X = 1;


    public static void main(String[] args) {
        double[] expected = new double[SAMPLE_SIZE];
        double[] conqruentDouble = new double[SAMPLE_SIZE];
        double[] macLarenDouble = new double[SAMPLE_SIZE];
        long[] conqruent = new long[SAMPLE_SIZE];
        long[] macLaren = new long[SAMPLE_SIZE];
        int mod = SAMPLE_SIZE;
        initV(mod);
        for (int i = 0; i < SAMPLE_SIZE; i++) {
            expected[i] = generateExpected(mod);
            conqruent[i] = generateCongruential(mod);
            conqruentDouble[i] = conqruent[i];
            macLaren[i] = generateMacLaren(mod);
            macLarenDouble[i] = macLaren[i];
        }
        System.out.println("*************** CHI-SQUARE ***************");
        System.out.println("P-value for conqruent: " + new ChiSquareTest().chiSquareTest(expected, conqruent));
        System.out.println("P-value for MacLaren: " + new ChiSquareTest().chiSquareTest(expected, macLaren));

        System.out.println("*************** Kolmogorov - Smirnov ***************");
        System.out.println("P-value for conqruent: " + new KolmogorovSmirnovTest().kolmogorovSmirnovTest(expected, conqruentDouble));
        System.out.println("P-value for MacLaren: " + new KolmogorovSmirnovTest().kolmogorovSmirnovTest(expected, macLarenDouble));

    }

    public static int generateExpected(int mod) {
        int rand = random.nextInt() % mod;
        return rand < 0 ? -rand + 1 : rand + 1;
    }

    public static int generateCongruential(int mod) {
        X = (X * A) % M;
        long temp = X % mod;
        return (int) temp + 1;
    }

    public static void initV(int mod) {
        V = new int[K];
        for (int i = 0; i < V.length; i++) {
            V[i] = generateExpected(mod);
        }
    }

    public static int generateMacLaren(int mod) {
        if (V == null) {
            throw new RuntimeException("Matrix V haven't been initialized!");
        }
        int x = generateExpected(mod);
        int y = generateCongruential(mod);
        int j = (int) Math.floor((double)K * y / mod);
        int z = V[j];
        V[j] = x;
        return z + 1;
    }


}
