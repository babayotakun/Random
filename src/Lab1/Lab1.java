package Lab1;

import java.util.Arrays;
import java.util.Random;
import org.apache.commons.math3.distribution.UniformRealDistribution;
import org.apache.commons.math3.stat.inference.ChiSquareTest;
import org.apache.commons.math3.stat.inference.KolmogorovSmirnovTest;

/**
 * Created by d.kalach on 23.03.2016.
 */
public class Lab1 {
    private static final Random random = new Random();
    private static final long M = (long) Math.pow(2.0, 24.0);
    private static final long A = 1997;
    private static final int K = 64;
    private static final int SAMPLE_SIZE = 100000;
    private static int[] V;
    private static long X = 1;


    public static void main(String[] args) {
        double[] conqruentDouble = new double[SAMPLE_SIZE];
        double[] macLarenDouble = new double[SAMPLE_SIZE];
        long[] conqruent = new long[SAMPLE_SIZE];
        long[] macLaren = new long[SAMPLE_SIZE];
        int mod = (int) A;
        initV(mod);
        for (int i = 0; i < SAMPLE_SIZE; i++) {
            conqruent[i] = generateCongruential(mod);
            conqruentDouble[i] = conqruent[i];
            macLaren[i] = generateMacLaren(mod);
            macLarenDouble[i] = macLaren[i];
        }
        System.out.println("*************** Chi-Square ***************");
        System.out.println("P-value for conqruent: " +
            new ChiSquareTest().chiSquareTest(theoreticalFrequencies(SAMPLE_SIZE, mod), frequencies(conqruent, mod)));
        System.out.println("P-value for MacLaren: " +
            new ChiSquareTest().chiSquareTest(theoreticalFrequencies(SAMPLE_SIZE, mod), frequencies(macLaren, mod)));
        System.out.println();
        System.out.println("*************** Kolmogorov - Smirnov *****");
        System.out.println("P-value for conqruent: " +
            new KolmogorovSmirnovTest().kolmogorovSmirnovTest(new UniformRealDistribution(1, mod), conqruentDouble));
        System.out.println("P-value for MacLaren: " +
            new KolmogorovSmirnovTest().kolmogorovSmirnovTest(new UniformRealDistribution(1, mod), macLarenDouble));

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
        V = new int[K + 1];
        for (int i = 0; i < V.length; i++) {
            V[i] = generateExpected(mod);
        }
    }

    private static long[] frequencies(long[] data, int mod) {
        long[] frequencies = new long[mod + 2];
        for (int j = 0; j < data.length; j++) {
            frequencies[(int) data[j]]++;
        }
        return frequencies;
    }

    private static double[] theoreticalFrequencies(int size, int mod) {
        double frequency = (double) size / mod;
        double[] frequencies = new double[mod + 2];
        Arrays.fill(frequencies, frequency);
        return frequencies;
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
