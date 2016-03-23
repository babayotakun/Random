package Lab2;

import java.util.Random;

/**
 * Created by d.kalach on 23.03.2016.
 */
public class MyHyperGeometric {
    private Random random = new Random();

    public int getNext(int N, int M, int n) {
        return get(init(N, M), n);
    }

    private int[] init(int N, int M) {
        int[] array = new int[N];
        for (int i = 0; i < M; i++) {
            int index = random.nextInt() % N;
            index = index < 0 ? -index : index;
            if (array[index] == 1) {
                i--;
            } else {
                array[index] = 1;
            }
        }
        return array;
    }

    private int get(int[] array, int n) {
        int sum = 0;
        for (int i = 0; i < n; i++) {
            int index = random.nextInt() % array.length;
            index = index < 0 ? -index : index;
            if (array[index] > -1) {
                sum += array[index];
                array[index] = -1;
            } else {
                i--;
            }
        }
        return sum;
    }
}
