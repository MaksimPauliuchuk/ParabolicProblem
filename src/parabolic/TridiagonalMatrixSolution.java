package parabolic;
/**
 * Created by ������ on 04.10.2015.
 */
public class TridiagonalMatrixSolution
{
    static double[] Solve(double[][] matrix, double[] f)
    {
        double[][] AdditionalValues = new double[2][matrix.length];
        AdditionalValues[0][0] = -matrix[0][1] / matrix[0][0];
        AdditionalValues[1][0] = f[0] / matrix[0][0];
        for (int i = 1; i < matrix.length; i++)
        {
            if (!(i == matrix.length - 1))
            {
                AdditionalValues[0][i] =
                        -matrix[i][i + 1] / (matrix[i][i] + matrix[i][i - 1] * AdditionalValues[0][i - 1]);
            }
            AdditionalValues[1][i] = (-matrix[i][i - 1] * AdditionalValues[1][i - 1] + f[i])
                    / (matrix[i][i] + matrix[i][i - 1] * AdditionalValues[0][i - 1]);
        }

        double[] x = new double[matrix.length];
        x[x.length - 1] = AdditionalValues[1][x.length - 1];
        for (int i = x.length - 2; i >= 0; i--)
        {
            x[i] = AdditionalValues[0][i] * x[i + 1] + AdditionalValues[1][i];
        }

        return x;
    }

    static double[] Solve(double[] left, double[] center, double[] right, double[] f)
    {
        double[][] AdditionalValues = new double[2][center.length];
        AdditionalValues[0][0] = -right[0] / center[0];
        AdditionalValues[1][0] = f[0] / center[0];
        for (int i = 1; i < center.length; i++)
        {
            if (!(i == center.length - 1))
            {
                AdditionalValues[0][i] =
                        -right[i] / (center[i] + left[i - 1] * AdditionalValues[0][i - 1]);
            }
            AdditionalValues[1][i] = (-left[i - 1] * AdditionalValues[1][i - 1] + f[i])
                    / (center[i] + left[i - 1] * AdditionalValues[0][i - 1]);
        }

        double[] x = new double[center.length];
        x[x.length - 1] = AdditionalValues[1][x.length - 1];
        for (int i = x.length - 2; i >= 0; i--)
        {
            x[i] = AdditionalValues[0][i] * x[i + 1] + AdditionalValues[1][i];
        }

        return x;
    }
    
    static void Print(double[][] matrix)
    {
        System.out.println("Matrix demension is " + matrix.length + "x" + matrix[0].length + ".");
        for (int i = 0; i < matrix.length; i++)
        {
            for (int j = 0; j < matrix[i].length; j++)
            {
                System.out.printf("%9.3f ", matrix[i][j]);
            }
            System.out.println();
        }
        System.out.println();
    }

    static void Print(double[] f)
    {
        System.out.println("Demension is " + f.length + ".");
        for (int i = 0; i < f.length; i++)
        {
            System.out.printf("%9.13f ", f[i]);
        }
        System.out.println();
        System.out.println();
    }

    static void Print(double[][] matrix, double[] f)
    {
        System.out.println("Matrix demension is " + matrix.length + "x" + matrix[0].length + ".");
        System.out.println("Free row demension is " + f.length + ".");
        for (int i = 0; i < matrix.length; i++)
        {
            for (int j = 0; j < matrix[i].length; j++)
            {
                System.out.printf("%6.2f ", matrix[i][j]);
            }
            System.out.print(" | ");
            System.out.printf("%7.2f ", f[i]);
            System.out.println();
        }
        System.out.println();
    }
}
