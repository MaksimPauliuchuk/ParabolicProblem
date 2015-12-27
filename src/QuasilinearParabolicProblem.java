import com.sun.org.apache.xerces.internal.impl.xpath.regex.Match;

import java.util.ArrayList;
import java.util.Scanner;
import java.util.StringTokenizer;

public class QuasilinearParabolicProblem
{
    double h, tao1, tao2, t, t_optim;
    double stroka1[], stroka2[], stroka[];
    int d, v;
    ArrayList<double[]> matriza;
    static Gui gui;

    QuasilinearParabolicProblem(Gui gui)
    {
        this.gui = gui;
    }

    public void initialization()
    {
        h = (gui.lengthX + 0.0) / gui.N;
        tao1 = Math.min(gui.lengthT, gui.lengthX * (1.0) / gui.N);
        t = 0;
        stroka = new double[gui.N + 1];
        stroka1 = new double[gui.N + 1];
        stroka2 = new double[gui.N + 1];
        matriza = new ArrayList<double[]>();
    }

    public void conditions()
    {
        // y(Xn,0)
        for (int i = 0; i < stroka.length; i++)
        {
            stroka[i] = gui.solveMu1(0, i * h);// <- запихнуть Mu1
        }
        matriza.add(stroka);
    }

    // правило Рунге
    public double RuleRunge()
    {
        double minTao = 1.0;
        while (true)
        {
            for (int i = 0; i < matriza.get(matriza.size() - 1).length; i++)
            {
                stroka[i] = matriza.get(matriza.size() - 1)[i];
            }
            while (true)
            {
                if (tao1 < minTao)
                {
                    minTao = tao1;
                }
                stroka1 = getstroka(stroka, 1, tao1);
                tao2 = tao1 / 2;
                stroka2 = getstroka(stroka, 2, tao2);
                if (norma(stroka1, stroka2) < gui.eps2)
                {
                    matriza.add(stroka1);
                    System.out.print("t = " + t);
                    t += tao1;
                    System.out.println(" tao = " + tao1);
                    tao1 = Math.min(gui.lengthT - t, gui.lengthX * (1.0) / gui.N);
                    break;
                }
                else
                {
                    tao1 = tao2;
                }
            }
            if (t == gui.lengthT)
            {
                break;
            }
        }
        return minTao;
    }

    public double RuleRunge(double tao)
    {
        tao1 = tao / 2;
        System.out.println(tao1);
        while (true)
        {
            for (int i = 0; i < matriza.get(matriza.size() - 1).length; i++)
            {
                stroka[i] = matriza.get(matriza.size() - 1)[i];
            }
            stroka1 = getstroka(stroka, 1, tao1);
            matriza.add(stroka1);
            System.out.println("t = " + t);
            if (t + tao1 > gui.lengthT)
            {
                tao1 = gui.lengthT - t;
            }
            t += tao1;
            if (t == gui.lengthT)
            {
                break;
            }
        }
        return tao1;
    }

    public double norma(double[] stroka1, double[] stroka2)
    {
        double norma = 0.0;
        for (int i = 0; i < stroka1.length; i++)
        {
            norma += Math.pow(stroka1[i] - stroka2[i], 2);
        }
        norma = Math.sqrt(norma);
        return norma;
    }

    public double[] getstroka(double[] stroka, int M, double tao)
    {
        for (int m = 1; m <= M; m++)
        {
            stroka = progonka(t, m, stroka, tao);
        }
        return stroka;
    }

    public void realFunctionAndNeviazka()
    {
        double func;
        System.out.println();
        double nev = 0;
        for (int j = 0; j <= gui.N; j++)
        {
            func = j * h * j * h + gui.lengthT * gui.lengthT + gui.lengthT * gui.lengthT * j * h * j * h;
            nev += Math.pow((Math.abs(func - matriza.get(matriza.size() - 1)[j])), 2);
            // System.out.println(func + " " + matriza.get(matriza.size() - 1)[j]);
        }
        nev = Math.sqrt(nev);
        System.out.printf("%20.15f", nev);
        gui.nevDate.setText(String.format("%20.15f", nev));
        System.out.println();
        System.out.println("d/v = " + (1.0 * d / v));
    }

    public double[] progonka(double t, int m, double[] stroka, double tao)
    {
        double beta = 0.1, matrixD_f[][], vectorDelta_Xn[], vectorF_Xn[], vectorY_M_iter[];
        matrixD_f = new double[gui.N + 1][gui.N + 1];
        vectorF_Xn = new double[gui.N + 1];
        vectorY_M_iter = new double[gui.N + 1];
        System.arraycopy(stroka, 0, vectorY_M_iter, 0, vectorY_M_iter.length);
        while (true)
        {
            matrixD_f[0][0] = matrixD_f[gui.N][gui.N] = 1;
            vectorF_Xn[0] = vectorY_M_iter[0] - gui.solveMu2(tao * m + t, m * h);
            vectorF_Xn[gui.N] = vectorY_M_iter[gui.N] - gui.solveMu3(tao * m + t, gui.lengthX);
            for (int n = 1; n < matrixD_f.length - 1; n++)
            {
                matrixD_f[n][n - 1] = ((gui.solveK2(tao * m + t, n * h, stroka[n])) / (2 * h * h))
                        * (vectorY_M_iter[n + 1] - vectorY_M_iter[n - 1])
                        - gui.solveK1(tao * m + t, n * h, stroka[n]) / (h * h);
                matrixD_f[n][n] = 1.0 / tao + (2 * gui.solveK1(tao * m + t, n * h, stroka[n])) / (h * h);
                matrixD_f[n][n + 1] = -(gui.solveK2(tao * m + t, n * h, stroka[n]) / (2 * h * h))
                        * (vectorY_M_iter[n + 1] - vectorY_M_iter[n - 1])
                        - (gui.solveK1(tao * m + t, n * h, stroka[n])) / (h * h);
                vectorF_Xn[n] = (vectorY_M_iter[n] - stroka[n]) / tao
                        - (gui.solveK2(tao * m + t, n * h, stroka[n]))
                                * Math.pow((vectorY_M_iter[n + 1] - vectorY_M_iter[n - 1]) / (2 * h), 2)
                        - gui.solveK1(tao * m + t, n * h, stroka[n]) / (h * h)
                                * ((vectorY_M_iter[n + 1] - 2 * vectorY_M_iter[n] + vectorY_M_iter[n - 1]))
                        - gui.solveG(tao * m + t, n * h, stroka[n]);
            }
            if (diagonalTest(matrixD_f))
                v++;
            else
                d++;
            double norma_Xn = 0.0;
            for (int n = 0; n < vectorF_Xn.length; n++)
            {
                norma_Xn += Math.pow(vectorF_Xn[n], 2);
                vectorF_Xn[n] *= -beta;
            }
            norma_Xn = Math.sqrt(norma_Xn);
            vectorDelta_Xn = TridiagonalMatrixSolution.Solve(matrixD_f, vectorF_Xn);

            for (int n = 0; n < vectorY_M_iter.length; n++)
            {
                vectorY_M_iter[n] += vectorDelta_Xn[n];
            }
            double norma_XnPlus = 0.0;
            norma_XnPlus += Math.pow(vectorY_M_iter[0] - gui.solveMu2(tao * m + t, m * h), 2);
            norma_XnPlus += Math.pow(vectorY_M_iter[gui.N] - gui.solveMu3(tao * m + t, gui.lengthX), 2);
            for (int n = 1; n < vectorY_M_iter.length - 1; n++)
            {
                norma_XnPlus +=
                        Math.pow(
                                (vectorY_M_iter[n] - stroka[n]) / tao
                                        - (gui.solveK2(tao * m + t, n * h, stroka[n])) * Math
                                                .pow((vectorY_M_iter[n + 1] - vectorY_M_iter[n - 1]) / (2 * h), 2)
                        - gui.solveK1(tao * m + t, n * h, stroka[n]) / (h * h)
                                * ((vectorY_M_iter[n + 1] - 2 * vectorY_M_iter[n] + vectorY_M_iter[n - 1]))
                        - gui.solveG(tao * m + t, n * h, stroka[n]), 2);
            }

            norma_XnPlus = Math.sqrt(norma_XnPlus);
            if (norma_XnPlus <= gui.eps1)
            {
                return vectorY_M_iter;
            }
            else
            {
                beta = Math.min(1.0, beta * norma_Xn / norma_XnPlus);
            }
        }
    }

    public boolean diagonalTest(double[][] matrix)
    {
        boolean check = false;
        double a, b, c;
        for (int i = 0; i < matrix.length; i++)
        {
            if (i == 0)
            {
                a = 0;
                c = matrix[i][i + 1];
            }
            else
            {
                if (i == matrix.length - 1)
                {
                    a = matrix[i][i - 1];
                    c = 0;
                }
                else
                {
                    a = matrix[i][i - 1];
                    c = matrix[i][i + 1];
                }
            }
            b = matrix[i][i];
            double sum = Math.abs(a) + Math.abs(c);
            if (Math.abs(b) < sum)
                return false;
            if (Math.abs(b) > sum)
                check = true;
        }
        return check;
    }
}
