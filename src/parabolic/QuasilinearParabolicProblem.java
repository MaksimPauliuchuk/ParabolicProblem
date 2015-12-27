package parabolic;
public class QuasilinearParabolicProblem
{
    int N, M;
    double lengthX, lengthT;
    double h, tao, epsRunge, epsProgonka;
    double vector[], vectorTao1[], vectorTao2[], t;
    String mu1, mu2, mu3, k1, k2, g;

    public QuasilinearParabolicProblem(Gui gui)
    {
        this.N = gui.N;
        this.M = N; // поменять
        this.lengthX = gui.lengthX;
        this.lengthT = gui.lengthT;
        this.mu1 = gui.textAreas[0].getText();
        this.mu2 = gui.textAreas[1].getText();
        this.mu3 = gui.textAreas[2].getText();
        this.k1 = gui.textAreas[3].getText();
        this.k2 = gui.textAreas[4].getText();
        this.g = gui.textAreas[5].getText();
        this.epsRunge = gui.eps2;
        this.epsProgonka = gui.eps1;
    }

    public void initialization()
    {
        h = lengthX / N;
        tao = lengthT / M;
        vector = new double[N + 1];
        vectorTao1 = new double[N + 1];
        vectorTao2 = new double[N + 1];
        t = 0;
    }

    public void conditions()
    {
        for (int i = 0; i < vector.length; i++)
        {
            vector[i] = solveMu1(0, i * h);
        }
    }

    public void ruleRunge()
    {
        int i = 0;
        double tao1 = tao, tao2;

        vectorTao1 = findAnswerVector(t, 2 ^ i, vector, tao1);
        while (true)
        {
            tao2 = tao1 / 2;
            vectorTao2 = findAnswerVector(t, 2 ^ (i + 1), vector, tao2);
            if (norma(vectorTao1, vectorTao2) < epsRunge)
            {
                vector = vectorTao1.clone();
                t += tao;
                tao = tao1;
                break;
            }
            else
            {
                vectorTao1 = vectorTao2.clone();
                tao1 = tao2;
                i++;
            }
        }
        
        while ( t <= lengthT)
        {
            vector = findAnswerVector(t, 1, vector, tao);
            t+=tao;
        }
        TridiagonalMatrixSolution.Print(vector);
        realFunctionAndNeviazka(vector);
        
    }

    private double[] findAnswerVector(double tBase, int M, double[] stroka, double tao)
    {
        double[] stroka1 = stroka.clone();
        for (int m = 1; m <= M; m++)
        {
            stroka1 = progonka(tBase, m, stroka1, tao);
        }
        return stroka1;
    }

    public double[] progonka(double t, int m, double[] stroka, double tao)
    {
        double beta = 0.1, matrixD_f[][], vectorDelta_Xn[], vectorF_Xn[], vectorY_M_iter[];
        matrixD_f = new double[N + 1][N + 1];
        vectorF_Xn = new double[N + 1];
        vectorY_M_iter = new double[N + 1];
        System.arraycopy(stroka, 0, vectorY_M_iter, 0, vectorY_M_iter.length);
        while (true)
        {
            matrixD_f[0][0] = matrixD_f[N][N] = 1;
            vectorF_Xn[0] = vectorY_M_iter[0] - solveMu2(tao * m + t, m * h);
            vectorF_Xn[N] = vectorY_M_iter[N] - solveMu3(tao * m + t, lengthX);
            for (int n = 1; n < matrixD_f.length - 1; n++)
            {
                matrixD_f[n][n - 1] = ((solveK2(tao * m + t, n * h, stroka[n])) / (2 * h * h))
                        * (vectorY_M_iter[n + 1] - vectorY_M_iter[n - 1])
                        - solveK1(tao * m + t, n * h, stroka[n]) / (h * h);
                matrixD_f[n][n] = 1.0 / tao + (2 * solveK1(tao * m + t, n * h, stroka[n])) / (h * h);
                matrixD_f[n][n + 1] = -(solveK2(tao * m + t, n * h, stroka[n]) / (2 * h * h))
                        * (vectorY_M_iter[n + 1] - vectorY_M_iter[n - 1])
                        - (solveK1(tao * m + t, n * h, stroka[n])) / (h * h);
                vectorF_Xn[n] = (vectorY_M_iter[n] - stroka[n]) / tao
                        - (solveK2(tao * m + t, n * h, stroka[n]))
                                * Math.pow((vectorY_M_iter[n + 1] - vectorY_M_iter[n - 1]) / (2 * h), 2)
                        - solveK1(tao * m + t, n * h, stroka[n]) / (h * h)
                                * ((vectorY_M_iter[n + 1] - 2 * vectorY_M_iter[n] + vectorY_M_iter[n - 1]))
                        - solveG(tao * m + t, n * h, stroka[n]);
            }
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
            norma_XnPlus += Math.pow(vectorY_M_iter[0] - solveMu2(tao * m + t, m * h), 2);
            norma_XnPlus += Math.pow(vectorY_M_iter[N] - solveMu3(tao * m + t, lengthX), 2);
            for (int n = 1; n < vectorY_M_iter.length - 1; n++)
            {
                norma_XnPlus +=
                        Math.pow(
                                (vectorY_M_iter[n] - stroka[n]) / tao
                                        - (solveK2(tao * m + t, n * h, stroka[n])) * Math
                                                .pow((vectorY_M_iter[n + 1] - vectorY_M_iter[n - 1]) / (2 * h), 2)
                        - solveK1(tao * m + t, n * h, stroka[n]) / (h * h)
                                * ((vectorY_M_iter[n + 1] - 2 * vectorY_M_iter[n] + vectorY_M_iter[n - 1]))
                        - solveG(tao * m + t, n * h, stroka[n]), 2);
            }

            norma_XnPlus = Math.sqrt(norma_XnPlus);
            if (norma_XnPlus <= epsProgonka)
            {
                return vectorY_M_iter;
            }
            else
            {
                beta = Math.min(1.0, beta * norma_Xn / norma_XnPlus);
            }
        }
    }

    private double norma(double[] stroka1, double[] stroka2)
    {
        double norma = 0.0;
        for (int i = 0; i < stroka1.length; i++)
        {
            norma += Math.pow(stroka1[i] - stroka2[i], 2);
        }
        norma = Math.sqrt(norma);
        return norma;
    }

    private void realFunctionAndNeviazka(double[] vect)
    {
        double func;
        System.out.println();
        double nev = 0;
        for (int j = 0; j <= N; j++)
        {
            func = j * h * j * h + lengthT * lengthT + lengthT * lengthT * j * h * j * h;
            nev += Math.pow((Math.abs(func - vect[j])), 2);
            System.out.println(func + " " + vect[j]);
        }
        nev = Math.sqrt(nev);
        System.out.printf("%20.15f", nev);
        System.out.println();
    }
    
    private double solveMu1(double t, double x)
    {
        MatchParser matchParser = new MatchParser();
        matchParser.setVariable("t", t);
        matchParser.setVariable("x", x);
        double result = 0.0;
        try
        {
            result = matchParser.Parse(mu1);
        }
        catch (Exception e)
        {
            e.printStackTrace();
        }
        return result;
    }

    private double solveMu2(double t, double x)
    {
        MatchParser matchParser = new MatchParser();
        matchParser.setVariable("t", t);
        matchParser.setVariable("x", x);
        double result = 0.0;
        try
        {
            result = matchParser.Parse(mu2);
        }
        catch (Exception e)
        {
            e.printStackTrace();
        }
        return result;
    }

    private double solveMu3(double t, double x)
    {
        MatchParser matchParser = new MatchParser();
        matchParser.setVariable("t", t);
        matchParser.setVariable("x", x);
        double result = 0.0;
        try
        {
            result = matchParser.Parse(mu3);
        }
        catch (Exception e)
        {
            e.printStackTrace();
        }
        return result;
    }

    private double solveK1(double t, double x, double u)
    {
        MatchParser matchParser = new MatchParser();
        matchParser.setVariable("t", t);
        matchParser.setVariable("x", x);
        matchParser.setVariable("u", u);
        double result = 0.0;
        try
        {
            result = matchParser.Parse(k1);
        }
        catch (Exception e)
        {
            e.printStackTrace();
        }
        return result;

    }

    private double solveK2(double t, double x, double u)
    {
        MatchParser matchParser = new MatchParser();
        matchParser.setVariable("t", t);
        matchParser.setVariable("x", x);
        matchParser.setVariable("u", u);
        double result = 0.0;
        try
        {
            result = matchParser.Parse(k2);
        }
        catch (Exception e)
        {
            e.printStackTrace();
        }
        return result;

    }

    private double solveG(double t, double x, double u)
    {
        MatchParser matchParser = new MatchParser();
        matchParser.setVariable("t", t);
        matchParser.setVariable("x", x);
        matchParser.setVariable("u", u);
        double result = 0.0;
        try
        {
            result = matchParser.Parse(g);
        }
        catch (Exception e)
        {
            e.printStackTrace();
        }
        return result;
    }
}
