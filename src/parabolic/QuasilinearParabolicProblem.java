package parabolic;

public class QuasilinearParabolicProblem
{
    private int N, M, numberBeta;
    private double lengthX, lengthT;
    private double h, tao, epsRunge, epsProgonka;
    double vector[], vectorTao1[], vectorTao2[], t;
    private String mu1, mu2, mu3, k1, k2, g, real;
    private Gui gui;

    public QuasilinearParabolicProblem(Gui gui)
    {
        this.gui = gui;
        this.N = gui.N;
        this.M = gui.M;
        this.numberBeta = Integer.parseInt(gui.textBeta.getText());
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
        this.real = gui.real.getText();
    }

    public void initialization()
    {
        h = lengthX / N;
        tao = lengthT / M;
        // tao = 0.001; // подобрал хорошее тау
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

    public double ruleRunge()
    {
        int i = 0;
        double tao1 = tao, tao2;
        long time = System.currentTimeMillis();
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
        System.out.println("Tao: \t" + tao);
        System.out.println("Time to find tao optim: " + (System.currentTimeMillis() - time));
        long max = 0;
        long itter = 0;
        long Maxitter = (long) ((lengthT - t) / tao) + 1;
        while (t <= lengthT)
        {

            time = System.currentTimeMillis();
            vector = findAnswerVector(t, 1, vector, tao);
            t += tao;
            itter++;
            if (itter == (long) (0.1 * Maxitter))
            {
                System.out.print("10% ");
            }
            else if (itter == (long) (0.2 * Maxitter))
            {
                System.out.print("20% ");
            }
            else if (itter == (long) (0.3 * Maxitter))
            {
                System.out.print("30% ");
            }
            else if (itter == (long) (0.4 * Maxitter))
            {
                System.out.print("40% ");
            }
            else if (itter == (long) (0.5 * Maxitter))
            {
                System.out.print("50% ");
            }
            else if (itter == (long) (0.6 * Maxitter))
            {
                System.out.print("60% ");
            }
            else if (itter == (long) (0.7 * Maxitter))
            {
                System.out.print("70% ");
            }
            else if (itter == (long) (0.8 * Maxitter))
            {
                System.out.print("80% ");
            }
            else if (itter == (long) (0.9 * Maxitter))
            {
                System.out.print("90%");
            }
            time = System.currentTimeMillis() - time;
            if (max < time)
                max = time;
        }
        System.out.println("\nGoing to end: " + max);
        System.out.println("Count of itterations = " + itter);
        TridiagonalMatrixSolution.Print(vector);
        vectorTao1 = vector.clone();
        realFunctionAndNeviazka(vector);
        return tao1;
    }

    public void findAnsverFromTao(double tao)
    {
        long time = System.currentTimeMillis();
        t = 0;
        System.out.println("Tao/2: \t" + tao);
        long max = 0;
        long itter = 0;
        long Maxitter = (long) ((lengthT - t) / tao) + 1;
        while (t <= lengthT)
        {

            time = System.currentTimeMillis();
            vector = findAnswerVector(t, 1, vector, tao);
            t += tao;
            itter++;
            if (itter == (long) (0.1 * Maxitter))
            {
                System.out.print("10% ");
            }
            else if (itter == (long) (0.2 * Maxitter))
            {
                System.out.print("20% ");
            }
            else if (itter == (long) (0.3 * Maxitter))
            {
                System.out.print("30% ");
            }
            else if (itter == (long) (0.4 * Maxitter))
            {
                System.out.print("40% ");
            }
            else if (itter == (long) (0.5 * Maxitter))
            {
                System.out.print("50% ");
            }
            else if (itter == (long) (0.6 * Maxitter))
            {
                System.out.print("60% ");
            }
            else if (itter == (long) (0.7 * Maxitter))
            {
                System.out.print("70% ");
            }
            else if (itter == (long) (0.8 * Maxitter))
            {
                System.out.print("80% ");
            }
            else if (itter == (long) (0.9 * Maxitter))
            {
                System.out.print("90%");
            }
            time = System.currentTimeMillis() - time;
            if (max < time)
                max = time;
        }
        System.out.println("\nGoing to end: " + max);
        System.out.println("Count of itterations = " + itter);
        TridiagonalMatrixSolution.Print(vector);
        vectorTao2 = vector.clone();
        Double a = realFunctionAndNeviazka(vector);
        gui.nevDate.setText(a.toString());
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
        double beta, gamma = 0.01, betaminus, vectorDelta_Xn[], vectorF_Xn[], vectorY_M_iter[];
        double[] left = new double[N], center = new double[N + 1], right = new double[N];
        vectorF_Xn = new double[N + 1];
        vectorY_M_iter = stroka.clone();

        switch (numberBeta)
        {
            case 1:
            {
                beta = 1;
                break;
            }
            case 2:
            {
                beta = 0.1;
                break;
            }
            case 3:
            {
                beta = 0.1;
                gamma = beta * beta;
                break;
            }
            default:
                beta = 0.1;
                break;
        }
        while (true)
        {
            center[0] = center[N] = 1;
            vectorF_Xn[0] = vectorY_M_iter[0] - solveMu2(tao * m + t, m * h);
            vectorF_Xn[N] = vectorY_M_iter[N] - solveMu3(tao * m + t, lengthX);
            for (int n = 1; n < center.length - 1; n++)
            {
                left[n - 1] = ((solveK2(tao * m + t, n * h, stroka[n])) / (2 * h * h))
                        * (vectorY_M_iter[n + 1] - vectorY_M_iter[n - 1])
                        - solveK1(tao * m + t, n * h, stroka[n]) / (h * h);
                center[n] = 1.0 / tao + (2 * solveK1(tao * m + t, n * h, stroka[n])) / (h * h);
                right[n] = -(solveK2(tao * m + t, n * h, stroka[n]) / (2 * h * h))
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
            vectorDelta_Xn = TridiagonalMatrixSolution.Solve(left, center, right, vectorF_Xn);

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
                switch (numberBeta)
                {
                    case 1:
                    {
                        beta = 1;
                        break;
                    }
                    case 2:
                    {
                        beta = Math.min(1.0, beta * norma_Xn / norma_XnPlus);
                        break;
                    }
                    case 3:
                    {
                        betaminus = beta;
                        beta = Math.min(1.0, (gamma * norma_Xn) / (norma_XnPlus * beta));
                        gamma = gamma * (norma_Xn / norma_XnPlus) * (beta / betaminus);
                        break;
                    }
                    default:
                        beta = Math.min(1.0, beta * norma_Xn / norma_XnPlus);
                        break;
                }
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

    private double realFunctionAndNeviazka(double[] vect)
    {
        double func;
        double nev = 0;
        for (int j = 0; j <= N; j++)
        {
            func = j * h * j * h + lengthT * lengthT + lengthT * lengthT * j * h * j * h;
            // func = j * h * j * h + lengthT;
            nev += Math.pow((Math.abs(solveReal(lengthT, j * h) - vect[j])), 2);
            // System.out.println(func + " " + vect[j]);
        }
        nev = Math.sqrt(nev);
        System.out.printf("Residual with concrete solve:\n%20.15f", nev);
        System.out.println();
        return nev;
    }

    void findResidual()
    {
        double nev = 0;
        for (int j = 0; j <= N; j++)
        {
            nev += Math.pow((Math.abs(vectorTao1[j] - vectorTao2[j])), 2);
            // System.out.println(func + " " + vect[j]);
        }
        nev = Math.sqrt(nev);
        Double n = nev;
        gui.nev2Date.setText(n.toString());
        System.out.printf("\nResidual: %20.15f", nev);
        System.out.println();
        System.out.println();
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

    private double solveReal(double t, double x)
    {
        MatchParser matchParser = new MatchParser();
        matchParser.setVariable("t", t);
        matchParser.setVariable("x", x);
        double result = 0.0;
        try
        {
            result = matchParser.Parse(real);
        }
        catch (Exception e)
        {
            e.printStackTrace();
        }
        return result;

    }
}
