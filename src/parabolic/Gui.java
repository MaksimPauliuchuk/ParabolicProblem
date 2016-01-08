package parabolic;

import javax.swing.*;

import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYSplineRenderer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import java.awt.*;
import java.awt.event.*;

public class Gui extends JFrame
{
    int width = 800, height = width * 9 / 16;
    int countLabels = 12, countTextAreas = 12;
    JLabel logo = new JLabel("Место для вашей рекламы");
    JLabel[] labels = new JLabel[countLabels];
    JLabel nevTitle = new JLabel("Невязка:");
    JLabel nevDate = new JLabel();
    JLabel nev2Title = new JLabel("Невязка Рунге:");
    JLabel nev2Date = new JLabel();
    JLabel realFunc = new JLabel("Функция");
    JLabel[] labelBeta = new JLabel[3];
    JTextArea textBeta = new JTextArea();
    JTextArea[] textAreas = new JTextArea[countTextAreas];
    JTextArea real = new JTextArea();
    ImageIcon[] pictures = new ImageIcon[countLabels];
    JButton start = new JButton("Пуск");
    double lengthX, lengthT, eps1, eps2;
    int N, M;

    public Gui()
    {
        super("<^_^>");
        // --------------
        setLayout(null);
        setSize(width, height);
        // --------------------
        createPict();
        createAllLabels();
        createButton();
        createNev();
        createBackground();
        // -------Создание формы
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        setResizable(false);
        setLocationRelativeTo(null);
        setVisible(true);
        // ---------------------
    }

    void createAllLabels()
    {
        logo = new JLabel(new ImageIcon("pict/logo.png"));
        logo.setBounds((width - 292) / 2, 20, 292, 75);
        logo.setBorder(BorderFactory.createLineBorder(Color.gray, 1));
        add(logo);
        int otstupHoriz = width * 5 / 100;
        int otstupVertic = height * 5 / 100;

        for (int i = 0; i < countLabels; i++)
        {
            labels[i] = new JLabel(pictures[i]);
            textAreas[i] = new JTextArea();
            textAreas[i].setBackground(Color.lightGray);
        }
        
        for(int i = 0; i < 3; i++)
        {
            labelBeta[i] = new JLabel();
        }
        labelBeta[0].setText("1. = Ньютона");
        labelBeta[1].setText("2. = Пузынина");
        labelBeta[2].setText("3. = Модифицированной невязки");
        textBeta.setText("2");
        textBeta.setBackground(Color.lightGray);
        real.setBackground(Color.lightGray);

        textAreas[0].setText("x*x");
        textAreas[1].setText("t*t");
        textAreas[2].setText("x*x+t*t+x*x*t*t");
        textAreas[3].setText("u*u*x+u*t");
        textAreas[4].setText("2*u*x+t");
        textAreas[5].setText(
                "2*t-2*t*u-2*t*t*t*u-2*u*u*x-2*t*t*u*u*x-2*t*x*x-8*t*t*t*x*x-4*t*t*t*t*t*x*x-8*u*x*x*x-16*t*t*u*x*x*x-8*t*t*t*t*u*x*x*x");
        real.setText("x*x+t*t+x*x*t*t");
        /*
         * textAreas[0].setText("x*x"); textAreas[1].setText("t"); textAreas[2].setText("x*x+t");
         * textAreas[3].setText("u*x"); textAreas[4].setText("x"); textAreas[5].setText("1-4*x*x*x-2*u*x");
         */

        textAreas[6].setText("1");
        textAreas[7].setText("1");
        textAreas[8].setText("1e-4");
        textAreas[9].setText("1e-4");
        textAreas[10].setText("10");
        textAreas[11].setText("10");

        int nowWidth = otstupHoriz, nowHeight = 75 + otstupVertic;

        for (int i = 0; i < 3; i++)
        {
            labels[i].setBounds(nowWidth, otstupVertic + nowHeight, pictures[i].getIconWidth(),
                    pictures[i].getIconHeight());
            nowWidth += pictures[i].getIconWidth();
            textAreas[i].setBounds(nowWidth, nowHeight + (pictures[i].getIconHeight() + 20) / 2, 90, 16);
            nowWidth += 90 + otstupHoriz;
        }
        nowHeight += 2 * otstupVertic;
        nowWidth = otstupHoriz;
        for (int i = 3; i < 5; i++)
        {
            labels[i].setBounds(nowWidth, otstupVertic + nowHeight, pictures[i].getIconWidth(),
                    pictures[i].getIconHeight());
            nowWidth += pictures[i].getIconWidth();
            textAreas[i].setBounds(nowWidth, nowHeight + (pictures[i].getIconHeight() + 20) / 2, 50, 16);
            nowWidth += 50 + otstupHoriz;
        }
        nowHeight += 2 * otstupVertic;
        nowWidth = otstupHoriz;

        int i = 5;
        labels[i].setBounds(nowWidth, otstupVertic + nowHeight, pictures[i].getIconWidth(),
                pictures[i].getIconHeight());
        nowWidth += pictures[i].getIconWidth();
        textAreas[i].setBounds(nowWidth, nowHeight + (pictures[i].getIconHeight() + 20) / 2, 300, 16);

        nowWidth += 335;
        nowHeight -= 2 * otstupVertic;

        for (i = 6; i < 7; i++)
        {
            labels[i].setBounds(nowWidth, otstupVertic + nowHeight, pictures[i].getIconWidth(),
                    pictures[i].getIconHeight());
            labels[i + 1].setBounds(nowWidth, 3 * otstupVertic + nowHeight, pictures[i + 1].getIconWidth(),
                    pictures[i + 1].getIconHeight());
            nowWidth += Math.max(pictures[i].getIconWidth(), pictures[i + 1].getIconWidth());
            textAreas[i].setBounds(nowWidth, nowHeight + (pictures[i].getIconHeight() + 20) / 2, 30, 16);
            textAreas[i + 1].setBounds(nowWidth,
                    2 * otstupVertic + nowHeight + (pictures[i + 1].getIconHeight() + 20) / 2, 30, 16);
        }

        nowWidth = otstupHoriz;
        nowHeight += 3 * otstupVertic + 10;
        i = 10;
        labels[i].setBounds(nowWidth, otstupVertic + nowHeight, pictures[i].getIconWidth(),
                pictures[i].getIconHeight());
        nowWidth += pictures[i].getIconWidth();
        textAreas[i].setBounds(nowWidth, nowHeight + (pictures[i].getIconHeight() + 20) / 2, 30, 16);

        i = 11;
        labels[i].setBounds(nowWidth - pictures[i].getIconWidth(), otstupVertic + 30 + nowHeight,
                pictures[i].getIconWidth(), pictures[i].getIconHeight());
        textAreas[i].setBounds(nowWidth, nowHeight + 30 + (pictures[i].getIconHeight() + 20) / 2, 30, 16);

        nowWidth += 4 * otstupHoriz - 10;
        for (i = 8; i < 10; i++)
        {
            labels[i].setBounds(nowWidth, otstupVertic + nowHeight, pictures[i].getIconWidth(),
                    pictures[i].getIconHeight());
            nowWidth += pictures[i].getIconWidth();
            textAreas[i].setBounds(nowWidth, nowHeight + (pictures[i].getIconHeight() + 20) / 2, 30, 16);
            nowWidth += 30 + otstupHoriz;
        }
        real.setBounds(130, 300, 80, 16);
        realFunc.setBounds(50, 300, 60, 16);
        add(real);
        add(realFunc);
        for (i = 0; i < labelBeta.length; i++)
        {
            labelBeta[i].setBounds(500, 300 + i*20, 250, 16);
            add(labelBeta[i]);
        }
        textBeta.setBounds(500, 360, 20, 16);
        add(textBeta);
        for (i = 0; i < countLabels; i++)
        {
            add(labels[i]);
            add(textAreas[i]);
        }
    }

    void createBackground()
    {
        JPanel panel = new JPanel();
        panel.setBounds(0, 0, width, height);
        panel.setBackground(Color.WHITE);
        add(panel);
    }

    void createPict()
    {
        for (int i = 0; i < countLabels; i++)
        {
            pictures[i] = new ImageIcon("pict/" + (i + 1) + ".png");
        }
    }

    void createButton()
    {
        int con = 128;
        start.setBounds((width - con) / 2, height - con * 9 / 16 - 40, con, con * 9 / 32);
        start.addActionListener(new ActionListener()
        {
            @Override
            public void actionPerformed(ActionEvent e)
            {
                readAllData();
            }
        });
        add(start);
    }

    void createNev()
    {
        nevTitle.setBounds(20, height - 80, 100, 20);
        nevDate.setBounds(120, height - 80, 200, 20);
        nev2Title.setBounds(20, height - 100, 100, 20);
        nev2Date.setBounds(120, height - 100, 200, 20);
        add(nevTitle);
        add(nevDate);
        add(nev2Title);
        add(nev2Date);
    }

    void readAllData()
    {
        eps1 = Double.parseDouble(textAreas[8].getText());
        eps2 = Double.parseDouble(textAreas[9].getText());
        lengthT = Double.parseDouble(textAreas[6].getText());
        lengthX = Double.parseDouble(textAreas[7].getText());
        N = Integer.parseInt(textAreas[10].getText());
        M = Integer.parseInt(textAreas[11].getText());
        start();
    }

    void start()
    {

        long time = System.currentTimeMillis();
        QuasilinearParabolicProblem obj3 = new QuasilinearParabolicProblem(this);
        obj3.initialization();
        obj3.conditions();
        double tao = obj3.ruleRunge();
        System.out.println("Time with tao: " + (System.currentTimeMillis() - time));
        System.out.println();
        time = System.currentTimeMillis();
        obj3.conditions();
        obj3.findAnsverFromTao(tao / 2.0);
        System.out.println("Time with tao/2: " + (System.currentTimeMillis() - time));
        obj3.findResidual();

        XYSeries series = new XYSeries("Real function");
        XYSeries series1 = new XYSeries("Approximate function");

        for (int i = 0; i <= N; i++)
        {
            series.add((i + 0.0) / N, lengthT * lengthT + (i + 0.0) / N * (i + 0.0) / N
                    + lengthT * lengthT * (i + 0.0) / N * (i + 0.0) / N);
            series1.add((i + 0.0) / N, obj3.vector[i]);
        }

        XYSeriesCollection dataset = new XYSeriesCollection();
        dataset.addSeries(series);
        dataset.addSeries(series1);

        NumberAxis domain = new NumberAxis("x");
        NumberAxis range = new NumberAxis("f(x)");
        XYSplineRenderer r = new XYSplineRenderer(3);
        XYPlot xyplot = new XYPlot(dataset, domain, range, r);
        JFreeChart chart = new JFreeChart(xyplot);
        ChartPanel chartPanel = new ChartPanel(chart){

            @Override
            public Dimension getPreferredSize() {
                return new Dimension(640, 480);
            }
        };
        JFrame frame = new JFrame("MinimalStaticChart");
        // Помещаем график на фрейм
        frame.getContentPane().add(chartPanel);
        frame.setSize(400, 300);
        frame.show();
        frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
    }
}
