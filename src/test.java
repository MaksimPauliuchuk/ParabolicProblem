import java.util.LinkedList;

public class test {
        public static void main(String[] args) {
            System.out.println(solveK2(1,2,5));
            MatchParser matchParser = new MatchParser();
            matchParser.setVariable("x", 2.0);
            matchParser.setVariable("t", 5.1);
            try
            {
                System.out.println(matchParser.Parse("t*x"));
            }
            catch (Exception e)
            {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
            
        }
    static double solveMu1(double t, double x){
        StringBuffer mu1 = new StringBuffer("x*x");
        for(int i = 0; i < mu1.length(); i++){
            if(mu1.charAt(i) == 't') mu1.replace(i, i + 1, t + "");
            if(mu1.charAt(i) == 'x') mu1.replace(i, i + 1, x + "");
        }
        return calculate(mu1.toString());
    }
    static double solveK2(double t, double x, double u){
        StringBuffer g = new StringBuffer("2*t-2*t*u-2*t*t*t*u-2*u*u*x-2*t*t*u*u*x-2*t*x*x-8*t*t*t*x*x-4*t*t*t*t*t*x*x-8*u*x*x*x-16*t*t*u*x*x*x-8*t*t*t*t*u*x*x*x");
        for(int i = 0; i < g.length(); i++){
            if(g.charAt(i) == 't') g.replace(i, i + 1, t + "");
            if(g.charAt(i) == 'x') g.replace(i, i + 1, x + "");
            if(g.charAt(i) == 'u') g.replace(i, i + 1, u + "");
        }
        return calculate(g.toString());
    }
    static double calculate(String pattern){
        StringBuffer s = new StringBuffer();
        for(int i = 0; i < pattern.length(); i++)
            if(pattern.charAt(i) != ' ') s.append(pattern.charAt(i));

        StringBuffer s1 = new StringBuffer(s);

        for (int i = 0; i < s1.length(); i++) {
            char c = s1.charAt(i);
            if(c == '-' || c == '+' || c == '*' || c == '/')
                s1.replace(i,i + 1," ");
        }
        String doubleString[] = s1.toString().split(" ");
        s = new StringBuffer("+" + s);
        int i = 0;
        int k = 0;
        LinkedList<Double> stack = new LinkedList<Double>();
        while (i < s.length()){
            char c = s.charAt(i);
            switch (c){
                case '+':  stack.add(0,Double.parseDouble(doubleString[k]));k++;break;
                case '-':  stack.add(0,-1*Double.parseDouble(doubleString[k]));k++;break;
                case '*':  stack.add(0,Double.parseDouble(doubleString[k]) * stack.remove(0));k++;break;
                case '/':  stack.add(0,stack.remove(0) / Double.parseDouble(doubleString[k]));k++;break;
            }
            i++;
        }
        double res = 0;
        for(Double d : stack){
            res += d;
        }
        return res;
    }
    }

