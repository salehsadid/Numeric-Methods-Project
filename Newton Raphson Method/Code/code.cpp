#include <bits/stdc++.h>
using namespace std;

long double f(long double x)
{
    return exp(-x)*cos(x);
}

long double df(long double x)
{
    return -exp(-x)*(sin(x)+cos(x));
}

int main()
{
    // Open input and output files
    ifstream input("../Input/input.txt");
    ofstream output("../Output/output.txt");
    
    if(!input.is_open())
    {
        return 1;
    }
    
    // Read initial guess and tolerances from input file
    long double x0, step, tol;
    int maxIter;
    input >> x0 >> step >> tol >> maxIter;
    
    output << "The equation is: f(x) = e^(-x) * cos(x)" << endl;
    output << endl;
    output << "Initial guess (x0): " << x0 << endl;
    output << "Step tolerance: " << step << endl;
    output << "Function value tolerance: " << tol << endl;
    output << "Maximum iterations: " << maxIter << endl;
    output << endl;
    output << fixed << setprecision(6);
    
    long double xn = x0;
    
    for(int i = 1; i <= maxIter; i++)
    {
        long double fx = f(xn);
        long double dfx = df(xn);
        
        if(dfx == 0)
        {
            output << "Error: Derivative is zero at iteration " << i << endl;
            output << "No solution found" << endl;
            input.close();
            output.close();
            return 0;
        }
        
        long double xn1 = xn - (fx / dfx);
        long double sc1 = fabs(xn1 - xn);
        long double sc2 = fabs(f(xn1));
        
       
        if(sc1 < step || sc2 < tol)
        {
            output << endl;
            output << "Root found: " << xn1 << endl;
            output << "Number of iterations: " << i << endl;
            break;
        }
        
        xn = xn1;
        
        if(i == maxIter)
        {
            output << endl;
            output << "Maximum iterations reached" << endl;
            output << "Approximate root: " << xn1 << endl;
        }
    }
    
    input.close();
    output.close();
    
    return 0;
}
