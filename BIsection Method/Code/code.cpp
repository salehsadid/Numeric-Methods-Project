#include <bits/stdc++.h>
using namespace std;

long double an1, an2, an3, an4;

long double xmax()
{
    long double val = ((an2 * an2) / (an1 * an1)) - 2 * (an3 / an1);
    return sqrt(val);
}

long double func(long double x)
{
    return an1*x*x*x*x + an2*x*x*x + an3*x*x + an4*x;
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
    
    // Read coefficients and error tolerance from input file
    input >> an1 >> an2 >> an3 >> an4;
    long double e;
    input >> e;
    
    output << "Equation: f(x) = (" << an1 << ")x^4 + (" << an2 << ")x^3 + (" << an3 << ")x^2 + (" << an4 << ")x" << endl;
    output << endl;
    output << "Coefficients are: " << an1 << " " << an2 << " " << an3 << " " << an4 << endl;
    output << "Error tolerance: " << e << endl;
    output << endl;
    
    // Find ranges where roots exist
    vector<pair<long double, long double>> range;
    long double range1 = -abs(xmax());
    long double range2 = abs(xmax());
    
    while(range1 <= range2)
    {
        if((func(range1) * func(range1 + 0.5)) < 0)
        {
            range.push_back({range1, range1 + 0.5});
            output << "Range found: " << range1 << " to " << range1 + 0.5 << endl;
        }
        range1 += 0.5;
    }
    output << endl;
    
    // Apply Bisection Method for each range
    for(int i = 0; i < (int)range.size(); i++)
    {
        int it = 0;
        long double a = range[i].first;
        long double b = range[i].second;
        long double c = (a + b) / 2;
        
        output << "Finding Root " << i + 1 << ":" << endl;
        output << fixed << setprecision(6);
        
        while(abs((b - a) / 2) > e)
        {
            it++;
            c = (a + b) / 2;
            
            if(func(c) == 0.0)
            {
                break;
            }
            else if((func(a) * func(c)) < 0)
            {
                b = c;
            }
            else
            {
                a = c;
            }
        }
        
        output << "Root " << i + 1 << ": " << c << endl;
        output << "Iterations: " << it << endl;
        output << endl;
    }
    
    input.close();
    output.close();
    
    return 0;
}
