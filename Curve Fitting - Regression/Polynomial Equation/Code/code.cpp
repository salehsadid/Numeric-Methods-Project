#include <bits/stdc++.h>
using namespace std;

// File Operation Class
class FileIO{
private:
    ifstream in;
    ofstream out;

public:
    FileIO(){}

    bool openInput(const string &filename = "../Input/input.txt"){
        in.open(filename);
        if(!in.is_open()){
            cerr << "Cannot open input file: " << filename << endl;
            return false;
        }
        return true;
    }

    bool openOutput(const string &filename = "../Output/output.txt"){
        out.open(filename);
        if(!out.is_open()){
            cerr << "Cannot open output file: " << filename << endl;
            return false;
        }
        return true;
    }

    template<typename T>
    void print(const T &val){
        cout << val;
        out  << val;
        out.flush();
    }

    template<typename T>
    void printToFile(const T &val){
        out << val;
        out.flush();
    }

    ifstream& input(){
        return in;
    }

    void closeFiles(){
        if(in.is_open())  in.close();
        if(out.is_open()) out.close();
    }
};

int main()
{
    FileIO io;
    if(!io.openInput()) return 1;
    if(!io.openOutput()) return 1;

    int n;
    io.input() >> n;

    vector<double> x(n), y(n);
    for (int i = 0; i < n; i++)
    {
        io.input() >> x[i] >> y[i];
    }

    double Sx = 0, Sy = 0, Sxx = 0, Sxy = 0, Sxxx = 0, Sxxxx = 0, Sxxy = 0;

    for (int i = 0; i < n; i++)
    {
        double X = x[i], Y = y[i];
        double X2 = X * X;

        Sx += X;
        Sy += Y;
        Sxx += X2;
        Sxy += X * Y;
        Sxxx += X2 * X;
        Sxxxx += X2 * X2;
        Sxxy += X2 * Y;
    }

    // 3x3 system
    // A * coeff = B
    vector<vector<double>> A =
    {
        { double(n), Sx, Sxx },
        { Sx, Sxx, Sxxx },
        { Sxx, Sxxx, Sxxxx }
    };

    vector<double> B = { Sy, Sxy, Sxxy }; // 1D vector for column

    int m = 3; // size of matrix

    // Apply Gauss elimination on augmented matrix [A | B]
    for (int i = 0; i < m; i++)
    {
        double pivot = A[i][i];
        if (fabs(pivot) < 1e-9)
        {
            io.print("Singular matrix\n");
            io.closeFiles();
            return 0;
        }

        // Normalize pivot row
        for (int j = 0; j < m; j++)
        {
            A[i][j] /= pivot;
        }
        B[i] /= pivot;

        // Eliminate other rows
        for (int r = 0; r < m; r++)
        {
            if (r == i)
            {
                continue;
            }

            double factor = A[r][i];

            for (int c = 0; c < m; c++)
            {
                A[r][c] -= factor * A[i][c];
            }
            B[r] -= factor * B[i];
        }
    }

    double a = B[0];
    double b = B[1];
    double c = B[2];

    stringstream ss;
    ss << fixed << setprecision(6);
    ss << "y = " << a << " + " << b << "x + " << c << "x^2\n";
    io.print(ss.str());

    io.closeFiles();
    return 0;
}
