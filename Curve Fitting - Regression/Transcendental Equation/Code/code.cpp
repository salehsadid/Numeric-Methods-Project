#include <bits/stdc++.h>
using namespace std;

int main()
{
    freopen("../Input/input.txt", "r", stdin);
    freopen("../Output/output.txt", "w", stdout);

    int n;
    cin >> n;

    vector<double> x(n), y(n);
    for (int i = 0; i < n; i++)
    {
        cin >> x[i] >> y[i];
    }

    double Sx = 0, Sy = 0, Sxx = 0, Sxy = 0;

    for (int i = 0; i < n; i++)
    {
        double X = x[i];
        double Y = log(y[i]); // Linearize: ln(y) = ln(a) + bx

        Sx += X;
        Sy += Y;
        Sxx += X * X;
        Sxy += X * Y;
    }

    // 2x2 system for Y = A + Bx
    // where A = ln(a), B = b
    vector<vector<double>> A_mat =
    {
        { double(n), Sx },
        { Sx, Sxx }
    };

    vector<double> B_vec = { Sy, Sxy };

    int m = 2; // size of matrix

    // Gauss elimination
    for (int i = 0; i < m; i++)
    {
        double pivot = A_mat[i][i];
        if (fabs(pivot) < 1e-9)
        {
            cout << "Singular matrix\n";
            return 0;
        }

        // Normalize pivot row
        for (int j = 0; j < m; j++)
        {
            A_mat[i][j] /= pivot;
        }
        B_vec[i] /= pivot;

        // Eliminate other rows
        for (int r = 0; r < m; r++)
        {
            if (r == i) continue;

            double factor = A_mat[r][i];
            for (int c = 0; c < m; c++)
            {
                A_mat[r][c] -= factor * A_mat[i][c];
            }
            B_vec[r] -= factor * B_vec[i];
        }
    }

    double A_val = B_vec[0];
    double B_val = B_vec[1];

    double a = exp(A_val);
    double b = B_val;

    cout << fixed << setprecision(6);
    cout << "y = " << a << "e^" << b << "x\n";

    return 0;
}
