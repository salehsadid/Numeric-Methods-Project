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
        Sx += x[i];
        Sy += y[i];
        Sxx += x[i] * x[i];
        Sxy += x[i] * y[i];
    }

    // 2x2 system
    vector<vector<double>> A =
    {
        { double(n), Sx },
        { Sx, Sxx }
    };

    vector<double> B = { Sy, Sxy };

    int m = 2; // size of matrix

    // Gauss elimination
    for (int i = 0; i < m; i++)
    {
        double pivot = A[i][i];
        if (fabs(pivot) < 1e-9)
        {
            cout << "Singular matrix\n";
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
            if (r == i) continue;

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

    cout << fixed << setprecision(6);
    cout << "y = " << a << " + " << b << "x\n";

    return 0;
}
