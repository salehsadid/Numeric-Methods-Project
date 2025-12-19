//author: tahmids55
#include <bits/stdc++.h>
using namespace std;

float dydx(float x, float y){
    return (x - y) / 2;
}

float rungeKutta(float x0, float y0, float x, float h){
    int n = (int)((x - x0) / h);
    float y = y0;

    for(int i = 0; i < n; i++){
        float k1 = h * dydx(x0, y);
        float k2 = h * dydx(x0 + 0.5f * h, y + 0.5f * k1);
        float k3 = h * dydx(x0 + 0.5f * h, y + 0.5f * k2);
        float k4 = h * dydx(x0 + h, y + k3);

        y += (k1 + 2*k2 + 2*k3 + k4) / 6.0f;
        x0 += h;
    }

    return y;
}

int main(){
    freopen("../Input/input.txt", "r", stdin);
    freopen("../Output/output.txt", "w", stdout);

    float x0, y0, x, h;
    
    // Reading from file
    if(cin >> x0 >> y0 >> x >> h){
        cout << "=== Runge Kutta Method ===\n";
        cout << "Initial Condition (x0, y0): (" << to_string(x0) << ", " << to_string(y0) << ")\n";
        cout << "Target x: " << to_string(x) << "\n";
        cout << "Step size h: " << to_string(h) << "\n";

        float result = rungeKutta(x0, y0, x, h);
        
        cout << "\nResult y(" << to_string(x) << ") = " << to_string(result) << "\n";
    } else {
        cout << "Error reading input data." << endl;
    }

    return 0;
}
