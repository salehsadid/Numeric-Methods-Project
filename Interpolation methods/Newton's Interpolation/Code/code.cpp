//Author :: tahmids55
#include<bits/stdc++.h>
using namespace std;
typedef long long ll;

class NewtonInterpolation{
private:
    int n;
    double y;
    vector<double> x , f;
    vector<vector<double>> dif;  // forward difference table
    vector<vector<double>> bif;  // backward difference table
    double h , u ;
    double fores , backres;
    bool forward_calculated;
    bool backward_calculated;

    // Recursive function for forward interpolation
    double recur(double i, double u){
        if(i == 0) return 1;
        if(i == 1) return u;
        return (u - i + 1) / i * recur(i - 1, u);
    }

    // Recursive function for backward interpolation
    double becur(double i, double u){
        if(i == 0) return 1;
        if(i == 1) return u;
        return (u + i - 1) / i * becur(i - 1, u);
    }

    // Calculate forward difference table
    void calculateForwardDifference(){
        dif.assign(n, vector<double>(n, 0));
        for(int i = 0; i < n; i++){
            dif[i][0] = f[i];
        }
        
        int cur = n;
        for(int i = 1; i < n; i++){
            for(int j = 0; j < cur - 1; j++){
                dif[j][i] = dif[j + 1][i - 1] - dif[j][i - 1];
            }
            cur--;
        }
    }

    // Calculate backward difference table
    void calculateBackwardDifference(){
        bif.assign(n, vector<double>(n, 0));
        for(int i = 0; i < n; i++){
            bif[i][0] = f[i];
        }
        
        int cur = 0;
        for(int i = 1; i < n; i++){
            for(int j = n - 1; j > cur; j--){
                bif[j][i] = bif[j][i - 1] - bif[j - 1][i - 1];
            }
            cur++;
        }
    }

    // Calculate forward interpolation result
    void calculateForwardResult() {
        h = (x[1] - x[0]);
        u = (y - x[0]) / h;
        
        fores = 0.0;
        for(int i = 0; i < n; i++){
            fores += dif[0][i] * recur(i, u);
        }
        forward_calculated = true;
    }

    // Calculate backward interpolation result
    void calculateBackwardResult() {
        h = (x[1] - x[0]);
        u = (y - x[n - 1]) / h;
    
        backres = 0.0;
        for(int i = 0; i < n; i++){
            backres += bif[n - 1][i] * becur(i, u);
        }
        backward_calculated = true;
    }

public:
    // Constructor
    NewtonInterpolation() : forward_calculated(false), backward_calculated(false) {}

    void setData(int n_val, double y_val, vector<double>& x_val, vector<double>& f_val){
        n = n_val;
        y = y_val;
        x = x_val;
        f = f_val;
    }

    //Perform Newton Forward Interpolation
    void newtonForward(){
        calculateForwardDifference();
        calculateForwardResult();
    }

    //Perform Newton Backward Interpolation
    void newtonBackward(){
        calculateBackwardDifference();
        calculateBackwardResult();
    }

    //Write forward difference table
    void writeForwardDifferenceTable(){
        cout << "\n=== Forward Difference Table ===\n";
        for(int i = 0; i < n; i++){
            for(int j = 0; j < n; j++){
                if(dif[i][j] != 0 || j == 0){
                    cout << fixed << setprecision(2) << dif[i][j] << " ";
                }
            }
            cout << "\n";
        }
    }

    //Write backward difference table
    void writeBackwardDifferenceTable(){
        cout << "\n=== Backward Difference Table ===\n";
        for(int i = 0; i < n; i++){
            for(int j = 0; j < n; j++){
                if(bif[i][j] != 0 || j == 0){
                    cout << fixed << setprecision(2) << bif[i][j] << " ";
                }
            }
            cout << "\n";
        }
    }

    // Write forward result
    void writeForwardResult(){
        if(!forward_calculated){
            return;
        }
        cout << "\n=== Newton Forward Result ===\n";
        cout << "f(" << y << ") = " << fixed << setprecision(2) << fores << "\n";
    }

    // Write backward result
    void writeBackwardResult(){
        if(!backward_calculated){
            return;
        }
        cout << "\n=== Newton Backward Result ===\n";
        cout << "f(" << y << ") = " << fixed << setprecision(6) << backres << "\n";
    }
};

// Controling function 
void Netwon_Interpolations_method(){
    NewtonInterpolation newton;

    // Read data from file
    int n;
    double y;
    vector<double> x, f;
    
    if(!(cin >> n >> y)) return;
    
    x.resize(n);
    f.resize(n);
    for(int i = 0; i < n; i++){
        cin >> x[i] >> f[i];
    }
    newton.setData(n, y, x, f);
    
    int mtd_choice;
    // cout << "\n=== Select Interpolation Method ===\n";
    // cout << "1. Newton Forward Interpolation\n";
    // cout << "2. Newton Backward Interpolation\n";
    // cout << "Enter your choice (1 or 2): ";
    
    // Read choice from input file
    if(cin >> mtd_choice){
        if(mtd_choice == 1){
            newton.newtonForward();
            newton.writeForwardDifferenceTable();
            newton.writeForwardResult();
        } else if(mtd_choice == 2){
            newton.newtonBackward();
            newton.writeBackwardDifferenceTable();
            newton.writeBackwardResult();
        } else {
            cout << "Invalid choice in input file!\n";
        }
    } else {
        // Default to Forward if no choice provided? Or error?
        // Let's assume Forward if not specified, or just exit.
        cout << "No method choice found in input file.\n";
    }

    cout << "\nThank you for using Newton Interpolation!\n";
}

int main(){
    freopen("../Input/input.txt", "r", stdin);
    freopen("../Output/output.txt", "w", stdout);

    cout << "\n === Methods === \n";
    Netwon_Interpolations_method();
    return 0;
}
