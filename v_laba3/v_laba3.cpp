#include <vector>
#include <fstream>
#include <iostream>

using namespace std;

struct Data {
    vector <double> X;
    vector <double> Y;
};

struct BigData{
    vector <double> a;
    vector <double> b;
    vector <double> c;
    vector <double> d;
};

void PrintVector(vector <double> tmp) {
    for (int i = 0; i < tmp.size(); i++) {
        cout << " " << tmp[i];
    }
    cout << ";\n";
}

void PrintData(Data data_for_print) {   //function for fast check data
    cout << "X:";
    PrintVector(data_for_print.X);

    cout << "Y:";
    PrintVector(data_for_print.Y);
}

void PrintBigData(BigData data_for_print) {   //function for fast check data
    cout << "a:";
    PrintVector(data_for_print.a);

    cout << "b:";
    PrintVector(data_for_print.b);

    cout << "c:";
    PrintVector(data_for_print.c);

    cout << "d:";
    PrintVector(data_for_print.d);

}

void CountA(int n, vector<double>& A, vector <double> Y) {
    for (int i = 0; i < n+1; i++) {
        A[i] = Y[i];
    }
}

void CountH(int n, vector <double>& H, vector <double> X) {
    for (int i = 0; i < n; i++) {
        H[i] = X[i + 1] - X[i];
    }
}

void CountAlpha(int n,vector <double>& alpha, vector <double> A, vector <double> H) {
    for (int i = 1; i < n; i++) {
        alpha[i] = (3 / H[i]) * (A[i + 1] - A[i]) - (3 / H[i - 1]) * (A[i] - A[i - 1]);
    }
}

void CountLmuZ(int n, vector <double>& l, vector <double> X, vector <double> H, vector <double>& mu, vector <double>& Z, vector <double> alpha) {
    l[0] = 1;
    for (int i = 1; i < n; i++) {
        l[i] = 2 * (X[i + 1] - X[i - 1]) - H[i - 1] * mu[i - 1];
        mu[i] = H[i] / l[i];
        Z[i] = (alpha[i] - H[i - 1] * Z[i - 1]) / l[i];
    }
    l[n] = 1; 
    Z[n] = 0;
}

void CountC(int n, vector<double>& c, vector<double> z, vector<double> mu) {
    for (int i = n - 1; i >= 0; i--) {
        c[i] = z[i] - mu[i] * c[i + 1];
    }
}

void CountBD(int n, vector <double>& b, vector <double>& d, vector <double> a, vector <double> h, vector <double> c) {
    for (int i = n - 1; i >= 0; i--) {
        b[i] = (a[i + 1] - a[i]) / h[i] - h[i] * (c[i + 1] + 2 * c[i]) / 3;
        d[i] = (c[i + 1] - c[i]) / (3 * h[i]);
    }
}

double FindPoint(double a, vector<double> X) {
    int i = 0;
    double num = 0;
    while (X[i] <= a) {
        num = i;
        i++;
    }
    return num;
}

class CubicSpline {
    Data Points;
    BigData Coefficients;
public:

    /*** Конструктор, в котором строится кубический сплайн ***/
    explicit CubicSpline(const Data& info) {
        Points.X = info.X;
        Points.Y = info.Y;

        int n = Points.X.size() - 1;
        Coefficients.a.resize(n + 1, 0.0);
        CountA(n, Coefficients.a, Points.Y);

        Coefficients.b.resize(n, 0.0);
        Coefficients.d.resize(n, 0.0);

        vector <double> h(n, 0.0);
        CountH(n, h, Points.X);

        vector <double> alpha (n,0.0);
        CountAlpha(n, alpha, Coefficients.a, h);

        Coefficients.c.resize(n + 1, 0.0);
        vector<double> l(n + 1, 0.0);
        vector<double> mu(n + 1, 0.0);
        vector<double> z(n + 1, 0.0);

        CountLmuZ(n, l, Points.X, h, mu, z, alpha);

        Coefficients.c[n] = 0;
        CountC(n, Coefficients.c, z, mu);

        CountBD(n, Coefficients.b, Coefficients.d, Coefficients.a, h, Coefficients.c);
    };

    /*** Метод, выполняющий подсчет интерполянта в точке ***/
    [[nodiscard]] double interpolate(double x) const {
        double num = FindPoint(x, Points.X);

        double a_p = Coefficients.a[num];
        double b_p = Coefficients.b[num];
        double c_p = Coefficients.c[num];
        double d_p = Coefficients.d[num];

        double ans = a_p + b_p * (x - Points.X[num]) + c_p * (x - Points.X[num]) * (x - Points.X[num]) + d_p * (x - Points.X[num]) * (x - Points.X[num]) * (x - Points.X[num]);
        return ans;
    };

    BigData GetCoefficient() {
        return Coefficients;
    }
};

vector <double> SplitSegment(int NumPoints, double A, double B) {
    vector <double> resault;
    double step = (B - A) / (NumPoints - 1);

    for (int i = 0; i < NumPoints; i++) {
        resault.push_back(A + i * step);
    }

    return resault;
}

double Mistake(const vector <double> points, const CubicSpline interpolant) {
    double resault = 0;
    double tmp = 0;

    for (int i = 0; i < points.size() - 1; i++) {
        tmp = abs( interpolant.interpolate(points[i]) - cos(points[i]) );
        if (tmp > resault) {
            resault = tmp;
        }
    }
    return resault;
}

Data CosValue(int NumPoints, int BeginSegment, double EndSegment) {
    Data resault;
    resault.X = SplitSegment(NumPoints, BeginSegment, EndSegment);

    for (int i = 0; i < NumPoints; i++) {
        resault.Y.push_back(cos(resault.X[i]));
    }

    return resault;
}

int main() {
    /*Data example;
    example.X = { 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1 };
    example.Y = { 0.000, 0.033, 0.067, 0.100, 0.133, 0.166, 0.199, 0.231, 0.264, 0.296, 0.327 };
    double dot = 0.95;

    CubicSpline Spline = CubicSpline(example);
    PrintBigData(Spline.GetCoefficient());

    cout << Spline.interpolate(dot);*/

    vector <double> Points = { 2,4,8,16,32,64,128 };

    ofstream fileX, fileY;
    fileX.open("X.txt"); 
    fileY.open("Y.txt");
    for (int i = 0; i < Points.size(); i++) {
        Data CosPoints = CosValue(Points[i], 0, 3);  //points for interpolating 
        CubicSpline CosInterpolator = CubicSpline(CosPoints);

        fileX << Points[i] << endl;
        fileY << Mistake(SplitSegment(1000, 0, 3), CosInterpolator) << endl;
    }
    fileX.close();
    fileY.close();
}