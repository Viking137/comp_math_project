#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iterator>
#include <string>

using namespace std;

struct Data {                         //create struct for more comfortable use points
    vector <double> X;
    vector <double> Y;
};

void PrintVector(vector <double> tmp) {
    for (int i = 0; i < tmp.size(); i++) {
        cout << " " << tmp[i];
    }
    cout << ";\n";
}

void PrintData (Data data_for_print){   //function for fast check data
    cout << "X:";
    PrintVector(data_for_print.X);

    cout << "Y:";
    PrintVector(data_for_print.Y);
}

void DividedDifference (Data& points) {   //function for count devided difference :D
    int n = points.X.size();
    int a = 1;

    for (int i = n - 1; i >= 0; i--) {
        for (int j = 1; j <= i; j++) {
            points.Y[n - j] = (points.Y[n - j] - points.Y[n - 1 - j]) / (points.X[n - j] - points.X[n - j - a]);
        }
        a++;
    }
}

vector <double> SplitSegment(int NumPoints, double A, double B) {
    vector <double> resault;
    double step = (B - A) / (NumPoints - 1);

    for (int i = 0; i < NumPoints; i++) {
        resault.push_back(A + i * step);
    }

    return resault;
}

class NewtonianInterpolator {
    Data nodes;
public:

    explicit NewtonianInterpolator(const Data info) {
        nodes.X = info.X;
        nodes.Y = info.Y;
        
        DividedDifference(nodes);
    };

    [[nodiscard]] double Interpolate(double a) const {
        int n = nodes.Y.size();
        double resault = 0;

        for (int i = 2; i <= n; i++) {
            if (i == 2) {
                resault = nodes.Y[n - 2] + nodes.Y[n - 1] * (a - nodes.X[n - 2]);
            }
            else {
                resault = resault* (a - nodes.X[n - i]) + nodes.Y[n - i];
            }
        }
        return resault;
    };
};

double Mistake(const vector <double> points, const NewtonianInterpolator interpolant) {
    double resault = 0;
    double tmp = 0;

    for (int i = 0; i < points.size(); i++) {
        tmp = abs(interpolant.Interpolate(points[i]) - cos(points[i]));
        if (tmp > resault) {
            resault = tmp;
        }
    }
    return resault;
}

vector<double> NodesCheb(int NumPoints, double Begin, double End) {
    vector<double> resault;

    for (int k = 1; k <= NumPoints; k++) {
        resault.push_back((Begin + End) / 2 + (End - Begin) / 2 * cos(((2 * k - 1) * M_PI) / (2 * NumPoints)));
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

Data CosValueCheb(int BeginSegment, double EndSegment) {
    Data resault;
    resault.X = NodesCheb(5, BeginSegment, EndSegment);

    for (int i = 0; i < 5; i++) {
        resault.Y.push_back(cos(resault.X[i]));
    }

    return resault;
}

int main()
{
    vector <double> EndSegments = { 8,4,2,1,0.5,0.25,0.125 };

    //3 points
    ofstream fileX3, fileY3;
    fileX3.open("X3.txt"); // окрываем файл для записи
    fileY3.open("Y3.txt");
    for (int i = 0; i < EndSegments.size(); i++) {
        Data CosPoints = CosValue(3, 0, EndSegments[i]);  //points for interpolating 
        NewtonianInterpolator CosInterpolator = NewtonianInterpolator(CosPoints);

        fileX3 << log(EndSegments[i]) << endl;
        fileY3 << log(Mistake(SplitSegment(1000, 0, EndSegments[i]), CosInterpolator)) << endl;
    }
    fileX3.close();
    fileY3.close();

    //5 points
    ofstream fileX5, fileY5;
    fileX5.open("X5.txt"); // окрываем файл для записи
    fileY5.open("Y5.txt");
    for (int i = 0; i < EndSegments.size(); i++) {
        Data CosPoints = CosValue(5, 0, EndSegments[i]);  //points for interpolating 
        NewtonianInterpolator CosInterpolator = NewtonianInterpolator(CosPoints);

        fileX5 << log(EndSegments[i]) << endl;
        fileY5 << log(Mistake(SplitSegment(1000, 0, EndSegments[i]), CosInterpolator)) << endl;
    }
    fileX5.close();
    fileY5.close();

    //10 points
    ofstream fileX10, fileY10;
    fileX10.open("X10.txt"); // окрываем файл для записи
    fileY10.open("Y10.txt");
    for (int i = 0; i < EndSegments.size(); i++) {
        Data CosPoints = CosValue(10, 0, EndSegments[i]);  //points for interpolating 
        NewtonianInterpolator CosInterpolator = NewtonianInterpolator(CosPoints);

        fileX10 << log(EndSegments[i]) << endl;
        fileY10 << log(Mistake(SplitSegment(1000, 0, EndSegments[i]), CosInterpolator)) << endl;
    }
    fileX10.close();
    fileY10.close();
    
    //5 points Cheb
    ofstream fileXCheb, fileYCheb;
    fileXCheb.open("XCheb.txt"); // окрываем файл для записи
    fileYCheb.open("YCheb.txt");
    for (int i = 0; i < EndSegments.size(); i++) {
        Data CosPoints = CosValueCheb(0, EndSegments[i]);  //points for interpolating 
        NewtonianInterpolator CosInterpolator = NewtonianInterpolator(CosPoints);

        fileXCheb << log(EndSegments[i]) << endl;
        fileYCheb << log(Mistake(SplitSegment(1000, 0, EndSegments[i]), CosInterpolator)) << endl;
    }
    fileXCheb.close();
    fileYCheb.close();
}
