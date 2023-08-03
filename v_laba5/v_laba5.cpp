#include <iostream>
#include <functional>
#include <vector>
#include <optional>
#include <fstream>


using namespace std;

//пробная функция
double f(double x) {
    return x * x - 4;
}

//уравнение Кеплера
double Kepler(double x) {
    return x - 0.1 * sin(x) - 3.14 / 4;
}

//уравнение 2
double eq2(double x) {
    return tan(x) - 4 * x / 3.14;
}

//уравнение 3
double eq3(double x) {
    return log(cosh(x));
}




vector<double> gauss(double** a, double* y, int n)
{
    double max;
    int k, index;
    const double eps = std::numeric_limits<double>::min();  // точность
    vector<double> x;
    k = 0;
    while (k < n)
    {
        // Поиск строки с максимальным a[i][k]
        max = abs(a[k][k]);
        index = k;
        for (int i = k + 1; i < n; i++)
        {
            if (abs(a[i][k]) > max)
            {
                max = abs(a[i][k]);
                index = i;
            }
        }

        // переносим строку с максимальным a[i][k] наверх
        for (int j = 0; j < n; j++)
        {
            double temp = a[k][j];
            a[k][j] = a[index][j];
            a[index][j] = temp;
        }

        // также переносим соответствующий этой строке элемент в столбце решений наверх
        double temp = y[k];
        y[k] = y[index];
        y[index] = temp;

        // Нормализация уравнений
        for (int i = k; i < n; i++)
        {
            double temp = a[i][k];
            if (abs(temp) < eps) continue; // для нулевого коэффициента пропустить
            for (int j = 0; j < n; j++)
                a[i][j] = a[i][j] / temp; // делим всю строку на максимальный элемент
            y[i] = y[i] / temp;
            if (i == k)  continue; // уравнение не вычитать само из себя
            for (int j = 0; j < n; j++)
                a[i][j] = a[i][j] - a[k][j]; // вычитаем верхнюю строку из всех нижних
            y[i] = y[i] - y[k];
        }
        k++;
    }
    // обратная подстановка (имеем верхнетреугольную матрицу, начинаем искать x(i) с нижней строки, двигаясь наверх
    for (k = n - 1; k >= 0; k--)
    {
        x.insert(x.begin(), y[k]);
        for (int i = 0; i < k; i++)
            y[i] = y[i] - a[i][k] * y[k];
    }
    return x;
}

[[nodiscard]] vector<double> firstDer(const vector<double>& points) noexcept { //nodiscard указывает, что возвращаемое функцией значение нельзя игнорировать и нужно сохранить в какую-либо переменную

    int n = points.size();
    // выделяем место под матрицу СЛАУ
    // заполняем матрицу
    double** A = new double* [n];
    for (int i = 0; i < n; i++) {
        A[i] = new double[n];
        for (int j = 0; j < n; j++) {
            A[i][j] = pow(points[j], i);
        }
    }

    // столбец решений
    double* y = new double[n];
    for (int i = 0; i < n; i++) {
        y[i] = 0;
    };
    y[1] = 1;

    //решаем СЛАУ методом Гаусса
    vector<double> x = gauss(A, y, n);

    delete[] y;

    for (int i = 0; i < n; i++) {
        delete[] A[i];
    };
    delete[] A;

    return x;
}

double dfdx(double h, vector<double> C, vector<double>& points, double x, const function<double(double)>& F) {
    double y = 0;
    for (int i = 0; i < points.size(); i++) {
        y += F(x + points[i] * h) * C[i] / h;
    }
    return y;
}
/**
* Функция для решения уравнения func(x) = 0 при помощи метода половинного деления
* a - левая граница отрезка локализации корня
* b - правая гранича отрезка локализации корня
* func - функция, корень которой нужно найти
* numberOfIterations - количество итераций, которое должен освершить метод
* ДОСТАТОЧНЫЕ УСЛОВИЯ СХОДИМОСТИ ПРОВЕРЯТЬ НЕ НУЖНО!
**/
[[nodiscard]] double bisectionMethod(double a, double b, const function<double(double)>& f, unsigned numberOfIterations) noexcept {
    double c; // корень
    for (int i = 0; i < numberOfIterations; i++) {
        c = (a + b) / 2;
        if (f(a) * f(c) > 0) a = c;
        else b = c;
    }
    return c;
}

/**
* Функция для решения уравнения func(x) = 0 при помощи метода простой итерации с релаксацией
* inital - начальное приближение для решения
* func - функция, корень которой нужно найти
* tau - параметр в методе простой итерации
* numberOfIterations - количество итераций, которое должен освершить метод
**/
[[nodiscard]] double simpleIterationMethod(double inital, const function<double(double)>& func, double tau, unsigned numberOfIterations) noexcept {
    double  x = inital;
    for (int i = 0; i < numberOfIterations; i++) {
        x = x + tau * func(x);
    }
    return x;
}

/**
* Функция для решения уравнения func(x) = 0 при помощи метода Ньютона
* inital - начальное приближение для решения
* func - функция, корень которой нужно найти
* numberOfIterations - количество итераций, которое должен освершить метод
* Производную стоит оценить численно (предложите схему не ниже 2 порядка самостоятельно)
**/
[[nodiscard]] double newtonMethod(double inital, const function<double(double)>& func, unsigned numberOfIterations) noexcept {
    vector<double> points = { -2, -1, 0, 1, 2 };
    vector<double> c = firstDer(points);
    double x = inital;
    for (int i = 0; i < numberOfIterations; i++) {
        x = x - (func(x) / dfdx(0.001, c, points, x, func));
    }
    return x;
}

int main()
{
    //Уравнение Кеплера
   // (метод Ньютона)
    ofstream file1_x, file1_y;
    file1_x.open("K_Nx"); // окрываем файл для записи
    file1_y.open("K_Ny");
    for (int i = 0; i < 10; i++) {
        double Kepler_Newton = newtonMethod(5, Kepler, i);
        file1_x << i << endl;
        file1_y << abs(0.861265 - Kepler_Newton) << endl;
    }
    file1_x.close();
    file1_y.close();

    //Уравнение Кеплера
    // (метод половинного деления)
    ofstream file2_x, file2_y;
    file2_x.open("K_mpix"); // окрываем файл для записи
    file2_y.open("K_mpiy");
    for (int i = 0; i < 10; i++) {
        double Kepler_mpi = simpleIterationMethod(1, Kepler, -1.8 + i * 0.05, 10);
        file2_x << -1.8 + i * 0.05 << endl;
        file2_y << abs(0.861265 - Kepler_mpi) << endl;
    }
    file1_x.close();
    file1_y.close();





    //Уравнение tg
    //с начальным приближением в окрестости 1 (все 3 метода)
      // (метод Ньютона)
    ofstream file3_x, file3_y;
    file3_x.open("tg_Nx"); // окрываем файл для записи
    file3_y.open("tg_Ny");
    for (int i = 0; i < 10; i++) {
        double eq2_Newton = newtonMethod(1, eq2, i);
        file3_x << i << endl;
        file3_y << abs(0.78539816 - eq2_Newton) << endl;
    }
    file3_x.close();
    file3_y.close();

    // (метод половинного деления)
    ofstream file4_x, file4_y;
    file4_x.open("tg_mpix"); // окрываем файл для записи
    file4_y.open("tg_mpiy");
    for (int i = 0; i < 10; i++) {
        double eq2_mpi = simpleIterationMethod(1, eq2, -1.2 + i * 0.05, 10);
        file4_x << -1.2 + i * 0.05 << endl;
        file4_y << abs(0.78539816 - eq2_mpi) << endl;
    }
    file4_x.close();
    file4_y.close();

    // (дихотомия)
    ofstream file5_x, file5_y;
    file5_x.open("tg_Dx"); // окрываем файл для записи
    file5_y.open("tg_Dy");
    for (int i = 0; i < 10; i++) {
        double eq2_D = bisectionMethod(0.5, 1, eq2, i + 1);
        file5_x << i + 1 << endl;
        file5_y << abs(0.78539816 - eq2_D) << endl;
    }
    file5_x.close();
    file5_y.close();






    // ln(ch) Уравнение(все 3 метода)
          // (метод Ньютона)
    ofstream file6_x, file6_y;
    file6_x.open("ln_Nx"); // окрываем файл для записи
    file6_y.open("ln_Ny");
    for (int i = 0; i < 10; i++) {
        double eq3_Newton = newtonMethod(0.5, eq3, i);
        file6_x << i << endl;
        file6_y << abs(0 - eq3_Newton) << endl;
    }
    file6_x.close();
    file6_y.close();

    // (метод половинного деления)
    ofstream file7_x, file7_y;
    file7_x.open("ln_mpix"); // окрываем файл для записи
    file7_y.open("ln_mpiy");
    for (int i = 0; i < 10; i++) {
        double eq3_mpi = simpleIterationMethod(0.5, eq3, -4.5 + i * 0.1, 10);
        file7_x << -4.5 + i * 0.1 << endl;
        file7_y << abs(0 - eq3_mpi) << endl;
    }
    file7_x.close();
    file7_y.close();

    // (дихотомия)
    ofstream file8_x, file8_y;
    file8_x.open("ln_Dx"); // окрываем файл для записи
    file8_y.open("ln_Dy");
    for (int i = 0; i < 10; i++) {
        double eq3_D = bisectionMethod(-0.5, 0.5, eq3, i + 1);
        file8_x << i + 1 << endl;
        file8_y << abs(0 - eq3_D) << endl;
    }
    file8_x.close();
    file8_y.close();

    double res2 = simpleIterationMethod(1, eq2, -0.001, 10000);
    cout << res2 << endl;

    double res3 = newtonMethod(0.1, eq3, 100);
    cout << res3 << endl;
}


