#include <iostream>  //realization
#include "linalg.h"
#include <cstring>
#include <math.h>
#include <vector>
#define eps 1.1102230246251565e-16
using namespace std;

namespace MyLinearAlgebra {
// Векторы

// Конструктор по умолчанию
TVector::TVector() : n(0), data(NULL) {}

// Конструктор по количеству элементов
TVector::TVector(int n) : n(0), data(NULL) {
    resize(n);
}

// Конструктор копирования
TVector::TVector(const TVector& rvalue) : n(0), data(NULL) {
    (*this) = rvalue;
}

// Оператор присваивания
TVector& TVector::operator = (const TVector& rvalue) {
    // Если левый операнд не совпадает с правым
    if (this != &rvalue) {
        // Если размер левого операнда не совпадает с правым, то память выделяется заново
        if (n != rvalue.n) {
            // Если память уже была выделена, удалить её
            if (data) { delete[] data;  }
            // Выделение новой памяти
            data = new double[rvalue.n];
            // Сохранение нового размера
            n = rvalue.n;
        }
        // Перенос данных из правого операнда в левый
        memcpy(data, rvalue.data, sizeof(double) * n);
    }
    // Возврат ссылки на левый операнд для возможности цепочки присваиваний
    return *this;
}

// Деструктор
TVector::~TVector() {
    // Если блок данных ранее был инициализирован, удаляем его
    if (data) {
        delete[] data;
        n = 0;
        data = NULL;
    }
}

// Функция задания кол-ва элементов вектора
void TVector::resize(int n) {
#ifdef _DEBUG
    if (n < 0)
        throw 1;
#endif
    // Если новый размер совпадает со старым - выходим
    if (n == this->n) return;
    // Новый блок памяти
    double *newData = new double[n];
    // Если блок данных ранее был инициализирован...
    if (data) {
        // Минимальный из старого и нового размера блока
        int min_n = (this->n < n) ? this->n : n;
        // Перенос данных из старого блока в новый
        memcpy(newData, data, sizeof(double)*min_n);
        // Удаление старого блока
        delete[] data;
    }
    // Прикрепление нового блока к объекту вектора
    data = newData;
    // Сохранение нового размера
    this->n = n;
}

// Оператор сложения векторов
TVector TVector::operator + (const TVector& arg) const {
#ifdef _DEBUG
    if (n != arg.n)
        throw 1;
#endif
    TVector V( n );
    for (int i = 0; i < n; i++)
        V[i] = data[i] + arg[i];
    return V;
}

//Операция унарный минус
TVector TVector::operator - () const {
    TVector V(n);
    for (int i=0; i<n; i++)
        V[i] = -data[i];
    return V;
}

//Операция вычитания векторов
TVector TVector::operator - (const TVector& arg) const {
    /* #ifdef _DEBUG
if (n != arg.n)
throw 1;
#endif */
    TVector V( n );
    for (int i = 0; i < n; i++)
        V[i] = data[i] - arg[i];
    return V;
}

//Операция умножения вектора на число
TVector TVector::operator * (double arg) const {
    TVector V( n );
    for (int i = 0; i < n; i++)
        V[i] = data[i]*arg;
    return V;
}

//Операция скалярного умножения векторов
double TVector::operator * (const TVector& arg) const {
    //TVector V(n);
    double R = 0;
    for (int i = 0; i < n; i++)
        R = R + data[i]*arg[i];
    return R;
}

// Оператор умножения вектора на матрицу
TVector TVector::operator * (const TMatrix& arg) const{
    TVector V(n);
#ifdef _DEBUG
    if (n != arg.n)  //сделаем вид что вектор по умолчанию - строка
        throw 1;
#endif
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++){
            V[i] = V[i] + data[j]*arg(j, i);
        }
    return V;
}

// Оператор векторного умножения векторов
TVector TVector::operator ^ (const TVector& arg) const {
#ifdef _DEBUG
    if ((n != arg.n) and (n != 3))
        throw 1;
#endif
    TVector V(n);
    V[0]= 0; V[1] = 0; V[2] = 0;
    V[0] = data[1]*arg[2] - data[2]*arg[1];
    V[1] = data[2]*arg[0] - data[0]*arg[2];
    V[2] = data[0]*arg[1] - data[1]*arg[0];
    return V;
}

//Дружественная функция - оператор умножения числа на вектор
TVector operator * (double lvalue, const TVector& rvalue) {
    return rvalue * lvalue;
}

// Функция получения модуля вектора
double TVector::length() const {
    double R;
    for (int i = 0; i < n; i++ ) {
        R = R + data[i]*data[i];
    }
    R = sqrt(R);
    return R;
}

// Функция нормирования вектора
TVector& TVector::norm() {
    double r = this->length();
    if (r > 0)
        for (int i = 0; i < n; i++)
            data[i] /= r;  // делим каждый "свой" элемент на модуль
    return *this; // возвращаем "себя". Т.к. this - это указатель, то для получения ссылки (т.е. типа TVector&) нужно перед ним поставить *.
}

// Поворот вектора вокруг заданной оси на заданный угол при помощи формулы Родрига
TVector TVector::rotateByRodrigFormula(double phi, const TVector& axis) const {
    TVector V = *this;
    TVector VR (n), d(n), q(n), b(n);
    double q0;
    d = axis;
    d.norm();     //d - отнормированный axis
    // VR = V*cos(phi) + d*(1-cos(phi))*d*V + d^V*sin(phi);
    q = d*sin(phi/2);
    b = q^V;
    q0 = cos(phi/2);
    VR = V + 2*q0*b + 2*(q^b);
    return VR;
}

// МАТРИЦЫ

// Конструктор по умолчанию
TMatrix::TMatrix() : n(0), m(0), data(NULL) {}

// Конструктор с заданной размерностью
TMatrix::TMatrix(int n, int m) : n(0), m(0), data(NULL) { resize(n, m); }

// Конструктор копий
TMatrix::TMatrix(const TMatrix& rvalue) : n(0), m(0), data(NULL) {
    (*this) = rvalue;
}

// Оператор присваивания
TMatrix& TMatrix::operator = (const TMatrix& rvalue) {
    // Если левый операнд не совпадает с правым
    if (this != &rvalue) {
        // Удаление ранее выделенной памяти
        this->~TMatrix();
        // Выделение новой памяти по размерам правого операнда
        resize(rvalue.n, rvalue.m);
        // Перенос данных из правого операнда в левый построчно
        for (int i = 0; i < n; i++)
            memcpy(data[i], rvalue.data[i], sizeof(double)*m);
    }
    // Возврат ссылки на левый операнд для возможности цепочки присваиваний
    return (*this);
}

// Деструктор объекта матрицы
TMatrix::~TMatrix() {
    if (data) {
        for (int i = 0; i < n; i++)
            delete[] data[i];
        delete[] data;
        data = NULL;
        n = m = 0;
    }
}

// Функция задания размерности матрицы
void TMatrix::resize(int n, int m) {
    // Кол-во строк, которые нужно перенести в новые блоки данных
    int min_n = this->n < n ? this->n : n;
    // Если кол-во столбцов не совпадает
    if (this->m != m) {
        // Кол-во столбцов, которые нужно перенести в новые блоки данных
        int min_m = this->m < m ? this->m : m;
        // Цикл построчного переноса данных в новые блоки
        for (int i = 0; i < min_n; i++) {
            // Создание нового блока-строки
            double *newDataRow = new double[m];
            // Перенос данных в новый блок-строку
            memcpy(newDataRow, data[i], sizeof(double)*min_m);
            // Удаление старого блока строки на этом месте
            delete[] data[i];
            // Прикрепление нового блока-строки на старое место
            data[i] = newDataRow;
        }
        // Сохранение нового размера
        this->m = m;
    }
    // Если кол-во строк не совпадает
    if (this->n != n) {
        // Создание нового блока-контейнера
        double **newData = new double*[n];
        // Перенос содержимого старого контейнера в новый
        memcpy(newData, data, sizeof(double*)*min_n);
        // Удаление лишних строк из старого контейнера
        for (int i = n; i < this->n; i++) { delete[] data[i]; }
        // Удаление старого контейнера
        if (data) { delete[] data; }
        // Создание недостающих строк в новом контейнере
        for (int i = this->n; i < n; i++) { newData[i] = new double[m]; }
        // Привязка старого контейнера к новому
        data = newData;
        this->n = n;
    }
}

// Оператор - унарный минус
TMatrix TMatrix::operator - () const {
    TMatrix M(n,m);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            M(i,j) = -data[i][j];
    return M;
}

// Оператор вычитания матриц
TMatrix TMatrix::operator - (const TMatrix& arg) const{
#ifdef _DEBUG
    if ((n != arg.n) || (m != arg.m))
        throw 1;
#endif
    TMatrix M(n, m);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            M(i,j) = data[i][j] - arg(i,j);
    return M;
}

// Оператор сложения матриц
TMatrix TMatrix::operator + (const TMatrix& arg) const {
#ifdef _DEBUG
    if ((n != arg.n) || (m != arg.m))
        throw 1;
#endif
    TMatrix M(n, m);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            M(i,j) = data[i][j] + arg(i,j);
    return M;
}

// Оператор умножения матрицы на число
TMatrix TMatrix::operator * (double arg) const{
    TMatrix M(n, m);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            M(i,j) = data[i][j]*arg;
    return M;
}

// Оператор умножения матриц
TMatrix TMatrix::operator * (const TMatrix& arg) const{
#ifdef _DEBUG
    if (m != arg.n)
        throw 1;
#endif
    int n1 = n;
    int m1 = m;
    int n2 = arg.colHigh();             //строго m1 = n2; не исключено m2 != n1
    int m2 = arg.colCount();
    TMatrix M(n1, m2);
    for (int i = 0; i < n1; i++)
        for (int j = 0; j < m2; j++){
            M(i, j) = 0;
            for (int k = 0; k < m1; k++)
                M(i, j) += data[i][k]*arg(k,j);
        }

    return M;
}

// Оператор умножения матрицы на вектор
TVector TMatrix::operator * (const TVector& arg) const{
#ifdef _DEBUG
    if (m != arg.n)   // предположим что вектор - столбец
        throw 1;
#endif
    TVector V(n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++) {
            V[i] = V[i] + arg[j]*data[i][j];
        }
    return V;
}
TMatrix operator * (double lvalue, const TMatrix& rvalue) {
    return rvalue * lvalue;
}

//функция приведения матрицы к треугольному виду
TMatrix TMatrix::triangle() const {
    int counter = 0;
    TMatrix M = *this;

    for (int i = 0; i < n-1; i++)
        for (int k = i+1; k < n; k++)
            for (int j = n-1; j >= i; j--) {
                if (M(i, i) != 0)
                    M(k, j) = M(k, j) - (M(k, i)/M(i, i))*M(i, j); //в общем случае
                else
                    for (int h = i+1; h < n; h++) {  //если M[i][i] = 0 ищем != в i-ом столбце и меняем строки с i-ой местами
                        if ( M(h, i) != 0 ) {
                            int nz = h; //индекс строки с ненулевым элементом
                            for (int ni = 0; ni < n; ni++) {
                                double reserve = M(i, ni);
                                M(i, ni) = M(nz, ni);
                                M(nz, ni) = reserve;
                                counter++;
                            }
                            M(k, j) = M(k, j) - (M(k, i)/M(i, i))*M(i, j);
                        }
                    }
            }
    return M;
}

// Функция вычисления детерминанта
double TMatrix::det() const {
    double D = 1488;                    //вычсляем методом Гаусса-Жордана
    TMatrix M = *this;                 //для этого приводим к треугольному виду
    int counter = 0;
    for (int i = 0; i < n-1; i++)
        for (int k = i+1; k < n; k++)
            for (int j = n-1; j >= i; j--) {
                if (M(i, i) != 0)
                    M(k, j) = M(k, j) - (M(k, i)/M(i, i))*M(i, j); //в общем случае
                else
                    for (int h = i+1; h < n; h++) { //если M[i][i] = 0 ищем != в i-ом столбце и меняем строки с i-ой местами
                        int col = 0;
                        if ( M(h, i) != 0 ) {
                            int nz = h; //индекс строки с ненулевым элементом
                            for (int ni = 0; ni < n; ni++) {
                                double reserve = M(i, ni);
                                M(i, ni) = M(nz, ni);
                                M(nz, ni) = reserve;
                                counter++;
                                col++;
                            }
                            M(k, j) = M(k, j) - (M(k, i)/M(i, i))*M(i, j);
                        }
                        else if (col == 0)
                            D = 0; //определитель равен нулю если найти строку с ненулевым значением не удалось
                    }
            }
    //непосредственно вычисление определителя
    if (D != 0) {
        D = 1;
        for (int i = 0; i < n; i++)
            D = D*pow(-1, counter)*M(i,i);
    }
    return D;
}

// Производящая функция для формирования единичной матрицы
TMatrix TMatrix::E(int n) {
    TMatrix E(n,n);
    for (int i = 0; i < n; i++) {
        E(i,i) = 1;
        for (int j = i+1; j < n; j++) { E(i,j) = E(j,i) = 0; }
    }
    return E;
}

// Оператор обращения матриц (метод Гаусса)
TMatrix TMatrix::operator ! () const {
#ifdef _DEBUG
    if (m != arg.n)
        throw 1;
#endif
    TMatrix A(*this), X( E(n) ); //A - исходная; X - результирующая, но пока что единичная
    for (int i = 0; i < n; i++) {
        double rez = A(i, i); //временная переменная для ведущего элемента
        if ( abs(rez) <= eps ) {        //проверка на вшивость A(i, i) != 0 при чем машинному
            double lrez, oldrez = rez;
            for (int l = i+1; l < n; l++) {    //поиск наибольшего элемента в i-ом столбце в строках [i+1; n)
                if (abs(A(l, i)) > rez) {
                    rez = A(l, i);
                    lrez = l;
                }
            }
            if (rez == oldrez) {        //проверка нашли ли наибольший элемент
                cout << "the determinant equals 0" << endl;
                throw 1;
            }
            for (int l = 0; l < n; l++) { //меняем i-ую и lrez-ую строки местами в обеих матрицах
                double arez = A(lrez, l), xrez = X(lrez, l);
                A(lrez, l) = A(i, l);
                A(i, l) = arez;
                X(lrez, l) = X(i, l);
                X(i, l) = xrez;
            }

        }
        for (int j = 0; j < n; j++) {
            A(i, j) /= rez;      //непосредственно сам процесс нормировки на A(i, i)
            X(i, j) /= rez;
        }

        for (int k = 0; k < n; k++) {
            double krez = A(k, i);  //временная переменная для A(k, i)
            for (int j = 0; j < n; j++)
                if (k != i) {
                    A(k, j) -= krez*A(i, j);
                    X(k, j) -= krez*X(i, j);
                }

        }

    }
    return X;
}

// Функция транспонирования
TMatrix TMatrix::t() const {
    if (n == m) {
        TMatrix tp = *this, a = *this;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++) {
                tp(i, j) = a(j, i);
            }
        return tp;
    }
    else {
        TMatrix tp (m, n);
        TMatrix a = *this;
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++) {
                tp(i, j) = a(j, i);
            }
        return tp;
    }
}

// Функция перестановки столбцов
TMatrix& TMatrix::swapColumns(int i, int j) {
    TMatrix a = *this;
    for (int k = 0; k < m; k++) {
        data[k][j] = a (k, i);
        data[k][i] = a (k, j);
    }
    return *this;
}

// Функция перестановки строк
TMatrix& TMatrix::swapRows(int i, int j) {
    TMatrix a = *this;
    for (int k = 0; k < n; k++) {
        data[i][k] = a (j, k);
        data[j][k] = a (i, k);
    }
    return *this;
}

}
