#ifndef LINALG_H
#define LINALG_H

#include <vector>
#include <stdio.h>

namespace MyLinearAlgebra {
// Опережающие объявления
class TMatrix;

// Объявление класса векторов
class TVector {
protected:
    int n;          // Размерность вектора
    double *data;  // Элементы вектора
public:
    TVector();  // Конструктор по умолчанию
    TVector(int n); // Конструктор с заданным кол-вом элементов
    TVector(const TVector& rvalue); // Конструктор копий
    TVector& operator = (const TVector& rvalue); // Оператор присваивания
    virtual ~TVector(); // Деструктор
    inline int size() const { return n; } // Функция получение кол-ва элементов вектора
    inline int high() const { return n - 1; } // Функция получения индекса последнего элемента
    void resize(int n); // Функция задания кол-ва элементов вектора
    inline double& operator[](int i) { return data[i]; } // Оператор доступа к элементам вектора
    inline const double& operator[](int i) const { return data[i]; } // Оператор константного доступа к элементам вектора
    TVector operator - () const; // Оператор - унарный минус
    TVector operator - (const TVector& arg) const; // Оператор вычитания векторов
    TVector operator + (const TVector& arg) const; // Оператор сложения векторов
    TVector operator * (double arg) const; // Оператор умножения вектора на число
    double operator * (const TVector& arg) const; // Оператор скалярного умножения векторов
    TVector operator * (const TMatrix& arg) const; // Оператор умножения вектора на матрицу
    TVector operator ^ (const TVector& arg) const; // Оператор векторного умножения векторов
    friend TVector operator * (double lvalue, const TVector& rvalue); // Дружественная функция - оператор умножения числа на вектор
    double length() const; // Функция получения модуля вектора
    TVector& norm(); // Функция нормирования вектора
    TVector rotateByRodrigFormula(double phi, const TVector& axis) const; // Поворот вектора вокруг заданной оси на заданный угол при помощи формулы Родрига
};

// Объявление класса матриц
class TMatrix {
protected:
    int n, m; // Размерность матрицы (число строк и столбцов)
    double **data; // Элементы матрицы
public:
    TMatrix(); // Конструктор по умолчанию
    TMatrix(int n, int m); // Конструктор с заданной размерностью
    TMatrix(const TMatrix& rvalue); // Конструктор копий
    TMatrix& operator = (const TMatrix& rvalue); // Оператор присваивания
    virtual ~TMatrix(); // Деструктор
    inline int rowCount() const { return n; } // Функция получения количества строк
    inline int colCount() const { return m; } // Функция получения кол-ва столбцов
    inline int rowHigh() const { return n-1; } // Функция получения индекса последней строки
    inline int colHigh() const { return m-1; } // Функция получения индекса последнего столбца
    void resize(int n, int m); // Функция задания размерности
    inline double& operator()(int i, int j) { return data[i][j]; } // Оператор доступа к элементам матрицы
    inline const double& operator()(int i, int j) const { return data[i][j]; } // Оператор константного доступа к элементам матрицы
    TMatrix operator - () const; // Оператор - унарный минус
    TMatrix operator - (const TMatrix& arg) const; // Оператор вычитания матриц
    TMatrix operator + (const TMatrix& arg) const; // Оператор сложения матриц
    TMatrix operator * (double arg) const; // Оператор умножения матрицы на число
    TMatrix operator * (const TMatrix& arg) const; // Оператор умножения матриц
    TVector operator * (const TVector& arg) const; // Оператор умножения матрицы на вектор
    friend TMatrix operator * (double lvalue, const TMatrix& rvalue); // Дружественная функция - оператор умножения числа на матрицу
    /*l*/   TMatrix operator ! () const; // Оператор обращения матриц (метод Гаусса)
    TMatrix triangle() const; //функция приведения матрицы к треугольному виду
    double det() const; // Функция вычисления детерминанта
    static TMatrix E(int n); // Функция формирования единичной матрицы
    TMatrix t() const; // Функция транспонирования
    TMatrix& swapRows(int i, int j); // Функция перестановки строк
    TMatrix& swapColumns(int i, int j); // Функция перестановки столбцов
};
}

#endif // LINALG_H
