базовая линейная алгебра
========================
данный пакет предоставляет возможность выполнения основных операций с векторами и матрицами.

Программный модуль представлен заголовочным файлом и файлом реализации. 
Заголовочный файл в рамках пространства имён MyLinearAlgebra содержит объявления двух классов: TVector и TMatrix.
Далее для каждого из классов приведен перечень методов с указанием назначения и примером использования.

TVector
-------

1) TVector(); - конструктор по умолчанию
Использование:​ TVector A;

2) TVector(int n); - конструктор с заданным кол-вом элементов
Использование:​ TVector B(3);

3) TVector(const TVector& rvalue); - конструктор копий
Использование:​ TVector C(B);

4) TVector& operator = (const TVector& rvalue); - оператор присваивания
Использование:​ C = B;

5) inline int size() const { return n; } - Функция получение кол-ва элементов вектора
Использование:​ B.size;

6) inline int high() const { return n - 1; } - Функция получения индекса последнего элемента
Использование:​ B.high();

7) void resize(int n); - Функция задания кол-ва элементов вектора
Использование:​ B.resize(2);

8) inline double& operator[](int i) { return data[i]; } - Оператор доступа к элементам вектора
Использование:​ B[0];

9) inline const double& operator[](int i) const { return data[i]; } - Оператор константного доступа к элементам вектора
Использование:​ B[0];

10) TVector operator - () const; - Оператор - унарный минус
Использование: ​ C = -B;

11) TVector operator - (const TVector& arg) const; - Оператор вычитания векторов
Использование: ​ A = B - C;

12) TVector operator + (const TVector& arg) const; - Оператор сложения векторов
Использование: ​ A = A = B + C;

13) TVector operator * (double arg) const; - Оператор умножения вектора на число
Использование:​ C = B*2;

14) friend TVector operator * (double lvalue, const TVector& rvalue); - Дружественная функция - оператор умножения числа на вектор
Использование:​ C = 0.5*C;

15) double operator * (const TVector& arg) const; - Оператор скалярного умножения векторов
Использование:​ double sm = C*B;

16) TVector operator ^ (const TVector& arg) const; - Оператор векторного умножения векторов
Использование:​ A = C^B;

17) TVector operator * (const TMatrix& arg) const; - Оператор умножения вектора на матрицу
Использование:​ TVector C = A*B;

18) double length() const; - Функция получения модуля вектора
Использование: ​ A.length();

19)TVector& norm(); - Функция нормирования вектора
Использование:​ A.norm();

20)TVector rotateByRodrigFormula(double phi, const TVector& axis) const; - Поворот вектора вокруг заданной оси на заданный угол при помощи формулы Родрига
Использование: ​ TVector S = A.rotateByRodrigFormula(-pi/2, e);

TMatrix
-------

1) TMatrix(); - Конструктор по умолчанию
Использование: ​ TMatrix A;

2) TMatrix(int n, int m); - Конструктор с заданной размерностью
Использование: ​ TMatrix B(3, 3), C(3, 3);

3) TMatrix(const TMatrix& rvalue); - Конструктор копий
Использование: ​ TMatrix D(C);

4) TMatrix& operator = (const TMatrix& rvalue); - Оператор присваивания
Использование: ​ A = B;

5) inline int rowCount() const { return n; } - Функция получения количества строк
Использование: ​ B.rowCount();

6) inline int colCount() const { return m; } - Функция получения кол-ва столбцов
Использование: ​ B.colCount();

7) inline int rowHigh() const { return n-1; } - Функция получения индекса последней строки
Использование: ​ B.rowHigh();

8) inline int colHigh() const { return m-1; } - Функция получения индекса последнего столбца
Использование: ​ B.colHigh();

9) void resize(int n, int m); - Функция задания размерности
Использование: ​ B.resize(2, 2);

10) inline double& operator()(int i, int j) { return data[i][j]; } // Оператор доступа к элементам матрицы
Использование: ​ B(0, 0);

11) inline const double& operator()(int i, int j) const { return data[i][j]; } // Оператор константного доступа к элементам матрицы
Использование: ​ B(0, 0);

12) TMatrix operator - () const; - Оператор - унарный минус
Использование: ​ A = -B;

13) TMatrix operator - (const TMatrix& arg) const; - Оператор вычитания матриц
Использование: ​ C = A - B;

14) TMatrix operator + (const TMatrix& arg) const; - Оператор сложения матриц
Использование: ​ C = A + B;

15) TMatrix operator * (double arg) const; - Оператор умножения матрицы на числоИспользование: ​ A = B*2;

16) TMatrix operator * (const TMatrix& arg) const; - Оператор умножения матриц
Использование: ​ C = A*B;

17) TVector operator * (const TVector& arg) const; - Оператор умножения матрицы на вектор
Использование: ​ V1 = B*V;

18) friend TMatrix operator * (double lvalue, const TMatrix& rvalue); - Дружественная функция - оператор умножения числа на матрицу
Использование: ​ A = 2*B;

19) TMatrix operator ! () const; - Оператор обращения матриц (метод Гаусса)
Использование: ​ Bobr = !B;

20) double det() const; - Функция вычисления детерминанта (определителя) матрицы
Использование: ​ d = B.det;

21) TMatrix t() const; - Функция транспонирования
Использование: ​ B.t();

22) TMatrix& swapRows(int i, int j); - Функция перестановки строк
Использование: ​ B.swapRows(1, 2);

23) TMatrix& swapColumns(int i, int j); - Функция перестановки столбцов
Использование: ​ B.swapColumns(1, 2);

Связь
-------
Если у вас остались вопросы относительно пакета или вы нашли ошибку, пожалуйста, напишите мне:

akulakova29@gmail.com
 
-----
ftn21
https://github.com/ftn21

base linear algebra
===================
this package provides an ability to use basic operations with vectors and matrices.

The program module is represented by a header file and an implementation file.
The header file within the MyLinearAlgebra namespace contains declarations of two classes: TVector and TMatrix.
Bellow you can find a list of methods for each of the classes with a description and an example of usage.

TVector
-------

1) TVector (); - default constructor
Usage:​ TVector A;

2) TVector(int n); - constructor with a specified number of elements.
Usage:​ TVector B(3);

3) TVector(const TVector& rvalue); - a copy constructor
Usage:​ TVector C(B);

4) TVector& operator = (const TVector& rvalue); - assignment operator
Usage:​ C = B;

5) inline int size() const { return n; } Function obtaining the number of elements of a vector
Usage: B. size;

6) inline int high () const { return n-1; } - Function for getting the index of the last element
Usage: B. high();

7) void resize(int n); - Function for setting the number of vector elements
Usage: B. resize(2);

8) inline double& operator[](int i) { return data[i]; } - Operator for accessing vector elements
Usage:​ B[0];

9) inline const double& operator[](int i) const { return data[i]; } - Operator for constant access to vector elements
Usage:​ B[0];

10) TVector operator - () const; - Operator-unary minus
Usage: C = -B;

11) TVector operator - (const TVector& arg) const; - Vector subtraction operator
Usage: A = B - C;

12) TVector operator + (const TVector& arg) const; - Vector addition operator
Usage: A = A = B + C;

13) TVector operator * (double arg) const; - Operator for multiplying a vector by a number
Usage:​ C = B*2;

14) friend TVector operator * (double lvalue, const TVector& rvalue); - Friendly function - operator for multiplying a number by a vector
Usage:​ C = 0.5*C;

15) double operator * (const TVector& arg) const; - Operator for scalar multiplication of vectors
Usage:​ double sm = C*B;

16) TVector operator ^ (const TVector& arg) const; - Vector multiplication operator of vectors
Usage:​ A = C^B;

17) TVector operator * (const TMatrix& arg) const; - Operator for multiplying a vector by a matrix
Usage:​ TVector C = A*B;

18) double length () const; - Function for getting the vector modulus
Usage: A. length();

19) TVector& norm (); - Vector normalization function
Usage: A. norm();

20) TVector rotateByRodrigFormula(double phi, const TVector& axis) const; - Rotation of the vector around a given axis by a given angle using the Rodrigue formula
Usage: TVector S = A. rotateByRodrigFormula(-pi/2, e);

TMatrix
-------

1) TMatrix (); - Default constructor
Usage: TMatrix A;

2) TMatrix(int n, int m); - Constructor with the specified dimension
Usage: TMatrix B (3, 3), C(3, 3);

3) TMatrix(const TMatrix& rvalue); - Copy Constructor
Usage: TMatrix D(C);

4) TMatrix& operator = (const TMatrix& rvalue); - The assignment operator
Usage: A = B;

5) inline int rowCount () const { return n; } - Function for getting the number of rows
Usage: B. rowCount();

6) inline int colCount() const { return m; } - Function for getting the number of columns
Usage: B. colCount();

7) inline int rowHigh () const { return n-1; } - Function for getting the index of the last row
Usage: B. rowHigh();

8) inline int colHigh () const { return m-1; } - Function for getting the index of the last column
Usage: B. colHigh();

9) void resize(int n, int m); - the Function of the task dimension
Usage: B. resize(2, 2);

10) inline double& operator () (int i, int j) { return data[i][j]; } / / Matrix element access operator
Usage: B(0, 0);

11) inline const double& operator () (int i, int j) const { return data[i][j]; } / / Operator of constant access to matrix elements
Usage: B(0, 0);

12) TMatrix operator - () const; - Operator-unary minus
Usage: A = -B;

13) TMatrix operator - (const TMatrix& arg) const; - Matrix subtraction operator
Usage: C = A - B;

14) TMatrix operator + (const TMatrix& arg) const; - Matrix addition operator
Usage: C = A + B;

15) TMatrix operator * (double arg) const; - Matrix multiplication operator by numberuse: A = B*2;

16) TMatrix operator * (const TMatrix& arg) const; - Matrix multiplication operator
Usage: C = A*B;

17) TVector operator * (const TVector& arg) const; - Operator for multiplying a matrix by a vector
Usage: V1 = B*V;

18) friend TMatrix operator * (double lvalue, const TMatrix& rvalue); - Friendly function - operator for multiplying a number by a matrix
Usage: A = 2*B;

19) TMatrix operator ! () const; - Matrix inversion operator (Gauss method)
Usage: Bobr = !B;

20) double det () const; - Function for calculating the determinant (determinant) of the matrix
Usage: d = B. det;

21) TMatrix t () const; - Transpose function
Usage: B. t();

22) TMatrix& swapRows(int i, int j); - String permutation function
Usage: B. swapRows(1, 2);

23) TMatrix& swapColumns(int i, int j); - Column permutation function
Usage: B. swapColumns(1, 2);

Contact
-------
If you have further questions regarding the package or you've found a mistake, please email me:

akulakova29@gmail.com
 
-----
ftn21
https://github.com/ftn21

