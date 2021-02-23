#include <iostream>
#include "linalg.h"
#include <cstring>
#include <math.h>
#define pi 3.141592653589793
using namespace std;

int main() {
    using namespace MyLinearAlgebra;
    cout << endl;
    TMatrix  B (2, 2);
    //A = B;


    B(0, 0) = 1361646033;
    B(0, 1) = 587867.64;
    B(1, 0) = 587867.64;                  //проверка вычисления определителя
    B(1, 1) = 254;

    cout << "Matrix B = " << endl;
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            cout << B(i, j) << " ";
        }
        cout << endl;
    }
    cout << endl;

    /*    TMatrix Bt = B.t();
    cout << "Matrix B^T = " << endl;
       for (int i = 0; i < 2; i++) {
           for (int j = 0; j < 2; j++) {
              cout << Bt(i, j) << " ";
           }
           cout << endl;
         }
    cout << endl;
    */

    TMatrix Bobr = !B;
    double Det = B.det();
    cout << "determinant = " << Det << endl << endl;

    cout << "Matrix !B = " << endl;
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            cout << Bobr(i, j) << " ";
        }
        cout << endl;
    }
    cout << endl;


}
