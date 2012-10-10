//Copyright (C) 2010 Pascal Monasse
//
//This program is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//This program is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "numerics.h"
#include <cmath>
#include <vector>
#include <limits>
#include <algorithm>

namespace libNumerics {

inline flnum ABS(flnum x)
{ return (x >= 0)? x: -x; }

/// Resolution by LU decomposition with pivot.
bool solveLU(const matrix<flnum>& A, const vector<flnum>& B, vector<flnum>& X)
{
    X = B;
    return solveLU(A, X);
}

/// Replace X by A^{-1}X, by LU solver.
bool solveLU(matrix<flnum> A, vector<flnum>& X)
{
    assert(A.nrow() == A.ncol());
    int	n = A.nrow();
    vector<flnum> rowscale(n); // Implicit scaling of each row
    std::vector<int> permut(n,0); // Permutation of rows

    // Get the implicit scaling information of each row
    for(int i=0; i< n; i++) {
        flnum max = 0.0;
        for(int j=0; j< n; j++) {
            flnum tmp = ABS(A(i,j));
            if (tmp> max)
                max = tmp;
        }
        if(max == 0.0)
            return false;
        rowscale(i) = 1.0/max;
    }
    // Perform the decomposition
    for(int k=0; k < n; k++) {
        // Search for largest pivot element
        flnum max = rowscale(k)*ABS(A(k,k));
        int imax = k;
        for(int i=k+1; i < n; i++) {
            flnum tmp = rowscale(i)*ABS(A(i,k));
            if(tmp > max) {
                max = tmp;
                imax = i;
            }
        }
        if(max == 0.0)
            return false;

        // Interchange rows if needed
        if(k != imax) {
            A.swapRows(k, imax);
            rowscale(imax) = rowscale(k); // Scale of row k no longer needed
        }
        permut[k] = imax; // permut(k) was not initialized before
        flnum Akk = 1/A(k,k);
        for(int i=k+1; i < n; i++) {
            flnum tmp = A(i,k) *= Akk; // Divide by pivot
            for (int j=k+1;j < n; j++) // Reduce the row
                A(i,j) -= tmp*A(k,j);
        }
    }
    // Forward substitution
    for (int k=0; k < n; k++) {
        flnum sum = X(permut[k]);
        X(permut[k]) = X(k);
        for(int j = 0; j < k; j++)
            sum -= A(k,j)*X(j);
        X(k) = sum;
    }
    // Backward substitution
    for(int k=n-1; k >= 0; k--) {
        flnum sum = X(k);
        for(int j=k+1; j < n; j++)
            sum -= A(k,j)*X(j);
        X(k) = sum/A(k,k);
    }
    return true;
}

} // namespace libNumerics
