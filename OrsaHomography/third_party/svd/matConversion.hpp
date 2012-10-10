//Copyright (C) 2011 Pierre Moulon
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

#ifndef LIBS_SVD_CONVERSION_H_
#define LIBS_SVD_CONVERSION_H_

#include <Eigen/Core>
#include <Eigen/SVD>
#include "extras/libNumerics/matrix.h"

namespace EigenWrapper {

typedef libNumerics::matrix<double> Mat;
typedef Eigen::MatrixXd MatEigen;

inline void EigenToMat(const MatEigen & matEigen, Mat & mat)
{
  mat = Mat(matEigen.rows(), matEigen.cols());
  for (int j=0; j<matEigen.rows();++j)
  {
    for (int i=0; i<matEigen.cols();++i)
      mat(j,i) = matEigen(j,i);
  }  
}

inline void MatToEigen(const Mat & mat, MatEigen & matEigen)
{
  matEigen.resize(mat.nrow(), mat.ncol());
  for (int j=0; j<mat.nrow();++j)
  {
    for (int i=0; i<mat.ncol();++i)
      matEigen(j,i) = mat(j,i);
  }  
}

} // namespace EigenWrapper


#endif // LIBS_SVD_CONVERSION_H_
