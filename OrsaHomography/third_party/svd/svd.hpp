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

#ifndef LIBS_SVD_H_
#define LIBS_SVD_H_

#include <Eigen/Core>
#include <Eigen/SVD>
#include "extras/libNumerics/matrix.h"
#include "third_party/svd/matConversion.hpp"
#include "third_party/svd/eigenWrapper.hpp"

namespace SVDWrapper {

  typedef libNumerics::matrix<double> Mat;
  typedef libNumerics::vector<double> Vec;
  typedef Eigen::MatrixXd MatEigen;

  /// Solve the linear system Ax = 0 with ||x|| = 1.0 via SVD.
  /// Return true if the ratio of singular values SV(0)/SV(N-2) < dRatio.
  /// The return value indicates whether the solution is unique.
  static bool Nullspace(const Mat & A, Vec *nullspace, double dRatio=1e-5)
  {
    typedef Eigen::MatrixXd MatEigen;
    MatEigen AEigen, null;
    EigenWrapper::MatToEigen(A, AEigen);
    bool bOk = EigenWrapper::Nullspace(&AEigen, &null, dRatio);
    (*nullspace).read( null.data() );
    return bOk;
  }

  /// Inverse of norm-2 condition value (ratio of extreme singular values)
  inline double InvCond(const Mat& A) {
    Eigen::MatrixXd AEigen, sing;
    EigenWrapper::MatToEigen(A, AEigen);
    EigenWrapper::SingularValues(&AEigen, &sing);
    return sing(A.nrow()-1)/sing(0);
  }

  /// Make rank<=2.
  inline void EnforceRank2_3x3(const Mat& A, Mat *ARank)
  {
    Eigen::MatrixXd MEigen, MRank2Eigen;
    EigenWrapper::MatToEigen(A, MEigen);
    EigenWrapper::EnforceRank2_3x3(MEigen, &MRank2Eigen);
    EigenWrapper::EigenToMat(MRank2Eigen, *ARank);
  }

  /// Save the two last nullspace vector as 3x3 matrices.
  /// It uses Eigen to compute the SVD decomposition.
  inline void Nullspace2_Remap33(const Mat &A, Mat& f1, Mat& f2) {
      using Eigen::Map;
      using namespace EigenWrapper;

      MatEigen f1E, f2E;
      MatEigen AEigen;
      MatToEigen(A, AEigen);
      Nullspace2(&AEigen, &f1E, &f2E);

      typedef Eigen::Matrix<double, 3, 3> Mat3;
      typedef Eigen::Matrix<double, 3, 3, Eigen::RowMajor> RMat3;
      // Update f1 and f2
      EigenToMat(Map<RMat3>(f1E.data()), f1);
      EigenToMat(Map<RMat3>(f2E.data()), f2);
  }

} // namespace SVDWrapper


#endif // LIBS_SVD_H_
