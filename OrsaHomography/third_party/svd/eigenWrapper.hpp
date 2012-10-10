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

#ifndef LIBS_SVD_EIGENWRAPPER_H_
#define LIBS_SVD_EIGENWRAPPER_H_

#include <Eigen/Core>
#include <Eigen/SVD>

namespace EigenWrapper  {

/// Solve the linear system Ax = 0 via SVD. Store the solution in x, such that
/// ||x|| = 1.0. Return true if the ratio of singular value SV(0)/SV(N-2) < ratio.
/// The return value allow to know if many solution are possible
/// Destroys A and resizes x if necessary.
template <typename TMat, typename TVec>
bool Nullspace(TMat *A, TVec *nullspace, double dRatio = 1e-5) {
  if (A->rows() >= A->cols()) {
    Eigen::JacobiSVD<TMat> svd(*A, Eigen::ComputeFullV);
    (*nullspace) = svd.matrixV().col(A->cols()-1);
    double a = svd.singularValues()(A->cols()-2);
    double c = svd.singularValues()(0);
    return !(a < dRatio * c);
  }
  // Extend A with rows of zeros to make it square. It's a hack, but is
  // necessary until Eigen supports SVD with more columns than rows.
  TMat A_extended(A->cols(), A->cols());
  A_extended.block(A->rows(), 0, A->cols() - A->rows(), A->cols()).setZero();
  A_extended.block(0,0, A->rows(), A->cols()) = (*A);
  return Nullspace(&A_extended, nullspace);
}

/// Singular values of square matrix
template <typename TMat, typename TVec>
void SingularValues(TMat *A, TVec *sing) {
  Eigen::JacobiSVD<TMat> svd(*A);
  *sing = svd.singularValues();
}

/// Solve the linear system Ax = 0 via SVD. Finds two solutions, x1 and x2, such
/// that x1 is the best solution and x2 is the next best solution (in the L2
/// norm sense). Store the solution in x1 and x2, such that ||x|| = 1.0. Return
/// the singular value corresponding to the solution x1.  Destroys A and resizes
/// x if necessary.
template <typename TMat, typename TVec1, typename TVec2>
inline double Nullspace2(TMat *A, TVec1 *x1, TVec2 *x2) {
  if (A->rows() >= A->cols()) {
    Eigen::JacobiSVD<TMat> svd(*A,Eigen::ComputeFullV);
    TMat V = svd.matrixV();
    *x1 = V.col(A->cols() - 1);
    *x2 = V.col(A->cols() - 2);
    return svd.singularValues()(A->cols()-1);
  }
  // Extend A with rows of zeros to make it square. It's a hack, but is
  // necessary until Eigen supports SVD with more columns than rows.
  TMat A_extended(A->cols(), A->cols());
  A_extended.block(A->rows(), 0, A->cols() - A->rows(), A->cols()).setZero();
  A_extended.block(0,0, A->rows(), A->cols()) = (*A);
  return Nullspace2(&A_extended, x1, x2);
}

template <typename TMat>
inline void EnforceRank2_3x3(const TMat & A, TMat * out)
{
  if (A.cols() == A.rows() && A.cols() == 3)
  {
    typedef Eigen::Matrix<double, 3, 3> Mat3;
    Eigen::JacobiSVD<Mat3> USV(A, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::VectorXd d = USV.singularValues();
    d[2] = 0.0;
    (*out) = USV.matrixU() * d.asDiagonal() * USV.matrixV().transpose();
  }
}

} // namespace EigenWrapper

#endif  // LIBS_SVD_EIGENWRAPPER_H_
