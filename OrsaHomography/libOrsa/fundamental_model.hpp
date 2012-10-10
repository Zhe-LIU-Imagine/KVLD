// Copyright (c) 2007-2011 libmv authors.
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to
// deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.

#ifndef FUNDAMENTAL_MODEL_H_
#define FUNDAMENTAL_MODEL_H_

#include <vector>
#include "libOrsa/orsa_model.hpp"
#include "extras/libNumerics/matrix.h"

namespace orsa {

///  Fundamental 7-point model, used for robust estimation.
///
/// See page 281 of book by Hartley-Zisserman.
/// The equation is \f$det(F_1 + \alpha F_2) = 0\f$.
class FundamentalModel : public OrsaModel {
public:
  FundamentalModel(const Mat &x1, int w1, int h1,
                   const Mat &x2, int w2, int h2,
                   bool symError=false);

  /// 7 points are required to compute a fundamental matrix.
  int SizeSample() const { return  7;}

  /// Up to 3 fundamental matrices are computed from a sample of 7 points.
  int NbModels() const { return 3;}

  /// Distance used to distinguish inlier/outlier is to a line
  virtual bool DistToPoint() const { return false; }

  void Fit(const std::vector<size_t> &indices, std::vector<Mat> *Fs) const;

  /// Square reprojection error for a given point through F.
  double Error(const Mat &F, size_t index, int* side=0) const;
  
  /// Unnormalize a given model (from normalized to image space).
  void Unnormalize(Model *model) const;
private:
    bool symError_;  ///< Use symmetric error or transfer error in image 2?
};

}  // namespace orsa

#endif
