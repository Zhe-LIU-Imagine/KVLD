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

#ifndef HOMOGRAPHY_MODEL_H_
#define HOMOGRAPHY_MODEL_H_

#include <vector>
#include <cmath>
#include "libOrsa/orsa_model.hpp"
#include "extras/libNumerics/matrix.h"

namespace orsa {

/// Homography model used for robust estimation with ORSA algorithm.
class HomographyModel : public OrsaModel {
public:
  HomographyModel(const Mat &x1, int w1, int h1,
                  const Mat &x2, int w2, int h2,
                  bool symError=false);

  /// 4 point correspondences required to compute a homography.
  int SizeSample() const { return  4; }

  /// Only 1 homography can be estimated from a sample of 4 points.
  int NbModels() const { return 1; }

  /// Distance used to distinguish inlier/outlier is to a point
  virtual bool DistToPoint() const { return true; }

  /// Estimated homography satisfies the equation y = H x.
  void Fit(const std::vector<size_t> &indices, std::vector<Mat> *H) const;

  /// Square reprojection error for a given point through the model H.
  double Error(const Mat &H, size_t index, int* side=0) const;

  /// Unnormalize a given model (from normalized to image space).
  void Unnormalize(Model *model) const;

private:
    bool symError_; ///< Use symmetric error or transfer error in image 2?
    bool IsOrientationPreserving(const std::vector<size_t> &indices,
                                 const Mat& H) const;
};

}  // namespace orsa

#endif
