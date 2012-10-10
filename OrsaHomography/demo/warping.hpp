//Copyright (C) 2011 Lionel Moisan, Pascal Monasse, Pierre Moulon
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

#ifndef WARPING_H
#define WARPING_H

#include <limits>
#include "extras/libNumerics/numerics.h"
#include "libImage/sample.hpp"
#include "demo/Rect.hpp"

/// Apply homography transform.
/// Indicate if \a H is orientation preserving around the point.
bool TransformH(const libNumerics::matrix<double> &H, double &x, double &y)
{
  libNumerics::vector<double> X(3);
  X(0)=x; X(1)=y; X(2)=1.0;
  X = H*X;
  bool positive = (X(2)*H(2,2)>0);
  X /= X(2);
  x = X(0); y = X(1);
  return positive;
}

// Compute the common area of warped by homography image1 and image2.
bool IntersectionBox(int w1, int h1, int w2, int h2,
                     const libNumerics::matrix<double>& H, Rect &inter)
{
  int xCoord[4] = {0, w1-1, w1-1,    0};
  int yCoord[4] = {0,    0, h1-1, h1-1};

  Rect rect1(std::numeric_limits<int>::max(),
             std::numeric_limits<int>::max(),
             std::numeric_limits<int>::min(),
             std::numeric_limits<int>::min());
  for(int i=0; i<4; ++i)
  {
    double xT=xCoord[i], yT=yCoord[i];
    TransformH(H, xT, yT);
    rect1.growTo(xT,yT);
  }

  Rect rect2(0,0,w2-1,h2-1);
  return rect2.intersect(rect1, inter);
}

/// Warp an image given a homography and parameter :
/// Parameter tx,ty are provided by Overlap_ComputeIntersectionBox(...).
template <class Image>
void warpImageTo(const Image &imageIn, Image &imageOut,
                 libNumerics::matrix<double> H, int tx, int ty)
{
  const int wOut=imageOut.Width(), hOut=imageOut.Height();
  H=H.inv();

#ifdef _OPENMP
    #pragma omp parallel for
#endif
  for(int j=0; j<hOut; ++j)
  {
    for(int i=0; i<wOut; ++i)
    {
      double xT=i-tx, yT=j-ty;
      if(TransformH(H, xT, yT) && imageIn.Contains(yT,xT))
          imageOut(j,i) = libs::SampleLinear(imageIn, yT, xT);
    }
  }
}

/// Warp imageA over imageB with homography.
/// Use backward mapping with bilinear sampling.
///- Algo :
/// For destination pixel search which pixels from ImageA and ImageB contribute.
/// Perform a mean blending in the overlap zone, transfered original
/// value in the other part.
template<class Image>
void Warp(const Image &imageA, // To be warped
          const Image &imageB, // Reference
          Image &warpImage,    // Output image
          libNumerics::matrix<double> H,
          int tx, int ty)
{
  const int wOut = warpImage.Width();
  const int hOut = warpImage.Height();
  H = H.inv();

#ifdef _OPENMP
  #pragma omp parallel for
#endif
  for(int j=0; j < hOut; ++j)
  {
    for(int i = 0; i < wOut; ++i)
    {
      const int xPos=i-tx, yPos=j-ty;
      double xA=xPos, yA=yPos;
      bool bAContrib = TransformH(H, xA, yA) && imageA.Contains(yA,xA);
      bool bBContrib = imageB.Contains(yPos,xPos);

      if(bAContrib && bBContrib) //mean blending between ImageA and ImageB
        warpImage(j,i) = typename Image::Tpixel(
          libs::SampleLinear(imageA,yA,xA)*.5f + imageB(yPos,xPos)*.5f);
      else if(bAContrib) //only ImageA
        warpImage(j,i) = libs::SampleLinear(imageA,yA,xA);
      else if(bBContrib) //only ImageB
        warpImage(j,i) = imageB(yPos,xPos);
    }
  }
}

#endif // WARPING_H
