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

#include <cstdlib>

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

#include "libImage/image_io.hpp"
#include "libImage/image_drawing.hpp"
#include "libImage/image_concat.hpp"

#include "libOrsa/homography_model.hpp"

#include "demo/warping.hpp"

#include "extras/libNumerics/numerics.h"
#include "extras/sift/library.h"
#include "demo/siftMatch.hpp"

#include <time.h>

/// Display average/max error of inliers of homography H.
static void display_stats(const std::vector<Match>& vec_matchings,
                          const std::vector<size_t>& vec_inliers,
                          libNumerics::matrix<double>& H) {
  std::vector<size_t>::const_iterator it=vec_inliers.begin();
  double l2=0, linf=0;
  for(; it!=vec_inliers.end(); ++it) {
    const Match& m=vec_matchings[*it];
    libNumerics::vector<double> x(3);
    x(0)=m.x1;
    x(1)=m.y1;
    x(2)=1.0;
    x = H*x;
    x /= x(2);
    double e = (m.x2-x(0))*(m.x2-x(0)) + (m.y2-x(1))*(m.y2-x(1));
    l2 += e;
    if(linf < e)
      linf = e;
  }
  std::cout << "Average/max error: "
            << sqrt(l2/vec_inliers.size()) << "/"
            << sqrt(linf) <<std::endl;
}

/// ORSA homography estimation
bool ORSA(const std::vector<Match>& vec_matchings, int w1,int h1, int w2,int h2,
          double precision,
          libNumerics::matrix<double>& H, std::vector<size_t>& vec_inliers)
{
  const size_t n = vec_matchings.size();
  if(n < 5)
  {
      std::cerr << "Error: ORSA needs 5 matches or more to proceed" <<std::endl;
      return false;
  }
  libNumerics::matrix<double> xA(2,n), xB(2,n);

  for (size_t i=0; i < n; ++i)
  {
    xA(0,i) = vec_matchings[i].x1;
    xA(1,i) = vec_matchings[i].y1;
    xB(0,i) = vec_matchings[i].x2;
    xB(1,i) = vec_matchings[i].y2;
  }

  orsa::HomographyModel model(xA, w1, h1, xB, w2, h2, true);
  //model.setConvergenceCheck(true);

  if(model.orsa(vec_inliers, 1000, &precision, &H, true)>0.0)
    return false;
  std::cout << "Before refinement: ";
  display_stats(vec_matchings, vec_inliers, H);
  if( model.ComputeModel(vec_inliers,&H) ) // Re-estimate with all inliers
  {
    std::cout << "After  refinement: ";
    display_stats(vec_matchings, vec_inliers, H);
  } else
    std::cerr << "Warning: error in refinement, result is suspect" <<std::endl;
  return true;
}

/// Output inlier and oulier matches in image files.
void display_match(const std::vector<Match>& vec_matchings,
                   std::vector<size_t>& vec_inliers,
                   const libNumerics::matrix<double>* H,
                   int w1,
                   Image<RGBColor>& in, const char* inFile,
                   Image<RGBColor>& out, const char* outFile) {
  std::sort(vec_inliers.begin(), vec_inliers.end());

  // For outliers, show vector (yellow) from prediction to observation
  const RGBColor col=YELLOW;
  std::vector<size_t>::const_iterator it = vec_inliers.begin();
  matchingslist::const_iterator m = vec_matchings.begin();
  if(H) // Otherwise, no prediction
    for(size_t i=0; m != vec_matchings.end(); ++m, ++i) {
      if(it != vec_inliers.end() && i==*it)
        ++it;
      else { //Outlier
        libNumerics::vector<double> v(3);
        v(0) = m->x1; v(1) = m->y1; v(2)=1.0;
        v = (*H)*v;
        v /= v(2);
        libs::DrawLine((int)v(0)+w1,(int)v(1),(int)m->x2+w1,(int)m->y2,
                       col,&out);
      }
    }

  // Show link for inliers (green) and outliers (red)
  it = vec_inliers.begin();
  m = vec_matchings.begin();
  for(size_t i=0; m != vec_matchings.end(); ++m, ++i)
  {
    Image<RGBColor>* im=&out;
    RGBColor col=RED;
    if(it != vec_inliers.end() && i==*it) {
      ++it;
      im=&in;
      col=GREEN;
    }
    libs::DrawLine((int)m->x1,(int)m->y1,(int)m->x2+w1, (int)m->y2, col, im);
  }
  libs::WriteImage(inFile, in);
  libs::WriteImage(outFile, out);
}

int main(int argc, char **argv)
{
  if(argc!=6 && argc!=8 && argc!=9 && argc!=11 && argc!=12) {
    std::cerr << "Usage: " << argv[0] << " imgInA imgInB precision "
              << "allMatches.txt orsaMatches.txt "
              << "[imgInliers imgOutliers [imgMosaic "
              << "[imgMosaicA imgMosaicB [siftRatio]]]"
              << std::endl;
    return 1;
  }

  // Init random seed
  time_t seed = time(0); // Replace by a fixed value to debug a reproducible run
  srand((unsigned int)seed);

  Image<RGBColor> image1, image2;
  if(! libs::ReadImage(argv[1], &image1))
    return 1;
  if(! libs::ReadImage(argv[2], &image2))
    return 1;

  std::istringstream f(argv[3]);
  double precision=0;
  if(! (f>>precision).eof()) {
    std::cerr << "Error reading precision parameter" << std::endl;
    return 1;
  }

  Image<unsigned char> image1Gray, image2Gray;
  libs::convertImage(image1, &image1Gray);
  libs::convertImage(image2, &image2Gray);

  const int w1=image1Gray.Width(), h1=image1Gray.Height();
  const int w2=image2Gray.Width(), h2=image2Gray.Height();

  // SIFT
  float fSiftRatio=0.6f;
  if(argc>11) {
    std::istringstream f(argv[11]);
    if(! (f>>fSiftRatio).eof()) {
      std::cerr << "Error reading SIFT ratio" << std::endl;
      return 1;
    }
  }
  std::vector<Match> vec_matchings;
  SIFT(image1Gray, image2Gray, vec_matchings, fSiftRatio);

  // Remove duplicates (frequent with SIFT)
  rm_duplicates(vec_matchings);

  // Save match files
  if(! Match::saveMatch(argv[4], vec_matchings)) {
    std::cerr << "Failed saving matches into " <<argv[4] <<std::endl;
    return 1;
  }

  // Estimation of homography with ORSA
  libNumerics::matrix<double> H(3,3);
  std::vector<size_t> vec_inliers;
  bool ok = ORSA(vec_matchings, w1, h1, w2, h2, precision, H, vec_inliers);
  if(ok)
  {
    H /= H(2,2);
    std::cout << "H=" << H <<std::endl;
  }

  std::vector<Match> good_match;
  std::vector<size_t>::const_iterator it = vec_inliers.begin();
  for(; it != vec_inliers.end(); it++)
    good_match.push_back(vec_matchings[*it]);
  if(! Match::saveMatch(argv[5], good_match)) {
    std::cerr << "Failed saving matches into " <<argv[5] <<std::endl;
    return 1;
  }

  // Sift de-duplicated output display
  if(argc>7)
  {
    Image<unsigned char> concat;
    ConcatH(image1Gray, image2Gray, concat);
    Image<RGBColor> in;
    libs::convertImage(concat, &in);
    Image<RGBColor> out(in);
    display_match(vec_matchings, vec_inliers, ok? &H: NULL,
                  w1, in,argv[6], out,argv[7]);
  }
  if(! ok)
  {
    std::cerr << "Failed to estimate a model" << std::endl;
    return 1;
  }

  // Mosaics
  if(argc>8)
  {
    std::cout << "-- Render Mosaic -- " << std::endl;

    Rect intersection;
    if(IntersectionBox(w1, h1, w2, h2, H, intersection) &&
      intersection.Width() > 0 && intersection.Height() > 0)
    {
      int xc=(intersection.left+intersection.right)/2;
      int yc=(intersection.top+intersection.bottom)/2;
      int wM = std::max(image1.Width(), image2.Width());
      int hM = std::max(image1.Height(), image2.Height());
      int xo=wM/2;
      int yo=hM/2;
      Image<RGBColor> imageMosaic(wM, hM);
      imageMosaic.fill(WHITE);
      Warp(image1, // Image (Identity)
        image2, // ImageB (will be warped)
        imageMosaic,    // Output image
        H,
        xo-xc, yo-yc);
      libs::WriteImage(argv[8], imageMosaic);

      if (argc>10) // Registered images
      {
        Image<RGBColor> WarpingA(wM, hM);
        Image<RGBColor> WarpingB(wM, hM);
        WarpingA.fill(WHITE);
        WarpingB.fill(WHITE);

        std::cout << "-- Render Mosaic - Image A -- " << std::endl;
        warpImageTo(image1, WarpingA, H, xo-xc, yo-yc);

        std::cout << "-- Render Mosaic - Image B -- " << std::endl;
        H = libNumerics::matrix<double>::eye(3);
        warpImageTo(image2, WarpingB, H, xo-xc, yo-yc);

        libs::WriteImage(argv[9], WarpingA);
        libs::WriteImage(argv[10], WarpingB);
      }
    }
  }
  return 0;
}
