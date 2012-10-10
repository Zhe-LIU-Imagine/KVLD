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

#ifndef LIBS_IMAGE_IMAGE_CONCAT_H_
#define LIBS_IMAGE_IMAGE_CONCAT_H_

#include "image.hpp"

/// Horizontal concatenation of images
template < class Image >
void ConcatH(const Image & imageA, const Image & imageB, Image & Out)
{
  // Compute new dimensions.
  int ww = imageA.Width() + imageB.Width();

  Out = Image(ww, std::max(imageA.Height(), imageB.Height()));

  // Fill with original data from imageA.
  for(size_t i = 0; i < imageA.Width(); ++i)
    for(size_t j = 0; j < imageA.Height(); ++j)
      Out(j,i) = imageA(j,i);

  // Fill with original data from imageB with the imageA Width offset.
  const size_t offset = imageA.Width();
  for(size_t i = 0; i < imageB.Width(); ++i)
    for(size_t j = 0; j < imageB.Height(); ++j)
      Out(j,i+offset) = imageB(j,i);
}

/// Vertical concatenation of images
template < class Image >
void ConcatV(const Image & imageA, const Image & imageB, Image & Out)
{
  // Compute new dimensions.
  int hh = imageA.Height() + imageB.Height();

  Out = Image(std::max(imageA.Width(), imageB.Width()), hh);

  // Fill with original data from imageA.
  for(size_t i = 0; i < imageA.Width(); ++i)
    for(size_t j = 0; j < imageA.Height(); ++j)
      Out(j,i) = imageA(j,i);

  // Fill with original data from imageB with the imageA Height offset.
  const size_t offset = imageA.Height();
  for(size_t i = 0; i < imageB.Width(); ++i)
    for(size_t j = 0; j < imageB.Height(); ++j)
      Out(j+offset,i) = imageB(j,i);
}

#endif // LIBS_IMAGE_IMAGE_CONCAT_H_
