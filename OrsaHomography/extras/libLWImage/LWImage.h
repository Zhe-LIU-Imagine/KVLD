// Copyright (c) 2010 Pascal Monasse
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

#ifndef LWIMAGE_H
#define LWIMAGE_H

#include <cstdlib>
#include <string.h>

/// A lightweight image container around an array. It is not responsible for
/// allocation/free of the array.
template <typename T>
class LWImage
{
public:
    LWImage();
    LWImage(T* i_data, int i_w, int i_h, int i_comps=1);
    ~LWImage() {}

    bool valid(int x, int y) const;
    int sizeBuffer() const { return w*h*comps; }
    int step() const;
    int stepComp() const;
    T* pixel(int x, int y);
    const T* pixel(int x, int y) const;
    const T* pixel_ext(int x, int y) const;

    T* data; ///< Array of pixels.
    int w; ///< Width of image.
    int h; ///< Height of image.
    int comps; ///< Components per pixel.
    bool planar; ///< Are components separated (true) or contiguous (false)
};

/// Utility function, avoiding the need to precise the type:
/// \code
///   void f(LWImage<unsigned char>&);
///   unsigned char* data = new unsigned char[10*10];
///   f( make_image(data, 10, 10) );
/// \endcode
template <typename T>
LWImage<T> make_image(T* data, int w, int h, int comps=1)
{
    return LWImage<T>(data, w, h, comps);
}

/// Do not forget free
template <typename T>
LWImage<T> alloc_image(int w, int h, int comps=1)
{
    return LWImage<T>((T*)malloc(w*h*comps*sizeof(T)), w, h, comps);
}

/// Do not forget free
template <typename T>
LWImage<T> alloc_image(const LWImage<T>& im)
{
    LWImage<T> out = alloc_image<T>(im.w, im.h, im.comps);
    out.planar = im.planar;
    //for(int i=im.sizeBuffer()-1; i>=0; i--)
    //    out.data[i] = im.data[i];
    memcpy(out.data, im.data, im.sizeBuffer()*sizeof(float));
    return out;
}

#include "LWImage.cpp"

#endif
