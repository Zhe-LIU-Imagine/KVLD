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

#ifdef LWIMAGE_H // Do nothing if not included from LWImage.h

/// Constructor.
template <typename T>
LWImage<T>::LWImage()
  : data(0), w(0), h(0), comps(0), planar(true)
{}

/// Constructor.
template <typename T>
LWImage<T>::LWImage(T* i_data, int i_w, int i_h, int i_comps)
  : data(i_data), w(i_w), h(i_h), comps(i_comps), planar(true)
{}

/// Is pixel inside image?
template <typename T>
bool LWImage<T>::valid(int x, int y) const
{
    return (0 <= x && x < w && 0 <= y && y < h);
}

/// Step between one pixel and next one.
template <typename T>
int LWImage<T>::step() const
{ return (planar? 1: comps); }

/// Step between one component of pixel and next one.
template <typename T>
int LWImage<T>::stepComp() const
{ return (planar? w*h: 1); }

/// Return pointer to data at pixel.
template <typename T>
T* LWImage<T>::pixel(int x, int y)
{
    return (data + step()*(y*w+x));
}

/// Return pointer to data at pixel.
template <typename T>
const T* LWImage<T>::pixel(int x, int y) const
{
    return (data + step()*(y*w+x));
}

/// Value of x in [0,w-1], with mirror and periodization.
inline void wrap(int& x, int w)
{
    if(x < 0)       x  =-(x+1);
    while(x >= 2*w) x -= 2*w;
    if(x >= w)      x  = 2*w-x-1;
}

/// Return pixel value even outside image: image is supposed infinite, with
/// mirror effect and periodization.
template <typename T>
const T* LWImage<T>::pixel_ext(int x, int y) const
{
    wrap(x,w);
    wrap(y,h);
    return pixel(x,y);
}

#endif // LWIMAGE_H
