/*
CubicSpline
Copyright (c) 2017 Josh Hooker

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

This program takes a set of arrays/vectors to perform a cubic spline. CubicSpline.h
has options to use std::array so it requires C++11.
*/

#ifndef CUBIC_SPLINE_H
#define CUBIC_SPLINE_H

#include <array>
#include <assert.h>
#include <vector>

#define ASSERT_WITH_MESSAGE(condition, message)                                \
  do {                                                                         \
    if (!(condition)) {                                                        \
      printf((message));                                                       \
    }                                                                          \
    assert((condition));                                                       \
  } while (false)

class CubicSpline {
public:
  CubicSpline();
  template <typename T> CubicSpline(const std::vector<T> &x, const std::vector<T> &y);
  template <typename T, int N, int M> CubicSpline(const T (&x) [N], const T (&y) [M]);
  template <typename T, std::size_t N, std::size_t M> CubicSpline(const std::array<T, N>& x, const std::array<T, M>& y);
  template <typename T> void SetPoints(const std::vector<T> &x, const std::vector<T> &y);
  template <typename T, int N, int M> void SetPoints(const T (&x) [N], const T (&y) [M]);
  template <typename T, std::size_t N, std::size_t M> void SetPoints(const std::array<T, N>& x, const std::array<T, M>& y);
  template <typename T> double operator()(T x) const;
  ~CubicSpline();

private:
  size_t size;
  std::vector<double> xVec, yVec;
  std::vector<double> bVec, cVec, dVec;

  void SetSpline();
};


inline CubicSpline::CubicSpline() {}

template<typename T> inline CubicSpline::CubicSpline(const std::vector<T> &x, const std::vector<T> &y) {
  ASSERT_WITH_MESSAGE(x.size() == y.size(),
    "In CubicSpline initialization, x vector size != y vector size\n");
  assert(x.size() == y.size());
  size = x.size();
  xVec = x; yVec = y;
  bVec.resize(size); cVec.resize(size); dVec.resize(size);

  SetSpline();
}

template <typename T, int N, int M> inline CubicSpline::CubicSpline(const T (&x) [N], const T (&y) [M]) {
  ASSERT_WITH_MESSAGE(N == M,
    "In CubicSpline initialization, x array size != y array size\n");
  assert(N == M);
  size = N;
  xVec.assign(x, x+N); yVec.assign(y, y+M);
  bVec.resize(size); cVec.resize(size); dVec.resize(size);

  SetSpline();
}

template <typename T, std::size_t N, std::size_t M> inline CubicSpline::CubicSpline(const std::array<T, N>& x, const std::array<T, M>& y) {
  ASSERT_WITH_MESSAGE(N == M,
    "In CubicSpline initialization, x array size != y array size\n");
  size = N;
  xVec.resize(size); yVec.resize(size);
  std::copy(x.begin(), x.begin()+size, xVec.begin());
  std::copy(y.begin(), y.begin()+size, yVec.begin());
  bVec.resize(size); cVec.resize(size); dVec.resize(size);

  SetSpline();
}

template <typename T> inline void CubicSpline::SetPoints(const std::vector<T> &x, const std::vector<T> &y) {
  ASSERT_WITH_MESSAGE(x.size() == y.size(),
    "In CubicSpline SetPoints, x vector size != y vector size\n");
  size = x.size();
  xVec = x; yVec = y;
  bVec.resize(size); cVec.resize(size); dVec.resize(size);

  SetSpline();
}

template <typename T, int N, int M> inline void CubicSpline::SetPoints(const T (&x) [N], const T (&y) [M]) {
  ASSERT_WITH_MESSAGE(N == M,
    "In CubicSpline SetPoints, x array size != y array size\n");
  size = N;
  xVec.assign(x, x+N); yVec.assign(y, y+M);
  bVec.resize(size); cVec.resize(size); dVec.resize(size);

  SetSpline();
}

template <typename T, std::size_t N, std::size_t M> inline void CubicSpline::SetPoints(const std::array<T, N>& x, const std::array<T, M>& y) {
  ASSERT_WITH_MESSAGE(N == M,
    "In CubicSpline SetPoints, x array size != y array size\n");
  size = N;
  xVec.resize(size); yVec.resize(size);
  std::copy(x.begin(), x.begin()+size, xVec.begin());
  std::copy(y.begin(), y.begin()+size, yVec.begin());
  bVec.resize(size); cVec.resize(size); dVec.resize(size);

  SetSpline();
}

void inline CubicSpline::SetSpline() {
  std::vector<double> h(size), alpha(size), l(size), z(size), u(size);

  l[0] = 1.;
  u[0] = 0.;
  z[0] = 0.;
  l[size-1] = 1.;
  u[size-1] = 0.;
  cVec[size-1] = 0.;
  for(unsigned int i=0; i<size-1; i++) {
    ASSERT_WITH_MESSAGE(xVec[i+1]>xVec[i],
      "In CubicSpline SetSpline, x array is not sorted from smallest to largest\n");
    assert(xVec[i+1]>xVec[i]);
    h[i] = xVec[i+1]-xVec[i];
    if(i>0) {
      alpha[i] = (3./h[i])*(yVec[i+1]-yVec[i]) - (3./h[i-1])*(yVec[i]-yVec[i-1]);
      l[i] = 2.*(xVec[i+1]-xVec[i-1]) - h[i-1]*u[i-1];
      u[i] = h[i]/l[i];
      z[i] = (alpha[i]-h[i-1]*z[i-1])/l[i];
    }
  }
  for(int i=size-2; i>-1; i--) {
    cVec[i] = z[i]-u[i]*cVec[i+1];
    bVec[i] = (yVec[i+1]-yVec[i])/h[i] - h[i]*(cVec[i+1]+2.*cVec[i])/3.;
    dVec[i] = (cVec[i+1]-cVec[i])/(3.*h[i]);
  }
}

template <typename T> inline double CubicSpline::operator()(T x) const{
  double xs = static_cast<double>(x);

  int l = 0;
  int h = size;
  while(l<h) {
    int mid = (l+h)/2;
    if(xs <= xVec[mid]) {
      h = mid;
    } else {
      l = mid+1;
    }
  }

  size_t idx = l-1;

  double xi = xs-xVec[idx];
  double result;
  if(idx==0) result = yVec[0]+bVec[0]*xi+cVec[0]*xi*xi;
  else if(idx==size-1) result = yVec[size-1]+bVec[size-1]*xi+cVec[size-1]*xi*xi;
  else result = yVec[idx]+bVec[idx]*xi+cVec[idx]*xi*xi+dVec[idx]*xi*xi*xi;
  return result;
}

inline CubicSpline::~CubicSpline() {}

#endif
