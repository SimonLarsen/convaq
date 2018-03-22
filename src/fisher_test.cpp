// The MIT License (MIT)
//
// Copyright (c) 2014 Julian Gehring
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include <boost/math/distributions/hypergeometric.hpp>
#include <algorithm>

double fisher_test(int a, int b, int c, int d) {
  unsigned N = a + b + c + d;
  unsigned r = a + c;
  unsigned n = c + d;
  unsigned max_for_k = std::min(r, n);
  unsigned min_for_k = (unsigned)std::max(0, int(r + n - N));
  boost::math::hypergeometric_distribution<> hgd(r, n, N);
  double cutoff = pdf(hgd, c);
  double tmp_p = 0.0;
  for(int k = min_for_k;k < max_for_k + 1;k++) {
    double p = pdf(hgd, k);
    if(p <= cutoff) tmp_p += p;
  }
  return tmp_p;
}
