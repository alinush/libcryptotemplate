#pragma once

#include <vector>
#include <algorithm>
#include <cmath>
#include <vector>
#include <iostream>

#include <cryptotemplate/NtlLib.h>
#include <libfqfft/polynomial_arithmetic/basic_operations.hpp>

#include <xassert/XAssert.h>
#include <xutils/Log.h>
#include <xutils/Utils.h>

using namespace std;
using namespace libfqfft;

namespace libcryptotemplate {

/**
 * Does an FFT: i.e., evaluates and returns p(\omega^i) for all i = 0, ..., N-1
 * 
 * @param   N   specifies the primitive Nth root of unity \omega (must be 2^k)
 */
template<class FieldT>
void poly_fft(const std::vector<FieldT>& p, size_t N, std::vector<FieldT>& vals) {
    assertTrue(Utils::isPowerOfTwo(N));
    FieldT omega = libff::get_root_of_unity<FieldT>(N);

    if(N < p.size())
        throw std::runtime_error("N has to be greater than polynomial degree");

    vals = p;
    vals.resize(N, FieldT::zero());

#ifdef MULTICORE
    _basic_parallel_radix2_FFT(vals, omega);
#else
    _basic_serial_radix2_FFT(vals, omega);
#endif
}

// WARNING: Slower than NTL::BuildFromRoots(), but not sure why
ZZ_pX poly_from_roots_ntl(const vec_ZZ_p& roots, long startIncl, long endExcl);

template<typename FieldT>
std::vector<FieldT> poly_from_roots(const vector<FieldT>& r) {
    std::vector<FieldT> a;
    vec_ZZ_p roots;
    roots.SetLength(static_cast<long>(r.size()));

    mpz_t rop;
    mpz_init(rop);
    ZZ mid;

    for(size_t i = 0; i < r.size(); i++) {
        libcryptotemplate::conv_single_fr_zp(r[i], roots[static_cast<long>(i)], rop, mid);
    }
    mpz_clear(rop);

    ZZ_pX acc(NTL::INIT_MONO, 0);
    NTL::BuildFromRoots(acc, roots);
    //acc = poly_from_roots_ntl(roots, 0, roots.length());  // WARNING: Slower than NTL::BuildFromRoots()

    libcryptotemplate::convNtlToLibff(acc, a);

    return a;
}

// interp(roots[0,n)) = interp(roots[0, n/2)) * interp(roots[n/2+1, n))
template<typename FieldT>
std::vector<FieldT> poly_from_roots_slow(const vector<FieldT>& roots, size_t startIncl, size_t endExcl) {
    // base case is startIncl == endExcl - 1, so return (x - root)
    if(startIncl == endExcl - 1) {
        std::vector<FieldT> monom(2);
        monom[0] = -roots[startIncl];
        monom[1] = 1;
        return monom;
    }

    std::vector<FieldT> p;

    size_t middle = (startIncl + endExcl) / 2;

    libfqfft::_polynomial_multiplication_on_fft(
        p,
        poly_from_roots_slow(roots, startIncl, middle),
        poly_from_roots_slow(roots, middle, endExcl)
        );

    return p;
}

template<typename FieldT>
std::vector<FieldT> poly_from_roots_slow(const vector<FieldT>& roots) {
    return poly_from_roots(roots, 0, roots.size());
}

/**
 Calculates derivative of polynomial f and stores in d
 */
template<typename FieldT>
void poly_differentiate(const vector<FieldT> &f, vector<FieldT> &d) {
    assertStrictlyGreaterThan(f.size(), 0);

    d.resize(f.size() - 1);
    for (size_t i = 0; i < d.size(); i++) {
        d[i] = f[i + 1] * (static_cast<int>(i) + 1);
    }
}

/**
 * O(n^2) naive polynomial multiplication
 */
template<typename FieldT>
void polynomial_multiplication_naive(vector<FieldT> &c, const vector<FieldT> &b, const vector<FieldT> &a) {
    c.resize(a.size() + b.size() - 1, 0);
    for (size_t i = 0; i < a.size(); i++)
        for (size_t j = 0; j < b.size(); j++)
            c[i+j] += a[i] * b[j];
}

/**
 * Returns true if the polynomial is of the type X^m - c, where m != 0
 */
template<class FieldT>
bool poly_is_xnc(const vector<FieldT>& p) {
    if(p.size() < 2)
        return false;

    if(p[p.size() - 1] != FieldT::one())
        return false;

    for(size_t j = 1; j < p.size() - 1; j++) {
        if(p[j] != FieldT::zero())
            return false;
    }

    return true;
}

template<class FieldT>
void poly_print(const std::vector<FieldT>& p) {
    size_t n = p.size();

    for(size_t i = 0; i < n - 1; i++) {
        auto deg = (n-1) - i;
        auto c = p[deg];

        if(deg == 1)
            std::cout << c << " x + ";
        else
            std::cout << c << " x^" << deg << " + ";
    }

    std::cout << p[0];

    std::cout << " of degree " << n-1 << std::flush;
}

/**
 * Prints a polynomial in a more-user friendly way if it has coefficients that
 * are w_n^k roots of unity, where n = allRoots.size().
 */
template<class FieldT>
void poly_print_wnk(const std::vector<FieldT>& p, const std::vector<FieldT>& allRoots) {
    size_t d = p.size();
    FieldT minusOne = FieldT(-1);
    for(size_t i = 0; i < d; i++) {
        auto deg = (d-1) - i;
        auto c = p[deg];

        // don't print (0 x^i) coeffs
        if(c == 0)
            continue;

        // identify if the coefficient is a root of unity (special treatment for 1 and -1)
        auto it = std::find(allRoots.begin(), allRoots.end(), c);

        // print the coefficient
        if(it != allRoots.end()) {
            std::string coeff;
            bool isPositive = (c != minusOne);

            if (c == FieldT::one() || c == minusOne) {
                if(deg == 0) {
                    coeff = "1";
                } else {
                    coeff = "";
                }
            } else {
                size_t pos = static_cast<size_t>(it - allRoots.begin());
                size_t n = allRoots.size();
                coeff = "w_" + std::to_string(n) + "^" + std::to_string(pos);
            }

            if(deg != d-1 || !isPositive) {
                std::cout << (isPositive ? " + " : " - ");
            }

            std::cout << coeff;

        } else {
            throw std::runtime_error("Non-root-of-unity coefficient found");
        }

        // print x^deg after the coefficient
        if(deg == 1) {
            std::cout << " x";
        } else if(deg == 0) {
            // for the last coefficient we don't print any x^deg (because it would be x^0)
        } else {
            std::cout << " x^" << deg; 
        }
    }

    std::cout << " of degree " << d-1 << std::flush;
}

template<class T>
bool poly_equal_slow(const std::vector<T>& a, const std::vector<T>& b) {
    std::vector<T> res;
    _polynomial_subtraction(res, a, b);
    
    return res.size() == 0;
}

}
