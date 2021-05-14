#include <cryptotemplate/Configuration.h>

#include <cryptotemplate/PolyOps.h>
#include <cryptotemplate/NtlLib.h>

#include <xassert/XAssert.h>

using namespace std;
using namespace NTL;

namespace libcryptotemplate {

ZZ_pX poly_from_roots_ntl(const vec_ZZ_p& roots, long startIncl, long endExcl) {
    // base case is startIncl == endExcl - 1, so return (x - root)
    if(startIncl == endExcl - 1) {
        ZZ_pX monom(NTL::INIT_MONO, 1);
        monom[0] = -roots[startIncl];
        monom[1] = 1;
        assertEqual(NTL::deg(monom), 1);
        return monom;
    }
    
    ZZ_pX p;

    long middle = (startIncl + endExcl) / 2;

    NTL::mul(
        p,
        poly_from_roots_ntl(roots, startIncl, middle),
        poly_from_roots_ntl(roots, middle, endExcl)
        );

    return p;
}

}
