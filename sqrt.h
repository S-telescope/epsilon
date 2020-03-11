#pragma once
// sqrt.h - sqrtonential function
#pragma once
#include "epsilon.h"
#include "multi_epsilon.h"
#include <vector>
#include <cmath>
namespace fms {
	//(1+x)^{1/2}=\sum_{k=0}^\infty\binom{1/2}{k}x^k for abs(x)<1
	template<class X>
	inline X sqrt_brute(const X& x)
	{
		
		auto ex = 0 * x + 1;
		X xn_(x-1); // (x-1)^n * nCr(1/2,n) 
		xn_ /= 2;
		int n = 0;
		while (fabs(xn_) + X(1) != X(1)) {
			
			ex += xn_;
			xn_ *= (x-1) * (0.5 - n) / ++n ;
		}

		return ex;
	}
	template<>
	inline double sqrt_brute(const double& x) {
		return ::sqrt(x);
	}

	template<size_t N, class X = double>
	inline epsilon<N, X> sqrt(const epsilon<N, X>& x) {
		assert(x[0] > 0);
		epsilon<N, X> res = ::sqrt(x[0]) * sqrt_brute(1 + (x - x[0])/x[0] );
		return res;
	}

	/*inline double sqrt(const double& x) {
		return ::sqrt(x);
	}*/

	inline multi_epsilon sqrt(const multi_epsilon& x) {
		assert(x[0] > 0);
		auto res = ::sqrt(x[0]) * sqrt_brute(1 + (x - x[0]) / x[0]);
		return res;
	}
}
