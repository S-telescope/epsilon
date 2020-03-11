#pragma once
#pragma once
// log.h - log function
#pragma once
#include "epsilon.h"
#include "multi_epsilon.h"
#include <vector>
#include <cmath>
namespace fms {
	//(1+x)^{1/2}=\sum_{k=0}^\infty\binom{1/2}{k}x^k for abs(x)<1
	template<class X>
	inline X log_brute(const X& x)
	{

		auto ex = 0 * x;
		X xn_(x - 1); // (1-x)^n / n
		int n = 1;
		while (fabs(xn_) + 1 != 1) {
			
			ex += xn_ / n;
			xn_ *= (1-x);
			++n;
		}

		return ex;
	}

	template<>
	inline double log_brute(const double& x) {
		return ::log(x);
	}

	template<size_t N, class X = double>
	inline epsilon<N, X> log(const epsilon<N, X>& x) {
		assert(x[0] > 0);
		epsilon<N, X> res = ::log(x[0]) + log_brute(1 + (x - x[0]) / x[0]);
		return res;
	}

	/*inline double log(const double& x) {
		return ::log(x);
	}*/

	inline multi_epsilon log(const multi_epsilon& x) {
		assert(x[0] > 0);
		auto res = ::log(x[0]) + log_brute(1 + (x - x[0]) / x[0]);
		return res;
	}
}

/*inline fms::multi_epsilon log(const fms::multi_epsilon& a)
{
	return fms::log(a);
}*/