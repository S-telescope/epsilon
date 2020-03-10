#pragma once
// exp.h - exponential function
#pragma once
#include "epsilon.h"
#include "multi_epsilon.h"
#include <vector>
#include <cmath>
namespace fms {

	template<class X>
	inline X exp_brute(const X& x)
	{
		X ex = 1;
		X xn_(x); // x^n/n! 1

		int n = 1;
		while (fabs(xn_) + X(1) != X(1)) {
			ex += xn_;
			xn_ *= x / ++n;
		}

		return ex;
	}
	template<>
	inline double exp_brute(const double& x) {
		return ::exp(x);
	}
	
	template<size_t N, class X = double>
	epsilon<N, X> exp(const epsilon<N, X>& x) {
		epsilon<N, X> res(exp_brute(x[0]));
		res *= exp_brute(x - epsilon<N, X>(x[0]));
		return res;
	}

	inline double exp(const double& x) {
		return ::exp(x);
	}

	multi_epsilon exp(const multi_epsilon& x) {
		auto res=::exp(x[0])+0*x;
		res *= exp_brute(x - res);
		return res;
	}
}
