#include <cassert>
#include "normal.h"
#include <iostream>

using namespace fms;

template<class X>
int test_erf()
{
	X x = 1;
	X ex = fms::erf_brute<X>(x);

	X dx = ex - ::erf(x);
	assert(fabs(dx) <= 1 * std::numeric_limits<X>::epsilon());

	return 0;
}
static int test_erf_double = test_erf<double>();


template<size_t N, class X>
int test_erfn()
{
	X x = 1;
	epsilon<N, X> x_(x, 1);
	auto ex = fms::erf_brute(x_);

	X erfx = ::erf(x);
	assert(fabs(ex[0] - erfx) <= 1 * std::numeric_limits<X>::epsilon());
	if constexpr (N > 1) {
		assert(fabs(ex[1] - M_2_SQRTPI * ::exp(-x * x)) <= 1 * std::numeric_limits<X>::epsilon());
		assert(fabs(ex[2] + 2 * x * M_2_SQRTPI * ::exp(-x * x)) <= 1 * std::numeric_limits<X>::epsilon());
	}

	X a = 1;
	X b = 1;
	X c = 1;
	x_ = epsilon<N, X>({ a,b,c });
	ex = fms::erf_brute(x_);

	assert(fabs(ex[0] - ::erf(a)) <= 1 * std::numeric_limits<X>::epsilon());
	if constexpr (N > 1) {
		assert(fabs(ex[1] - M_2_SQRTPI * ::exp(-a * a) * b) <= 1 * std::numeric_limits<X>::epsilon());
		assert(fabs(ex[2] - M_2_SQRTPI * ::exp(-a * a) * (c - 2 * a * b * b)) <= 1 * std::numeric_limits<X>::epsilon());
	}

	auto x__ = multi_epsilon({ a,b,c/2 },1,2);
	auto ex_ = fms::erf_brute(x__);

	assert(fabs(ex_[0] - ::erf(a)) <= 1 * std::numeric_limits<X>::epsilon());
	if constexpr (N > 1) {
		assert(fabs(ex_[1] - M_2_SQRTPI * ::exp(-a * a) * b) <= 1 * std::numeric_limits<X>::epsilon());
		assert(fabs(ex_[2] - M_2_SQRTPI * ::exp(-a * a) * (c - 2 * a * b * b)/2) <= 1 * std::numeric_limits<X>::epsilon());
	}
	return 0;
}
static int test_erfn_double = test_erfn<3, double>();

template<size_t N, class X>
int test_erf2n()
{
	X x = 1;
	epsilon<N, X> x_(x, 1);
	auto ex = fms::erf(x_);

	X erfx = ::erf(x);
	assert(fabs(ex[0] - erfx) <= 0 * std::numeric_limits<X>::epsilon());
	if constexpr (N > 1) {
		assert(fabs(ex[1] - M_2_SQRTPI * ::exp(-x * x)) <= 0 * std::numeric_limits<X>::epsilon());
		assert(fabs(ex[2] + 2 * x * M_2_SQRTPI * ::exp(-x * x)) <= 0 * std::numeric_limits<X>::epsilon());
	}

	X a = 1;
	X b = 1;
	X c = 1;
	x_ = epsilon<N, X>({ a,b,c });
	ex = fms::erf(x_);


	assert(fabs(ex[0] - ::erf(a)) <= 0 * std::numeric_limits<X>::epsilon());
	if constexpr (N > 1) {
		assert(fabs(ex[1] - M_2_SQRTPI * ::exp(-a * a) * b) <= 0 * std::numeric_limits<X>::epsilon());
		assert(fabs(ex[2] - M_2_SQRTPI * ::exp(-a * a) * (c - 2 * a * b * b)) <= 0 * std::numeric_limits<X>::epsilon());
	}

	auto x__ = multi_epsilon({ a,b,c/2 }, 1,2);
	auto ex_ = fms::erf(x__);


	assert(fabs(ex_[0] - ::erf(a)) <= 0 * std::numeric_limits<X>::epsilon());
	if constexpr (N > 1) {
		assert(fabs(ex_[1] - M_2_SQRTPI * ::exp(-a * a) * b) <= 0 * std::numeric_limits<X>::epsilon());
		assert(fabs(ex_[2] - M_2_SQRTPI * ::exp(-a * a) * (c - 2 * a * b * b) / 2) <= 0 * std::numeric_limits<X>::epsilon());
	}
	return 0;
}
static int test_erf2n_double = test_erf2n<3, double>();

template<class X>
int test_cdf()
{
	X x = 1;
	X ex = fms::cdf_brute<X>(x);

	X dx = ex - 0.5 * erfc(-x * M_SQRT1_2);
	assert(fabs(dx) <= 0 * std::numeric_limits<X>::epsilon());


	return 0;
}
static int test_cdf_double = test_cdf<double>();

template<size_t N, class X>
int test_cdfn()
{
	X x = 1;
	epsilon<N, X> x_(x, 1);
	auto ex = fms::cdf_brute(x_);

	X cdfx = 0.5 * erfc(-x * M_SQRT1_2);
	assert(fabs(ex[0] - cdfx) <= 1 * std::numeric_limits<X>::epsilon());
	if constexpr (N > 1) {
		assert(fabs(ex[1] - (1 / (sqrt(2 * M_PI))) * exp(-0.5 * x * x)) <= 1 * std::numeric_limits<X>::epsilon());
		assert(fabs(ex[2] + (1 / (sqrt(2 * M_PI))) * x * exp(-0.5 * x * x)) <= 1 * std::numeric_limits<X>::epsilon());
	}

	X a = 1;
	X b = 2;
	X c = 1;
	x_ = epsilon<N, X>({ a,b,c });
	ex = fms::cdf_brute(x_);

	assert(fabs(ex[0] - 0.5 * erfc(-a * M_SQRT1_2)) <= 1 * std::numeric_limits<X>::epsilon());
	if constexpr (N > 1) {
		assert(fabs(ex[1] - 0.5 * M_SQRT1_2 * M_2_SQRTPI * b * ::exp(-0.5 * a * a)) <= 1 * std::numeric_limits<X>::epsilon());
		assert(fabs(ex[2] - 0.5 * M_2_SQRTPI * M_SQRT1_2 * c * ::exp(-0.5 * a * a) + 0.5 * M_2_SQRTPI * M_SQRT1_2 * a * b * b * ::exp(-0.5 * a * a)) <= 2 * std::numeric_limits<X>::epsilon());
	}

	auto x__ = multi_epsilon({ a,b,c/2 }, 1,2);
	auto ex_ = fms::cdf_brute(x__);

	assert(fabs(ex_[0] - 0.5 * erfc(-a * M_SQRT1_2)) <= 1 * std::numeric_limits<X>::epsilon());
	if constexpr (N > 1) {
		assert(fabs(ex_[1] - 0.5 * M_SQRT1_2 * M_2_SQRTPI * b * ::exp(-0.5 * a * a)) <= 1 * std::numeric_limits<X>::epsilon());
		assert(fabs(ex_[2]*2 - 0.5 * M_2_SQRTPI * M_SQRT1_2 * c * ::exp(-0.5 * a * a) + 0.5 * M_2_SQRTPI * M_SQRT1_2 * a * b * b * ::exp(-0.5 * a * a)) <= 2 * std::numeric_limits<X>::epsilon());
	}
	return 0;
}

static int test_cdfn_double = test_cdfn<3, double>();

template<size_t N, class X>
int test_cdf2n()
{
	X x = 1;
	epsilon<N, X> x_(x, 1);
	auto ex = fms::cdf(x_);

	X cdfx = 0.5 * erfc(-x * M_SQRT1_2);
	//X cdfx = 0.5 + 0.5 * ::erf(x * M_SQRT1_2);
	assert(fabs(ex[0] - cdfx) <= 1 * std::numeric_limits<X>::epsilon());
	if constexpr (N > 1) {
		assert(fabs(ex[1] - (1 / (sqrt(2 * M_PI))) * exp(-0.5 * x * x)) <= 1 * std::numeric_limits<X>::epsilon());
		assert(fabs(ex[2] + (1 / (sqrt(2 * M_PI))) * x * exp(-0.5 * x * x)) <= 1 * std::numeric_limits<X>::epsilon());
	}

	X a = 1;
	X b = 2;
	X c = 5;
	x_ = epsilon<N, X>({ a,b,c });
	ex = fms::cdf(x_);

	assert(fabs(ex[0] - 0.5 * erfc(-a * M_SQRT1_2)) <= 0 * std::numeric_limits<X>::epsilon());
	if constexpr (N > 1) {
		assert(fabs(ex[1] - 0.5 * M_SQRT1_2 * M_2_SQRTPI * b * ::exp(-0.5 * a * a)) <= 1 * std::numeric_limits<X>::epsilon());
		assert(fabs(ex[2] - 0.5 * M_2_SQRTPI * M_SQRT1_2 * c * ::exp(-0.5 * a * a) + 0.5 * M_2_SQRTPI * M_SQRT1_2 * a * b * b * ::exp(-0.5 * a * a)) <= 1 * std::numeric_limits<X>::epsilon());
	}

	auto x__ = multi_epsilon({ a,b,c/2 },1,2);
	auto ex_ = fms::cdf(x__);

	assert(fabs(ex_[0] - 0.5 * erfc(-a * M_SQRT1_2)) <= 0 * std::numeric_limits<X>::epsilon());
	if constexpr (N > 1) {
		assert(fabs(ex_[1] - 0.5 * M_SQRT1_2 * M_2_SQRTPI * b * ::exp(-0.5 * a * a)) <= 1 * std::numeric_limits<X>::epsilon());
		assert(fabs(ex_[2]*2 - 0.5 * M_2_SQRTPI * M_SQRT1_2 * c * ::exp(-0.5 * a * a) + 0.5 * M_2_SQRTPI * M_SQRT1_2 * a * b * b * ::exp(-0.5 * a * a)) <= 1 * std::numeric_limits<X>::epsilon());
	}
	return 0;
}

static int test_cdf2n_double = test_cdf2n<3, double>();

int test_multi_erf() {
	using X=double;
	X a = 1;
	X b = 2;
	X c = 3;
	X d = 4;
	auto x__ = multi_epsilon({ a,b,c,d }, 2, 1);
	auto ex_ = fms::erf(x__);
	assert(fabs(ex_[0] - ::erf(a)) <= 0 * std::numeric_limits<X>::epsilon());
	
	assert(fabs(ex_[1] - M_2_SQRTPI * ::exp(-a * a) * b) <= 0 * std::numeric_limits<X>::epsilon());
	assert(fabs(ex_[2] - M_2_SQRTPI * ::exp(-a * a) * c) <= 1 * std::numeric_limits<X>::epsilon());
	assert(fabs(ex_[3] - M_2_SQRTPI * ::exp(-a * a) * (d - 2* a * b * c)) <= 2 * std::numeric_limits<X>::epsilon());
	
	return 0;
}

static int test_multi_erf_double = test_multi_erf();

int test_multi_cdf() {
	using X=double;
	X a = 1;
	X b = 2;
	X c = 3;
	X d = 4;
	auto x__ = multi_epsilon({ a,b,c,d }, 2, 1);
	auto ex_ = fms::cdf(x__);
	assert(fabs(ex_[0] - 0.5 * erfc(-a * M_SQRT1_2)) <= 0 * std::numeric_limits<X>::epsilon());

	assert(fabs(ex_[1] - 0.5 * M_SQRT1_2 * M_2_SQRTPI * ::exp(-a * a/2) * b) <= 2 * std::numeric_limits<X>::epsilon());
	assert(fabs(ex_[2] - 0.5 * M_SQRT1_2 * M_2_SQRTPI * ::exp(-a * a/2) * c) <= 1 * std::numeric_limits<X>::epsilon());
	assert(fabs(ex_[3] - 0.5 * M_SQRT1_2 * M_2_SQRTPI * ::exp(-a * a/2) * (d - a * b * c)) <= 1 * std::numeric_limits<X>::epsilon());

	return 0;
}

static int test_multi_cdf_double = test_multi_cdf();