#include <cassert>
#include "exp.h"
#include <iostream>
using namespace fms;

template<class X>
int test_exp()
{
	X x = 1;
	X ex = fms::exp_brute<X>(x);

	X dx = ex - ::exp(x);
	assert(fabs(dx) <= 2 * std::numeric_limits<X>::epsilon());

	return 0;
}
static int test_exp_double = test_exp<double>();

template<size_t N, class X>
int test_expn()
{
	X x = 1;
	epsilon<N, X> x_(x, 1);
	auto ex = fms::exp_brute(x_);

	X expx = ::exp(x);
	assert(fabs(ex[0] - expx) <= 2 * std::numeric_limits<X>::epsilon());
	if constexpr (N > 1) {
		assert(fabs(ex[1] - expx) <= 2 * std::numeric_limits<X>::epsilon());
	}

	return 0;
}
static int test_expn_double = test_expn<2, double>();

template<size_t N, class X>
int test_expn_exp2n_at_small()
{
	X x = 0.1;
	epsilon<N, X> x_(x, 1);
	auto ex = fms::exp_brute(x_);

	X expx = ::exp(x);
	assert(fabs(ex[0] - expx) <= 2 * std::numeric_limits<X>::epsilon());
	if constexpr (N > 1) {
		assert(fabs(ex[1] - expx) <= 2 * std::numeric_limits<X>::epsilon());
	}


	return 0;
}
static int test_expn_exp2n_at_small_double = test_expn_exp2n_at_small<2, double>();

template<size_t N, class X>
int test_expn_exp2n_at_large()
{
	X x = 10;
	epsilon<N, X> x_(x, 1);
	auto ex = fms::exp_brute(x_);

	X expx = ::exp(x);

	assert(fabs(ex[0] - expx) <= 16384 * std::numeric_limits<X>::epsilon());
	if constexpr (N > 1) {
		assert(fabs(ex[1] - expx) <= 16384 * std::numeric_limits<X>::epsilon());
	}

	

	return 0;
}
static int test_expn_exp2n_at_large_double = test_expn_exp2n_at_large<2, double>();

template<size_t N, class X>
int test_expn_higher_order()
{
	epsilon<N, X> x_({ 10,2,0.1 });
	
	std::vector<X> expx({ ::exp(10),::exp(10) * 2, ::exp(10) * 4.1 });

	

	auto ex = fms::exp(x_);
	assert(fabs(ex[0] - expx[0]) <= 0 * std::numeric_limits<X>::epsilon());
	if constexpr (N > 1) {
		assert(fabs(ex[1] - expx[1]) <= 0 * std::numeric_limits<X>::epsilon());
		assert(fabs(ex[2] - expx[2]) <= 0 * std::numeric_limits<X>::epsilon());
	}
	return 0;
}
static int test_expn_higher_order_double = test_expn_higher_order<3, double>();

int test_multi_exp() {
	using X= double;
	X x = 10;
	auto x_ = fms::multi_epsilon::identity(2,2)*x;
	auto ex = fms::exp(x_);

	auto dx = ex - ::exp(x);
	assert(fabs(dx) <= 2 * std::numeric_limits<X>::epsilon());

	auto x__=fms::multi_epsilon({ 10,2,0.1 }, 1, 2);

	std::vector<X> expx({ ::exp(10),::exp(10) * 2, ::exp(10) * 2.1 });

	auto ex_ = fms::exp(x__);
	assert(fabs(ex_[0] - expx[0]) <= 0 * std::numeric_limits<X>::epsilon());
	assert(fabs(ex_[1] - expx[1]) <= 0 * std::numeric_limits<X>::epsilon());
	assert(fabs(ex_[2] - expx[2]) <= 0 * std::numeric_limits<X>::epsilon());
	
	auto y_ = fms::multi_epsilon({ 10,2,0.1,10 }, 2, 1);

	std::vector<X> expy({ ::exp(10),::exp(10) * 2, ::exp(10)*0.1, ::exp(10) * 10.2 });

	auto ey_ = fms::exp(y_);
	assert(fabs(ey_[0] - expy[0]) <= 0 * std::numeric_limits<X>::epsilon());
	assert(fabs(ey_[1] - expy[1]) <= 0 * std::numeric_limits<X>::epsilon());
	assert(fabs(ey_[2] - expy[2]) <= 0 * std::numeric_limits<X>::epsilon());
	assert(fabs(ey_[3] - expy[3]) <= 0 * std::numeric_limits<X>::epsilon());

	
	
	return 0;
}

static int test_multi_exp_double = test_multi_exp();

