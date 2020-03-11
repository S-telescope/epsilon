#include "sqrt.h"
#include <cassert>
using namespace fms;
template<class X=double>
int test_sqrt() {
	X x = 1.5;
	assert(fms::sqrt_brute(x) == ::sqrt(x));

	auto x_ = multi_epsilon::identity(2, 2) * x;
	assert(fabs(fms::sqrt_brute(x_) - multi_epsilon::identity(2, 2) * ::sqrt(x)) <= 2 * std::numeric_limits<X>::epsilon());
	
	x = 4.0;
	x_ = multi_epsilon::identity(2, 2) * x;
	assert(fms::sqrt(x_) == multi_epsilon::identity(2, 2) * ::sqrt(x));

	X a = 4.0;
	X b = 40.0;
	x_ = multi_epsilon({ a,b }, 1, 1);
	assert(fms::sqrt(x_)[0] == ::sqrt(a));
	assert(fms::sqrt(x_)[1] == b/2/::sqrt(a));
	
	
	X c = 5.0;
	X d = 1.0;
	x_ = multi_epsilon({ a,b,c,d }, 2, 1);
	auto sqrtx_ = fms::sqrt(x_);
	assert(sqrtx_[0] == ::sqrt(a));
	assert(sqrtx_[1] == b / 2 / ::sqrt(a));
	assert(sqrtx_[2] == c / 2 / ::sqrt(a));
	assert(sqrtx_[3] == (2*a*d-b*c)/4/a / ::sqrt(a));

	auto y_ = epsilon<4, X>({ a,b,c,d });
	assert(fms::sqrt(y_)[0] == ::sqrt(a));
	assert(fms::sqrt(y_)[1] == b/2/::sqrt(a));
	assert(fms::sqrt(y_)[2] == (2*a*c-b*b)/4/a / ::sqrt(a));
	assert(fms::sqrt(y_)[3] == (4*a*a*d-6*a*b*c+3*b*b*b)/8/a/a/::sqrt(a));

	x_ = multi_epsilon({ a,b,c/2,d/6 },1,3);
	assert(fms::sqrt(x_)[0] == ::sqrt(a));
	assert(fms::sqrt(x_)[1] == b / 2 / ::sqrt(a));
	assert(fms::sqrt(x_)[2]*2 == (2 * a * c - b * b) / 4 / a / ::sqrt(a));
	assert(fms::sqrt(x_)[3]*6 == (4 * a * a * d - 6 * a * b * c + 3 * b * b * b) / 8 / a / a / ::sqrt(a));

	return 0;
}
static int test_sqrt_double = test_sqrt();