#include "log.h"
#include <cassert>
using namespace fms;
template<class X = double>
int test_log() {
	X x = 1.5;
	assert(fms::log_brute(x) == ::log(x));

	auto x_ = multi_epsilon::identity(2, 2) * x;
	assert(fabs(fms::log_brute(x_) - multi_epsilon::identity(2, 2) * ::log(x)) <= 2 * std::numeric_limits<X>::epsilon());

	x = 4.0;
	x_ = multi_epsilon::identity(2, 2) * x;
	assert(fms::log(x_) == multi_epsilon::identity(2, 2) * ::log(x));

	X a = 4.0;
	X b = 40.0;
	x_ = multi_epsilon({ a,b }, 1, 1);
	assert(fms::log(x_)[0] == ::log(a));
	assert(fms::log(x_)[1] == b/a);


	X c = 5.0;
	X d = 1.0;
	x_ = multi_epsilon({ a,b,c,d }, 2, 1);
	auto logx_ = fms::log(x_);
	assert(logx_[0] == ::log(a));
	assert(logx_[1] == b / a);
	assert(logx_[2] == c / a);
	assert(logx_[3] == (a * d - b * c) / a / a);

	auto y_ = epsilon<4, X>({ a,b,c,d });
	assert(fms::log(y_)[0] == ::log(a));
	assert(fms::log(y_)[1] == b / a);
	assert(fms::log(y_)[2] == (a * c - b * b) / a / a);
	assert(fms::log(y_)[3] == (a * a * d - 3 * a * b * c + 2 * b * b * b) / a / a / a);

	x_ = multi_epsilon({ a,b,c / 2,d / 6 }, 1, 3);
	assert(fms::log(x_)[0] == ::log(a));
	assert(fms::log(x_)[1] == b / a);
	assert(fms::log(x_)[2] * 2 == (a * c - b * b) / a / a);
	assert(fms::log(x_)[3] * 6 == (a * a * d - 3 * a * b * c + 2 * b * b * b) / a / a / a);

	return 0;
}
static int test_log_double = test_log();