#include "sqrt.h"
#include <cassert>
using namespace fms;
template<class X=double>
int test_sqrt() {
	X x = 1.0;
	//assert(fms::sqrt_brute(x) == ::sqrt(x));

	auto x_ = multi_epsilon::identity(2, 2) * x;
	//assert(fms::sqrt_brute(x_) == multi_epsilon::identity(2, 2) * ::sqrt(x));

	x = 4.0;
	x_ = multi_epsilon::identity(2, 2) * x;
	//assert(fms::sqrt(x_) == multi_epsilon::identity(2, 2) * ::sqrt(x));

	X a = 4.0;
	X b = 40.0;
	x_ = multi_epsilon({ a,b }, 1, 1);
	assert(fms::sqrt(x_)[0] == ::sqrt(a));
	assert(fms::sqrt(x_)[1] == b/2/::sqrt(a));
	
	return 0;
}
int test_sqrt_double = test_sqrt();