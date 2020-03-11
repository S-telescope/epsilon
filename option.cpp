#include "normal.h"
template <class X>
X bs_call(X S, X T, X r, X q, X v, X K) {
	auto d1= (log(S / K) + (r +  0.5 * v * v) * T) / (v * sqrt(T));
	auto d2 = (log(S / K) + (r - 0.5 * v * v) * T) / (v * sqrt(T));
	return S * exp(-q * T) * cdf(d1) + K * exp(-r * T) * cdf(d2);
}