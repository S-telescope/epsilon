#include "normal.h"
#include "log.h"
#include "sqrt.h"
#include<iostream>
template <class X>
X bs_call(X S, X T, X r, X q, X v, X K) {
	X d1= (log(S / K) + (r - q +  0.5 * v * v) * T) / (v * sqrt(T));
	X d2 = (log(S / K) + (r - q - 0.5 * v * v) * T) / (v * sqrt(T));
	X e1 = exp(-q * T);
	X e2 = cdf(d1);
	return S * exp(-q * T) * cdf(d1) - K * exp(-r * T) * cdf(d2);
}

template <class X>
X call_delta(X S, X T, X r, X q, X v, X K) {
	X d1 = (log(S / K) + (r - q + 0.5 * v * v) * T) / (v * sqrt(T));
	return exp(-q * T) * cdf(d1);
}

template <class X>
X call_theta(X S, X T, X r, X q, X v, X K) {
	X d1 = (log(S / K) + (r - q + 0.5 * v * v) * T) / (v * sqrt(T));
	X d2 = (log(S / K) + (r - q - 0.5 * v * v) * T) / (v * sqrt(T));
	X phi_d1 = M_2_SQRTPI / 2 * M_SQRT1_2 * exp(-d1 * d1 / 2);
	return -exp(-q * T) * S * phi_d1 * v / 2 / sqrt(T) - r * K * exp(-r * T) * cdf(d2) + q * S * exp(-q * T) * cdf(d1);
}

template <class X>
X call_rho(X S, X T, X r, X q, X v, X K) {
	X d2 = (log(S / K) + (r - q - 0.5 * v * v) * T) / (v * sqrt(T));
	return K * T * exp(-r * T) * cdf(d2);
}

template <class X>
X call_epsilon(X S, X T, X r, X q, X v, X K) {
	X d1 = (log(S / K) + (r - q + 0.5 * v * v) * T) / (v * sqrt(T));
	X d2 = (log(S / K) + (r - q - 0.5 * v * v) * T) / (v * sqrt(T));
	X phi_d1 = M_2_SQRTPI / 2 * M_SQRT1_2 * exp(-d1 * d1 / 2);
	X phi_d2 = M_2_SQRTPI / 2 * M_SQRT1_2 * exp(-d2 * d2 / 2);
	return -T * S * exp(-q * T) * cdf(d1) - S * exp(-q * T) * phi_d1 * sqrt(T) / v + K * exp(-r * T) * phi_d2 * sqrt(T) / v;
}

template <class X>
X call_vega(X S, X T, X r, X q, X v, X K) {
	X d1 = (log(S / K) + (r - q + 0.5 * v * v) * T) / (v * sqrt(T));
	X phi_d1 = M_2_SQRTPI / 2 * M_SQRT1_2 * exp(-d1 * d1 / 2);
	return S*exp(-q*T)*phi_d1*sqrt(T);
}

template <class X>
X call_strike(X S, X T, X r, X q, X v, X K) {
	X d2 = (log(S / K) + (r - q - 0.5 * v * v) * T) / (v * sqrt(T));
	return -exp(-r * T) * cdf(d2);
}

static int test_greeks() {
	double S = 100;
	double K = 100;
	double r = 0.02;
	double q = 0.005;
	double v = 0.25;
	double T = 0.25;
	auto params = fms::multi_epsilon::add_epsilon({ S,T,r,q,v,K }, 2);
	auto S_ = params[0];
	auto T_ = params[1];
	auto r_ = params[2];
	auto q_ = params[3];
	auto v_ = params[4];
	auto K_ = params[5];
	auto res = bs_call<fms::multi_epsilon>(S_, T_, r_, q_, v_, K_);
	
	

	assert(fabs(res[fms::multi_epsilon::rep({ 0,0,0,0,0,0 }, 2)] - bs_call(S, T, r, q, v, K)) <= 2 * std::numeric_limits<double>::epsilon());
	assert(fabs(res[fms::multi_epsilon::rep({ 1,0,0,0,0,0 }, 2)] - call_delta(S, T, r, q, v, K)) <= 2 * std::numeric_limits<double>::epsilon());
	assert(fabs(res[fms::multi_epsilon::rep({ 0,1,0,0,0,0 }, 2)] + call_theta(S, T, r, q, v, K)) <= 8 * std::numeric_limits<double>::epsilon());
	assert(fabs(res[fms::multi_epsilon::rep({ 0,0,1,0,0,0 }, 2)] - call_rho(S, T, r, q, v, K)) <= 56 * std::numeric_limits<double>::epsilon());
	assert(fabs(res[fms::multi_epsilon::rep({ 0,0,0,1,0,0 }, 2)] - call_epsilon(S, T, r, q, v, K)) <= 0 * std::numeric_limits<double>::epsilon());
	assert(fabs(res[fms::multi_epsilon::rep({ 0,0,0,0,1,0 }, 2)] - call_vega(S, T, r, q, v, K)) <= 32 * std::numeric_limits<double>::epsilon());
	assert(fabs(res[fms::multi_epsilon::rep({ 0,0,0,0,0,1 }, 2)] - call_strike(S, T, r, q, v, K)) <= 1 * std::numeric_limits<double>::epsilon());

	
	return 0;
}
static int test_greeks_double = test_greeks();