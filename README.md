﻿# epsilon

This library can be used to compute derivatives of arbitrary order to machine precision.
The class `epsilon` implements Toeplitz matrices together with standard numerical operations on them.

If ε is not 0 and ε<sup>2</sup> = 0 then for any sufficiently differentiable function _f_,
_f_(_x_ + ε) = _f_(_x_) + _f_'(_x_) ε + _f_''(_x_) ε<sup>2</sup>/2 +  =  _f_(_x_) + _f_'(_x_) ε.
The coefficient of ε is the derivative of _f_.

For example, if _f_(_x_) = _x_<sup>2</sup> 
then (_x_ + ε)<sup>2</sup> = _x_<sup>2 + 2_x_ε + ε<sup>2</sup> = _x_<sup>2 + 2_x_ε,
so the derivative of _x_<sup>2</sup> is 2_x_. No need to take limits of difference
quotients.

In order to use this library you must define functions that take generic arguments. For example:
```
template<class X>
X f(X x)
{
    return x*x + 2*x + 3;
}
```
Now you can use `epsilon` to compute derivatives.
```
double x = 1.23;
auto y = f(epsilon<3>(x));
assert (y[0] == f(x));
assert (y[1] == 2*x + 2);
assert (y[2] == 2);
```

