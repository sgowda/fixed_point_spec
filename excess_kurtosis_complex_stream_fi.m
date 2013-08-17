function [kappa_x, num, den, mean_power] = ...
    excess_kurtosis_complex_stream_fi2(x_re_fi, x_im_fi)
% excess_kurtosis_complex
%
%   kappa_x = excess_kurtosis_complex(x)
%
%   Calculates the complex random variable version of the excess kurtosis
%   of a set of sample data. The set of sample data can be either real or
%   complex. Here, the biased estimate of the kurtosis is calculated.
%
%   Inputs:
%       x           =   vector of input sample data (real or complex)
%
%   Outputs:
%       kappa_x     =   biased sample excess kurtosis of the vector x
%
verbose = 0;

% Configure
acc_len = log2(length(x_re_fi));
type_x = fi_dtype(1, 18, 17);

% Data types of cross products
single_type_x = type_x;
double_type_x = type_x^2 + type_x^2;
triple_type_x = type_x^3;
quad_type_x = double_type_x^2;

% Data types of accumulator outputs
single_type_x_acc = fi_dtype(1, single_type_x.WordLength+acc_len, single_type_x.FractionLength+acc_len);
double_type_x_acc = fi_dtype(1, double_type_x.WordLength+acc_len, double_type_x.FractionLength);
triple_type_x_acc = fi_dtype(1, triple_type_x.WordLength+acc_len, triple_type_x.FractionLength);
quad_type_x_acc = fi_dtype(1, quad_type_x.WordLength+acc_len, quad_type_x.FractionLength);

%-- Calculate cross products: 5 reg multipliers, 1 large multiplier, and 2 adders
x_re_sq = x_re_fi.*x_re_fi;
x_im_sq = x_im_fi.*x_im_fi;

abs_x_sq = x_re_sq + x_im_sq;
abs_x_4th = abs_x_sq .* abs_x_sq;

x_sq_re = x_re_sq - x_im_sq; % TODO reuse hardware with the abs_x_sq computation!
x_sq_im_unscaled = x_re_fi.*x_im_fi;
x_sq_im = scale_fi('', x_sq_im_unscaled, 1);

x_3rd_re = x_re_fi.*abs_x_sq;
x_3rd_im = x_im_fi.*abs_x_sq;

%-- Simulate accumulators
single_acc_dtype = fi_dtype(1, 96, single_type_x.FractionLength);
double_acc_dtype = fi_dtype(1, 96, double_type_x.FractionLength);
triple_acc_dtype = fi_dtype(1, 96, triple_type_x.FractionLength);
quad_acc_dtype = fi_dtype(1, 96, quad_type_x.FractionLength);

x_acc_re = sum(x_re_fi);
x_acc_im = sum(x_im_fi);
x_sq_acc_re = sum(x_sq_re);
x_sq_acc_im = sum(x_sq_im);
abs_x_sq_acc = sum(abs_x_sq);
abs_x_4th_acc = sum(abs_x_4th);
x_3rd_acc_re = sum(x_3rd_re);
x_3rd_acc_im = sum(x_3rd_im);

m_x_re = scale_fi('', x_acc_re, -acc_len);
m_x_im = scale_fi('', x_acc_im, -acc_len);

% Simulate rounding
single_acc_quantizer = fixed.Quantizer('WordLength', 25, ...
    'FractionLength', 25-(type_x.WordLength - type_x.FractionLength), ...
    'RoundingMethod', 'round', 'OverflowAction', 'saturate');

double_acc_quantizer = fixed.Quantizer('WordLength', 35, ...
    'FractionLength', 35-(double_type_x_acc.WordLength - double_type_x_acc.FractionLength), ...
    'RoundingMethod', 'round', 'OverflowAction', 'saturate');

triple_acc_quantizer = fixed.Quantizer('WordLength', 35, ...
    'FractionLength', 35-(triple_type_x_acc.WordLength - triple_type_x_acc.FractionLength), ...
    'RoundingMethod', 'round', 'OverflowAction', 'saturate');

m_x_re_rounded       = quantize(single_acc_quantizer, m_x_re);
m_x_im_rounded       = quantize(single_acc_quantizer, m_x_im);
x_sq_acc_re_rounded  = quantize(double_acc_quantizer, x_sq_acc_re);
x_sq_acc_im_rounded  = quantize(double_acc_quantizer, x_sq_acc_im);
abs_x_sq_acc_rounded = quantize(double_acc_quantizer, abs_x_sq_acc);
x_3rd_acc_re_rounded = quantize(triple_acc_quantizer, x_3rd_acc_re);
x_3rd_acc_im_rounded = quantize(triple_acc_quantizer, x_3rd_acc_im);

[num, den, mean_power] = kurtosis_moment_calc_fi(m_x_re_rounded, m_x_im_rounded, x_sq_acc_re_rounded, x_sq_acc_im_rounded, abs_x_sq_acc_rounded, abs_x_4th_acc, x_3rd_acc_re_rounded, x_3rd_acc_im_rounded, 'acc_len', acc_len);
kappa_x = double(num)/double(den);