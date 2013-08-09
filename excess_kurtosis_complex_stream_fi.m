function [kappa_x, kurtosis_num, second_abs_moment_sq, fourth_abs_moment, ...
    mean_power] = excess_kurtosis_complex_stream_fi(x_re_fi, x_im_fi)
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

% x_acc_re = accumulate_fi(x_re_fi, single_acc_dtype);
% x_acc_im = accumulate_fi(x_im_fi, single_acc_dtype);
% x_sq_acc_re = accumulate_fi(x_sq_re, double_acc_dtype);
% x_sq_acc_im = accumulate_fi(x_sq_im, double_acc_dtype);
% abs_x_sq_acc = accumulate_fi(abs_x_sq, double_acc_dtype);
% abs_x_4th_acc = accumulate_fi(abs_x_4th, quad_acc_dtype);
% x_3rd_acc_re = accumulate_fi(x_3rd_re, triple_acc_dtype);
% x_3rd_acc_im = accumulate_fi(x_3rd_im, triple_acc_dtype);

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

% fourth abs moment
m_x_re = scale_fi('', x_acc_re, -acc_len);
m_x_im = scale_fi('', x_acc_im, -acc_len);


round = 1;
if round
    m_x_re_rounded       = quantize(single_acc_quantizer, m_x_re);
    m_x_im_rounded       = quantize(single_acc_quantizer, m_x_im);
    x_sq_acc_re_rounded  = quantize(double_acc_quantizer, x_sq_acc_re);
    x_sq_acc_im_rounded  = quantize(double_acc_quantizer, x_sq_acc_im);
    abs_x_sq_acc_rounded = quantize(double_acc_quantizer, abs_x_sq_acc);
    x_3rd_acc_re_rounded = quantize(triple_acc_quantizer, x_3rd_acc_re);
    x_3rd_acc_im_rounded = quantize(triple_acc_quantizer, x_3rd_acc_im);
else
    m_x_re_rounded       = m_x_re;
    m_x_im_rounded       = m_x_im;
    x_sq_acc_re_rounded  = x_sq_acc_re;
    x_sq_acc_im_rounded  = x_sq_acc_im;
    abs_x_sq_acc_rounded = abs_x_sq_acc;
    x_3rd_acc_re_rounded = x_3rd_acc_re;
    x_3rd_acc_im_rounded = x_3rd_acc_im;
end

a = round_inf_and_saturate_fi('', abs_x_4th_acc, fi_dtype(1, 74+acc_len, 68)); % 1e-14 unrounded error

x_3rd_acc = double(x_3rd_acc_re_rounded) + 1j*double(x_3rd_acc_im_rounded);
b_unsc = x_3rd_acc_re_rounded*m_x_re_rounded + x_3rd_acc_im_rounded*m_x_im_rounded;
b = scale_fi('', b_unsc, 2); % 1e-18 unrounded error

m_x_sq_re = m_x_re_rounded*m_x_re_rounded - m_x_im_rounded*m_x_im_rounded;
m_x_sq_im = 2*m_x_re_rounded*m_x_im_rounded;
m_x_sq_re_rounded = round_inf_and_saturate_fi('', m_x_sq_re, fi_dtype(1, 25, 22));
m_x_sq_im_rounded = round_inf_and_saturate_fi('', m_x_sq_im, fi_dtype(1, 25, 22));

m_x_sq_re_error = double(m_x_sq_re_rounded) - double(m_x_sq_re);
m_x_sq_im_error = double(m_x_sq_im_rounded) - double(m_x_sq_im);
x_sq_acc_re_error = double(x_sq_acc_re_rounded) - double(x_sq_acc_re);
x_sq_acc_im_error = double(x_sq_acc_im_rounded) - double(x_sq_acc_im);

% NOTE: rounding |m_x|^2 too much has a significant effect on performance!
c_unsc = x_sq_acc_re_rounded*m_x_sq_re_rounded + x_sq_acc_im_rounded*m_x_sq_im_rounded;
% c_unsc = x_sq_acc_re_rounded*m_x_sq_re_rounded + x_sq_acc_im_rounded*m_x_sq_im_rounded;
% c_unsc = x_sq_acc_re_rounded*m_x_sq_re_rounded + x_sq_acc_im_rounded*m_x_sq_im_rounded;
c = scale_fi('', c_unsc, 1); % 1e-19 unrounded error

% TODO some of these multipliers are duplicated!
abs_m_x_sq = m_x_re_rounded*m_x_re_rounded + m_x_im_rounded*m_x_im_rounded; 
abs_m_x_sq_rounded = round_inf_and_saturate_fi('', abs_m_x_sq, fi_dtype(1, 25, 22));
d_unsc = abs_x_sq_acc_rounded * abs_m_x_sq_rounded;
d = scale_fi('', d_unsc, 2); % 1e-18 unrounded error


e = square_35x25_fi('', abs_m_x_sq, 5); % TODO fewer integer bits?!

% e = abs_m_x_sq_rounded * abs_m_x_sq_rounded; % 1e-24 unrounded error
f = scale_fi('', e, 2); % 1e-24 unrounded error

fourth_abs_moment = scale_fi('', (a - b) + (c + d), -acc_len) + (e - f);

x = double(x_re_fi) + 1j*double(x_im_fi);
m_x = mean(x);
z = x - m_x;
E_abs_z_4 = mean((abs(z)).^4);
fourth_abs_moment_error = E_abs_z_4 - fourth_abs_moment


a_fl = sum(abs(x).^4);
b_fl = 4*real( sum(x .* (abs(x).^2)) .* conj(m_x) );
c_fl = 2*real( sum(x.^2) * conj(m_x^2) );
c_fl = 2*(real(sum(x.^2))*real(m_x^2) + imag(sum(x.^2))*imag(m_x^2));
d_fl = 4*sum( abs(x).^2 * abs(m_x)^2 );
e_fl = abs(m_x)^4;
f_fl = 4*abs(m_x)^4;
fourth_abs_moment_fl = (2^-acc_len)*(a_fl - b_fl + c_fl + d_fl) + (e_fl - f_fl);

if verbose

    fourth_abs_moment_error_fl = E_abs_z_4 - fourth_abs_moment_fl
end

if verbose
    fprintf('a_error\n', a_fl - double(a));
    fprintf('b_error\n', b_fl - double(b));
    fprintf('c_error\n', c_fl - double(c));
    fprintf('d_error\n', d_fl - double(d));
    fprintf('e_error\n', e_fl - double(e));
    fprintf('f_error\n', f_fl - double(f));
end

% second abs moment
sum_power = round_inf_and_saturate_fi('', abs_x_sq_acc, fi_dtype(1, 35, 17));
mean_power = scale_fi('', sum_power, -acc_len);

second_abs_moment = mean_power - abs_m_x_sq;
second_abs_moment_sq = square_35x25_fi('', second_abs_moment, 5);
% second_abs_moment_25bit = round_inf_and_saturate_fi('', second_abs_moment, fi_dtype(1, 25, 20));
% second_abs_moment_35bit = round_inf_and_saturate_fi('', second_abs_moment, fi_dtype(1, 35, 30));
% second_abs_moment_sq = second_abs_moment_35bit * second_abs_moment_25bit;
second_abs_moment_sq_fl = (mean(abs(x).^2) - abs(m_x).^2).^2;
second_moment_sq_error = second_abs_moment_sq - second_abs_moment_sq_fl


% third term
m_x_sq_re_rounded2 = round_inf_and_saturate_fi('', m_x_sq_re, fi_dtype(1, 34, 31));
m_x_sq_im_rounded2 = round_inf_and_saturate_fi('', m_x_sq_im, fi_dtype(1, 34, 31));
third_term_sqrt_re = scale_fi('', x_sq_acc_re_rounded, -acc_len) - m_x_sq_re_rounded2;
third_term_sqrt_im = scale_fi('', x_sq_acc_im_rounded, -acc_len) - m_x_sq_im_rounded2;
third_term = third_term_sqrt_re*third_term_sqrt_re + third_term_sqrt_im*third_term_sqrt_im;
m_x = mean(x);
third_term_fl = abs( mean(x.^2) - m_x.^2 ).^2;
third_term_error = third_term_fl - third_term

kurtosis_num = double(fourth_abs_moment - third_term);
kappa_x = kurtosis_num/double(second_abs_moment_sq) - 2;
