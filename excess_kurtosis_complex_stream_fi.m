function [kappa_x, num, den, mean_power] = ...
    excess_kurtosis_complex_stream_fi(x_re_fi, x_im_fi)
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

% Calculate origin moments
[x_re_fi, x_im_fi, x_sq_re, x_sq_im, abs_x_sq, abs_x_4th, x_3rd_re, x_3rd_im] = kurtosis_origin_moments_fi(x_re_fi, x_im_fi);

%-- Simulate accumulators
[x_acc_re, x_acc_re_dtype]           = accumulate_fi(x_re_fi);
[x_acc_im, x_acc_im_dtype]           = accumulate_fi(x_im_fi);
[x_sq_acc_re, x_sq_acc_re_dtype]     = accumulate_fi(x_sq_re);
[x_sq_acc_im, x_sq_acc_im_dtype]     = accumulate_fi(x_sq_im);
[abs_x_sq_acc, abs_x_sq_acc_dtype]   = accumulate_fi(abs_x_sq);
[abs_x_4th_acc, abs_x_4th_acc_dtype] = accumulate_fi(abs_x_4th);
[x_3rd_acc_re, x_3rd_acc_re_dtype]   = accumulate_fi(x_3rd_re);
[x_3rd_acc_im, x_3rd_acc_im_dtype]   = accumulate_fi(x_3rd_im);

[x_mean_re, x_mean_re_dtype]           = scale_fi('scale1', x_acc_re, -acc_len);
[x_mean_im, x_mean_im_dtype]           = scale_fi('scale2', x_acc_im, -acc_len);
[x_sq_mean_re, x_sq_mean_re_dtype]     = scale_fi('scale3', x_sq_acc_re, -acc_len);
[x_sq_mean_im, x_sq_mean_im_dtype]     = scale_fi('scale4', x_sq_acc_im, -acc_len);
[abs_x_sq_mean, abs_x_sq_mean_dtype]   = scale_fi('scale5', abs_x_sq_acc, -acc_len);
[abs_x_4th_mean, abs_x_4th_mean_dtype] = scale_fi('scale6', abs_x_4th_acc, -acc_len);
[x_3rd_mean_re, x_3rd_mean_re_dtype]   = scale_fi('scale7', x_3rd_acc_re, -acc_len);
[x_3rd_mean_im, x_3rd_mean_im_dtype]   = scale_fi('scale8', x_3rd_acc_im, -acc_len);

% calculate central moments from origin moments
[num, den, mean_power] = kurtosis_moment_calc_norm_inputs_fi(x_mean_re, x_mean_im, x_sq_mean_re, x_sq_mean_im, abs_x_sq_mean, abs_x_4th_mean, x_3rd_mean_re, x_3rd_mean_im);

% use floating-point division to get the kurtosis ratio (not used in hardware)
kappa_x = double(num)/double(den) - 2;
