clear; clc; close all

stdev = 2^-2;
n_samples = 2^16;
n_trials = 100;
x = cplx_randn([n_samples, n_trials], stdev);

x_quant = round_inf_and_saturate_cplx_fi('', x, fi_dtype(1, 18, 17), 'logging', 0);
x_fi = fi(x_quant, 1, 18, 17);

T = 100;
for t=1:T
    t
    x_re_fi = real(x_fi(:,t));
    x_im_fi = imag(x_fi(:,t));
%     [x_re_fi, x_im_fi, x_sq_re, x_sq_im, abs_x_sq, abs_x_4th, x_3rd_re, ...
%         x_3rd_im, origin_moment_error] = kurtosis_origin_moments_fi(x_re_fi, x_im_fi);

    [kappa_fi(t), num_fi(t), den_fi(t)] = excess_kurtosis_complex_stream_fi(x_re_fi, x_im_fi);
    [kappa(t), num(t), den(t)] = excess_kurtosis_complex_stream(x_quant(:,t));
end

num_error = double(num_fi) - num;
den_error = double(den_fi) - den;