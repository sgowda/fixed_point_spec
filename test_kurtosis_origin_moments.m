clear; clc; close all

stdev = 2^-2;
n_samples = 2^10;
n_trials = 100;
x = cplx_randn([n_samples, n_trials], stdev);

x_quant = round_inf_and_saturate_cplx_fi('', x, fi_dtype(1, 18, 17));
x_fi = fi(x_quant, 1, 18, 17);

T = 100;
for t=1:1
    x_re_fi = real(x_fi);
    x_im_fi = imag(x_fi);
    [x_re_fi, x_im_fi, x_sq_re, x_sq_im, abs_x_sq, abs_x_4th, x_3rd_re, ...
        x_3rd_im, origin_moment_error] = kurtosis_origin_moments_fi(x_re_fi(:,t), x_im_fi(:,t));

    [kappa_fi, num_fi, den_fi] = excess_kurtosis_complex_stream_fi(x_re_fi(:,t), x_im_fi(:,t))
    [kappa, num, den] = excess_kurtosis_complex_stream(x_quant(:,t))
end

origin_moment_error
