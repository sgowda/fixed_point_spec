
function [d, f, a, c, b, abs_mx_sq, e, abs_mean_x_sq] = kurtosis_cross_products_fi(x_mean_re, x_mean_im, x_sq_mean_re, x_sq_mean_im, abs_x_sq_mean, abs_x_4th_mean, x_3rd_mean_re, x_3rd_mean_im, varargin)
defaults = {'total_latency', 15, 'type_x', fi_dtype(1,18,17), 'acc_len', 14, 'logging', 1};
type_x = get_var('type_x', 'defaults', defaults, varargin{:});
acc_len = get_var('acc_len', 'defaults', defaults, varargin{:});
total_latency = get_var('total_latency', 'defaults', defaults, varargin{:});
logging = get_var('logging', 'defaults', defaults, varargin{:});
[m_x_dtype, x_sq_dtype, x_3rd_dtype, x_4th_dtype] = kurtosis_mean_types(type_x, acc_len);

% round inputs
[x_mean_re, x_mean_re_dtype]         = round_to_wordlength_fi('round_in1', x_mean_re, 25, m_x_dtype, 'logging', logging);
[x_mean_im, x_mean_im_dtype]         = round_to_wordlength_fi('round_in2', x_mean_im, 25, m_x_dtype, 'logging', logging);
[x_sq_mean_re, x_sq_mean_re_dtype ]  = round_to_wordlength_fi('round_in3', x_sq_mean_re, 35, x_sq_dtype, 'logging', logging);
[x_sq_mean_im, x_sq_mean_im_dtype ]  = round_to_wordlength_fi('round_in4', x_sq_mean_im, 35, x_sq_dtype, 'logging', logging);
[abs_x_sq_mean, abs_x_sq_mean_dtype] = round_to_wordlength_fi('round_in5', abs_x_sq_mean, 35, x_sq_dtype, 'logging', logging);
[x_3rd_mean_re, x_3rd_mean_re_dtype] = round_to_wordlength_fi('round_in6', x_3rd_mean_re, 35, x_3rd_dtype, 'logging', logging);
[x_3rd_mean_im, x_3rd_mean_im_dtype] = round_to_wordlength_fi('round_in7', x_3rd_mean_im, 35, x_3rd_dtype, 'logging', logging);

[x_mean_re_sq, x_mean_re_sq_dtype] = mult_fi('mult1', x_mean_re, x_mean_re, 'latency', 3, 'type_a', m_x_dtype, 'type_b', m_x_dtype);
[x_mean_im_sq, x_mean_im_sq_dtype] = mult_fi('mult2', x_mean_im, x_mean_im, 'latency', 3, 'type_a', m_x_dtype, 'type_b', m_x_dtype);
[mx_sq_re, mx_sq_re_dtype] = subtract_fi('add1', x_mean_re_sq, x_mean_im_sq, 'latency', 2, 'type_a', x_mean_re_sq_dtype, 'type_b', x_mean_im_sq_dtype);
[mx_re_times_mx_im, mx_re_times_mx_im_dtype] = mult_fi('mult3', x_mean_re, x_mean_im, 'latency', 5, 'type_a', m_x_dtype, 'type_b', m_x_dtype);
[mx_sq_im, mx_sq_im_dtype] = scale_fi('scale_mx_sq', mx_re_times_mx_im, 1, 'type_x', mx_re_times_mx_im_dtype);

[mx_sq_re_rounded, mx_sq_re_rounded_dtype] = round_to_wordlength_fi('round_mx_sq_re', mx_sq_re, 25, mx_sq_re_dtype, 'logging', logging);
[mx_sq_im_rounded, mx_sq_im_rounded_dtype] = round_to_wordlength_fi('round_mx_sq_im', mx_sq_im, 25, mx_sq_im_dtype, 'logging', logging);

x_sq_mean_re_del = delay_srl_fi('delay_sq1', x_sq_mean_re, 6);
x_sq_mean_im_del = delay_srl_fi('delay_sq2', x_sq_mean_im, 6);
[alpha, alpha_dtype] = mult_fi('mult4', x_sq_mean_re_del, mx_sq_re_rounded, 'latency', 3, 'type_a', x_sq_dtype, 'type_b', mx_sq_re_rounded_dtype);
[beta, beta_dtype] = mult_fi('mult5', x_sq_mean_im_del, mx_sq_im_rounded, 'latency', 3, 'type_a', x_sq_dtype, 'type_b', mx_sq_im_rounded_dtype);
[c_unscaled, c_unscaled_dtype] = add_fi('add3', alpha, beta, 'latency', 2, 'type_a', alpha_dtype, 'type_b', beta_dtype);

[gamma, gamma_dtype] = mult_fi('mult6', x_3rd_mean_re, x_mean_re, 'latency', 6, 'type_a', x_3rd_dtype, 'type_b', m_x_dtype);
[delta, delta_dtype] = mult_fi('mult7', x_3rd_mean_im, x_mean_im, 'latency', 6, 'type_a', x_3rd_dtype, 'type_b', m_x_dtype);
[b_unscaled, b_unscaled_dtype] = add_fi('add4', gamma, delta, 'latency', 2, 'type_a', gamma_dtype, 'type_b', delta_dtype);

[abs_m_x_sq, abs_m_x_sq_dtype] = add_fi('add2', x_mean_re_sq, x_mean_im_sq, 'latency', 2, 'type_a', x_mean_re_sq_dtype, 'type_b', x_mean_im_sq_dtype);
[abs_m_x_sq_25bit, abs_m_x_sq_25bit_dtype] = round_to_wordlength_fi('conv_25bit', abs_m_x_sq, 25, abs_m_x_sq_dtype, 'latency', 3, 'logging', logging);
[abs_m_x_sq_35bit, abs_m_x_sq_35bit_dtype] = round_to_wordlength_fi('conv_35bit', abs_m_x_sq, 35, abs_m_x_sq_dtype, 'latency', 3, 'logging', logging);
[e_adv, e_dtype] = mult_fi('mult8', abs_m_x_sq_25bit, abs_m_x_sq_35bit, 'latency', 5, 'type_a', abs_m_x_sq_25bit_dtype, 'type_b', abs_m_x_sq_35bit_dtype);

% calc fourth moment terms
[b_adv, b_dtype] = scale_fi('Scale', b_unscaled, 2, 'type_x', b_unscaled_dtype);
b_sig = delay_srl_fi('delay_b', b_adv, 7);

[c_adv, c_dtype] = scale_fi('Scale1', c_unscaled, 1, 'type_x', c_unscaled_dtype);
c_sig = delay_srl_fi('c_del', c_adv, 4);

% TODO check latency!
[abs_m_x_sq_rounded, abs_m_x_sq_rounded_dtype] = round_to_wordlength_fi('round_abs_mx_sq', abs_m_x_sq, 25, abs_m_x_sq_dtype, 'latency', 5, 'logging', logging); 
abs_x_sq_mean_del = delay_srl_fi('del1', abs_x_sq_mean, 10);
[d_unscaled, d_unscaled_dtype] = mult_fi('Mult4',  abs_x_sq_mean_del, abs_m_x_sq_rounded, 'latency', 5, 'type_a', x_sq_dtype, 'type_b', abs_m_x_sq_rounded_dtype);
[d_sig, d_dtype] = scale_fi('Scale3', d_unscaled, 2, 'type_x', d_unscaled_dtype);
% d_sig = delay_srl('delay_sq3', d_adv, 5);

e_sig = delay_srl_fi('del_e', e_adv, 2);

[f_adv, f_dtype] = scale_fi('Scale4', e_adv, 2, 'type_x', e_dtype);
f_sig = delay_srl_fi('del_f', f_adv, 2);

% third term in complex kurtosis
add_latency = 2;
mx_sq_re_del = delay_srl_fi('del2', mx_sq_re, 1);
mx_sq_im_del = delay_srl_fi('del3', mx_sq_im, 1);
[mx_sq_re_rounded_mb, mx_sq_re_rounded_mb_dtype] = round_to_wordlength_fi('conv_mx_sq_re', mx_sq_re_del, 35, mx_sq_re_dtype, 'latency', 0, 'logging', logging);
[mx_sq_im_rounded_mb, mx_sq_im_rounded_mb_dtype] = round_to_wordlength_fi('conv_mx_sq_im', mx_sq_im_del, 35, mx_sq_im_dtype, 'latency', 0, 'logging', logging);
[abs_mean_x_re, abs_mean_x_im, abs_mean_x_dtype] = cplx_sub_fi('cplx_sub', {x_sq_mean_re, x_sq_mean_im}, ...
    {mx_sq_re_rounded_mb, mx_sq_im_rounded_mb}, 'latency', add_latency, 'full_precision', 1, 'type_a', x_sq_mean_re_dtype, 'type_b', mx_sq_im_rounded_mb_dtype);

% mean-squared for denominator
[abs_mean_x_re_sq, abs_mean_x_re_sq_dtype] = mult_fi('mult9', abs_mean_x_re, abs_mean_x_re, 'latency', 4, 'type_a', abs_mean_x_dtype, 'type_b', abs_mean_x_dtype);
[abs_mean_x_im_sq, abs_mean_x_im_sq_dtype] = mult_fi('mult10', abs_mean_x_im, abs_mean_x_im, 'latency', 4, 'type_a', abs_mean_x_dtype, 'type_b', abs_mean_x_dtype);
[abs_mean_x_sq_sig, abs_mean_x_sq_type] = add_fi('add5', abs_mean_x_re_sq, abs_mean_x_im_sq, 'latency', 3, 'type_a', abs_mean_x_re_sq_dtype, 'type_b', abs_mean_x_im_sq_dtype);

a = abs_x_4th_mean;
b = b_sig;
c = c_sig;
d = d_sig;
e = e_sig;
f = f_sig;
abs_mean_x_sq = abs_mean_x_sq_sig;
abs_mx_sq = abs_m_x_sq;

end
