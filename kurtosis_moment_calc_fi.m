function [num, den, abs_X_sq_mean] = kurtosis_moment_calc_fi(m_x_re, m_x_im, X_sq_acc_re, X_sq_acc_im, abs_X_sq_acc, abs_X_4th_acc, X_3rd_acc_re, X_3rd_acc_im, varargin)
defaults = {'acc_len', 14, 'type_x', fi_dtype(1, 18, 17)};
acc_len = get_var('acc_len', 'defaults', defaults, varargin{:});
type_x = get_var('type_x', 'defaults', defaults, varargin{:});
[m_x_type, x_sq_type, x_3rd_type, x_4th_type] = kurtosis_acc_rounding_types(type_x, acc_len);

bit_width_power = 32;
bit_width_num = 96;
bit_width_den = 64;

rounding_latency = 1;
add_latency_96bit = 4;
cross_product_latency = 15;
adder_tree_latency = add_latency_96bit * 3 + rounding_latency;
abs_X_sq_acc_del = delay_srl_fi('cross_product_del', abs_X_sq_acc, cross_product_latency);

abs_X_sq_acc_del2 = round_to_wordlength_fi('adder_tree_del', abs_X_sq_acc_del, bit_width_power, x_sq_type, 'latency', adder_tree_latency);
abs_X_sq_mean = scale_fi('scale', abs_X_sq_acc_del2, -acc_len);

% cross-products
[d, f, a, c, b, abs_mx_sq, e, h] = kurtosis_cross_products(m_x_re, m_x_im, X_sq_acc_re, X_sq_acc_im, abs_X_sq_acc, abs_X_4th_acc, X_3rd_acc_re, X_3rd_acc_im, 'acc_len', acc_len, 'total_latency', cross_product_latency, 'type_x', type_x);

% cross-product data types
[a_dtype, b_dtype, c_dtype, d_dtype, e_dtype, f_dtype, abs_m_x_sq_dtype, abs_mean_x_sq_dtype] = ...
    kurtosis_cross_product_dtypes(m_x_type, m_x_type, ...
    x_sq_type, x_sq_type, x_sq_type, x_4th_type, x_3rd_type, x_3rd_type, 'acc_len', acc_len);

% denominator
den = kurtosis_den(abs_X_sq_acc_del, abs_mx_sq, 'acc_len', acc_len, 'bit_width', bit_width_den, 'abs_m_x_sq_dtype', abs_m_x_sq_dtype, 'abs_X_sq_acc_del_dtype', x_sq_type);

% numerator
num = kurtosis_num(d, f, a, c, b, e, h, 'acc_len', acc_len, 'bit_width', bit_width_num, 'a_dtype', a_dtype, ...
        'b_dtype', b_dtype, 'c_dtype', c_dtype, 'd_dtype', d_dtype, ...
        'e_dtype', e_dtype, 'f_dtype', f_dtype, 'h_dtype', abs_mean_x_sq_dtype);
    
end

function [d, f, a, c, b, abs_mx_sq, e, abs_mean_x_sq] = kurtosis_cross_products(m_x_re, m_x_im, X_sq_acc_re, X_sq_acc_im, abs_X_sq_acc, abs_X_4th_acc, X_3rd_acc_re, X_3rd_acc_im, varargin)
defaults = {'acc_len', 14, 'total_latency', 15, 'type_x', fi_dtype(1,18,17)};
acc_len = get_var('acc_len', 'defaults', defaults, varargin{:});
type_x = get_var('type_x', 'defaults', defaults, varargin{:});
total_latency = get_var('total_latency', 'defaults', defaults, varargin{:});
[m_x_dtype, x_sq_dtype, x_3_dtype, x_4th_type] = kurtosis_acc_rounding_types(type_x, acc_len);

[m_x_re_sq, m_x_re_sq_dtype] = mult_fi('mult1', m_x_re, m_x_re, 'latency', 3, 'type_a', m_x_dtype, 'type_b', m_x_dtype);
[m_x_im_sq, m_x_im_sq_dtype] = mult_fi('mult2', m_x_im, m_x_im, 'latency', 3, 'type_a', m_x_dtype, 'type_b', m_x_dtype);
[mx_sq_re, mx_sq_re_dtype] = subtract_fi('add1', m_x_re_sq, m_x_im_sq, 'latency', 2, 'type_a', m_x_re_sq_dtype, 'type_b', m_x_im_sq_dtype);
[mx_re_times_mx_im, mx_re_times_mx_im_dtype] = mult_fi('mult3', m_x_re, m_x_im, 'latency', 5, 'type_a', m_x_dtype, 'type_b', m_x_dtype);
[mx_sq_im, mx_sq_im_dtype] = scale_fi('scale_mx_sq', mx_re_times_mx_im, 1, 'type_x', mx_re_times_mx_im_dtype);

[mx_sq_re_rounded, mx_sq_re_rounded_dtype] = round_to_wordlength_fi('Convert', mx_sq_re, 25, mx_sq_re_dtype);
[mx_sq_im_rounded, mx_sq_im_rounded_dtype] = round_to_wordlength_fi('Convert1', mx_sq_im, 25, mx_sq_im_dtype);

X_sq_acc_re_del = delay_srl_fi('delay_sq1', X_sq_acc_re, 6);
X_sq_acc_im_del = delay_srl_fi('delay_sq2', X_sq_acc_im, 6);
[alpha, alpha_dtype] = mult_fi('mult4', X_sq_acc_re_del, mx_sq_re_rounded, 'latency', 3, 'type_a', x_sq_dtype, 'type_b', mx_sq_re_rounded_dtype);
[beta, beta_dtype] = mult_fi('mult5', X_sq_acc_im_del, mx_sq_im_rounded, 'latency', 3, 'type_a', x_sq_dtype, 'type_b', mx_sq_im_rounded_dtype);
[c_unscaled, c_unscaled_dtype] = add_fi('add3', alpha, beta, 'latency', 2, 'type_a', alpha_dtype, 'type_b', beta_dtype);

[gamma, gamma_dtype] = mult_fi('mult6', X_3rd_acc_re, m_x_re, 'latency', 6, 'type_a', x_3_dtype, 'type_b', m_x_dtype);
[delta, delta_dtype] = mult_fi('mult7', X_3rd_acc_im, m_x_im, 'latency', 6, 'type_a', x_3_dtype, 'type_b', m_x_dtype);
[b_unscaled, b_unscaled_dtype] = add_fi('add4', gamma, delta, 'latency', 2, 'type_a', gamma_dtype, 'type_b', delta_dtype);

[abs_m_x_sq, abs_m_x_sq_dtype] = add_fi('add2', m_x_re_sq, m_x_im_sq, 'latency', 2, 'type_a', m_x_re_sq_dtype, 'type_b', m_x_im_sq_dtype);
[abs_m_x_sq_25bit, abs_m_x_sq_25bit_dtype] = round_to_wordlength_fi('conv_25bit', abs_m_x_sq, 25, abs_m_x_sq_dtype, 'latency', 3);
[abs_m_x_sq_35bit, abs_m_x_sq_35bit_dtype] = round_to_wordlength_fi('conv_35bit', abs_m_x_sq, 35, abs_m_x_sq_dtype, 'latency', 3);
[e_adv, e_dtype] = mult_fi('mult8', abs_m_x_sq_25bit, abs_m_x_sq_35bit, 'latency', 5, 'type_a', abs_m_x_sq_25bit_dtype, 'type_b', abs_m_x_sq_35bit_dtype);

% calc fourth moment terms
[b_adv, b_dtype] = scale_fi('Scale', b_unscaled, 2, 'type_x', b_unscaled_dtype);
b_sig = delay_srl_fi('delay_b', b_adv, 7);

[c_adv, c_dtype] = scale_fi('Scale1', c_unscaled, 1, 'type_x', c_unscaled_dtype);
c_sig = delay_srl_fi('c_del', c_adv, 4);

% TODO check latency!
[abs_m_x_sq_rounded, abs_m_x_sq_rounded_dtype] = round_to_wordlength_fi('Convert5', abs_m_x_sq, 25, abs_m_x_sq_dtype, 'latency', 5); 
abs_X_sq_acc_del = delay_srl_fi('del1', abs_X_sq_acc, 10);
[d_unscaled, d_unscaled_dtype] = mult_fi('Mult4',  abs_X_sq_acc_del, abs_m_x_sq_rounded, 'latency', 5, 'type_a', x_sq_dtype, 'type_b', abs_m_x_sq_rounded_dtype);
[d_sig, d_dtype] = scale_fi('Scale3', d_unscaled, 2, 'type_x', d_unscaled_dtype);
% d_sig = delay_srl('delay_sq3', d_adv, 5);

e_sig = delay_srl_fi('del_e', e_adv, 2);

[f_adv, f_dtype] = scale_fi('Scale4', e_adv, 2, 'type_x', e_dtype);
f_sig = delay_srl_fi('del_f', f_adv, 2);

% third term in complex kurtosis
add_latency = 2;
mx_sq_re_del = delay_srl_fi('del2', mx_sq_re, 1);
mx_sq_im_del = delay_srl_fi('del3', mx_sq_im, 1);
[mx_sq_re_rounded_mb, mx_sq_re_rounded_mb_dtype] = round_to_wordlength_fi('conv_mx_sq_re', mx_sq_re_del, 35, mx_sq_re_dtype, 'latency', 0);
[mx_sq_im_rounded_mb, mx_sq_im_rounded_mb_dtype] = round_to_wordlength_fi('conv_mx_sq_im', mx_sq_im_del, 35, mx_sq_im_dtype, 'latency', 0);
[X_sq_mean_re, X_sq_mean_re_dtype] = scale_fi('scale_m2_re', X_sq_acc_re_del, -acc_len, 'type_x', x_sq_dtype);
[X_sq_mean_im, X_sq_mean_im_dtype] = scale_fi('scale_m2_im', X_sq_acc_im_del, -acc_len, 'type_x', x_sq_dtype);
[abs_mean_x_re, abs_mean_x_im, abs_mean_x_dtype] = cplx_sub_fi('cplx_sub', {X_sq_mean_re, X_sq_mean_im}, ...
    {mx_sq_re_rounded_mb, mx_sq_im_rounded_mb}, 'latency', add_latency, 'full_precision', 1, 'type_a', X_sq_mean_re_dtype, 'type_b', mx_sq_im_rounded_mb_dtype);

% mean-squared for denominator
[abs_mean_x_re_sq, abs_mean_x_re_sq_dtype] = mult_fi('mult9', abs_mean_x_re, abs_mean_x_re, 'latency', 4, 'type_a', abs_mean_x_dtype, 'type_b', abs_mean_x_dtype);
[abs_mean_x_im_sq, abs_mean_x_im_sq_dtype] = mult_fi('mult10', abs_mean_x_im, abs_mean_x_im, 'latency', 4, 'type_a', abs_mean_x_dtype, 'type_b', abs_mean_x_dtype);
[abs_mean_x_sq_sig, abs_mean_x_sq_type] = add_fi('add5', abs_mean_x_re_sq, abs_mean_x_im_sq, 'latency', 3, 'type_a', abs_mean_x_re_sq_dtype, 'type_b', abs_mean_x_im_sq_dtype);

a = abs_X_4th_acc;
b = b_sig;
c = c_sig;
d = d_sig;
e = e_sig;
f = f_sig;
abs_mean_x_sq = abs_mean_x_sq_sig;
abs_mx_sq = abs_m_x_sq;

end


function [den] = kurtosis_den(sum_abs_x_sq, abs_m_x_sq, varargin)
defaults = {'acc_len', 14, 'bit_width', 64, 'abs_m_x_sq_dtype'};
acc_len = get_var('acc_len', 'defaults', defaults, varargin{:});
bit_width = get_var('bit_width', 'defaults', defaults, varargin{:});
abs_m_x_sq_dtype = get_var('abs_m_x_sq_dtype', 'defaults', defaults, varargin{:});
abs_X_sq_acc_del_dtype = get_var('abs_X_sq_acc_del_dtype', 'defaults', defaults, varargin{:});

[mean_abs_x_sq, mean_abs_x_sq_dtype] = scale_fi('scale', sum_abs_x_sq, -acc_len, 'type_x', abs_X_sq_acc_del_dtype);
[sec_moment_unr, sec_moment_unr_dtype] = subtract_fi('Sub', mean_abs_x_sq, abs_m_x_sq, 'latency', 4, 'type_a', mean_abs_x_sq_dtype, 'type_b', abs_m_x_sq_dtype);
[den_unr, den_unr_dtype] = mult_35x25_fi('square', sec_moment_unr, sec_moment_unr, 'type_a', sec_moment_unr_dtype, 'type_b', sec_moment_unr_dtype);
[den_sig, den_dtype] = round_to_wordlength_fi('conv', den_unr, bit_width, den_unr_dtype);
den = den_sig;

end

function num = kurtosis_num(d, f, a, c, b, e, h, varargin)
defaults = {'acc_len', 14, 'bit_width', 96, 'conv_latency', 1};
acc_len = get_var('acc_len', 'defaults', defaults, varargin{:});
bit_width = get_var('bit_width', 'defaults', defaults, varargin{:});
conv_latency = get_var('conv_latency', 'defaults', defaults, varargin{:});
a_dtype = get_var('a_dtype', 'defaults', defaults, varargin{:});
b_dtype = get_var('b_dtype', 'defaults', defaults, varargin{:});
c_dtype = get_var('c_dtype', 'defaults', defaults, varargin{:});
d_dtype = get_var('d_dtype', 'defaults', defaults, varargin{:});
e_dtype = get_var('e_dtype', 'defaults', defaults, varargin{:});
f_dtype = get_var('f_dtype', 'defaults', defaults, varargin{:});
h_dtype = get_var('h_dtype', 'defaults', defaults, varargin{:});

add_latency = 4;

% scale factors
[a_scale, a_dtype] = scale_fi('scale_a', a, -acc_len, 'type_x', a_dtype);
[b_scale, b_dtype] = scale_fi('scale_b', b, -acc_len, 'type_x', b_dtype);
[c_scale, c_dtype] = scale_fi('scale_c', c, -acc_len, 'type_x', c_dtype);
[d_scale, d_dtype] = scale_fi('scale_d', d, -acc_len, 'type_x', d_dtype);

% adder layer 1
[e_minus_f, e_minus_f_dtype] = subtract_fi('sub_ef', e, f, 'latency', add_latency, 'full_precision', 1, 'type_a', e_dtype, 'type_b', f_dtype);
[c_minus_b, c_minus_b_dtype] = subtract_fi('sub_cb', c_scale, b_scale, 'latency', add_latency, 'full_precision', 1, 'type_a', c_dtype, 'type_b', b_dtype);
[d_minus_h, d_minus_h_dtype] = subtract_fi('sub_dh', d_scale, h, 'latency', add_latency, 'full_precision', 1, 'type_a', d_dtype, 'type_b', h_dtype);
a_del = delay_srl_fi('a_del', a_scale, add_latency);

% adder layer 2
[a_plus_e_minus_f, a_plus_e_minus_f_dtype] = add_fi('add_aef', e_minus_f, a_del, 'latency', add_latency, 'full_precision', 1, 'type_a', e_minus_f_dtype, 'type_b', a_dtype);
[d_plus_c_minus_b, d_plus_c_minus_b_dtype] = add_fi('add_dcb', d_minus_h, c_minus_b, 'latency', add_latency, 'full_precision', 1, 'type_a', d_minus_h_dtype, 'type_b', c_minus_b_dtype);

[num_unr, num_unr_dtype] = add_fi('add', a_plus_e_minus_f, d_plus_c_minus_b, 'latency', add_latency', 'full_precision', 1, 'type_a', a_plus_e_minus_f_dtype, 'type_b', d_plus_c_minus_b_dtype);
num_sig = round_to_wordlength_fi('conv', num_unr, bit_width, num_unr_dtype);


num = num_sig;

end
