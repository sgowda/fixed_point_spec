function [num, den, abs_x_sq_mean] = kurtosis_moment_calc_fi(x_mean_re, x_mean_im, x_sq_mean_re, x_sq_mean_im, abs_x_sq_mean, abs_x_4th_mean, x_3rd_mean_re, x_3rd_mean_im, varargin)
defaults = {'type_x', fi_dtype(1, 18, 17), 'acc_len', 14};
type_x = get_var('type_x', 'defaults', defaults, varargin{:});
acc_len = get_var('acc_len', 'defaults', defaults, varargin{:});
[m_x_type, x_sq_type, x_3rd_type, x_4th_type] = kurtosis_mean_types(type_x, acc_len);

bit_width_power = 32;
bit_width_num = 96;
bit_width_den = 64;

rounding_latency = 1;
add_latency_96bit = 4;
cross_product_latency = 15;
adder_tree_latency = add_latency_96bit * 3 + rounding_latency;
abs_x_sq_mean_del = delay_srl_fi('cross_product_del', abs_x_sq_mean, cross_product_latency);

abs_x_sq_mean = round_to_wordlength_fi('adder_tree_del', abs_x_sq_mean_del, bit_width_power, x_sq_type, 'latency', adder_tree_latency);

% cross-products
[d, f, a, c, b, abs_mx_sq, e, h] = kurtosis_cross_products_fi(x_mean_re, x_mean_im, x_sq_mean_re, x_sq_mean_im, abs_x_sq_mean, abs_x_4th_mean, x_3rd_mean_re, x_3rd_mean_im, 'total_latency', cross_product_latency, 'type_x', type_x, 'acc_len', acc_len);

% cross-product data types
[a_dtype, b_dtype, c_dtype, d_dtype, e_dtype, f_dtype, abs_m_x_sq_dtype, abs_mean_x_sq_dtype] = ...
    kurtosis_cross_products_fi(m_x_type, m_x_type, ...
    x_sq_type, x_sq_type, x_sq_type, x_4th_type, x_3rd_type, x_3rd_type);

% denominator
den = kurtosis_den(abs_x_sq_mean_del, abs_mx_sq, 'bit_width', bit_width_den, 'abs_m_x_sq_dtype', abs_m_x_sq_dtype, 'abs_x_sq_mean_del_dtype', x_sq_type);

% numerator
num = kurtosis_num(d, f, a, c, b, e, h, 'bit_width', bit_width_num, 'a_dtype', a_dtype, ...
        'b_dtype', b_dtype, 'c_dtype', c_dtype, 'd_dtype', d_dtype, ...
        'e_dtype', e_dtype, 'f_dtype', f_dtype, 'h_dtype', abs_mean_x_sq_dtype);
    
end

function [den] = kurtosis_den(mean_abs_x_sq, abs_m_x_sq, varargin)
defaults = {'bit_width', 64, 'abs_m_x_sq_dtype'};
bit_width = get_var('bit_width', 'defaults', defaults, varargin{:});
abs_m_x_sq_dtype = get_var('abs_m_x_sq_dtype', 'defaults', defaults, varargin{:});
mean_abs_x_sq_dtype = get_var('abs_x_sq_mean_del_dtype', 'defaults', defaults, varargin{:});

[sec_moment_unr, sec_moment_unr_dtype] = subtract_fi('Sub', mean_abs_x_sq, abs_m_x_sq, 'latency', 4, 'type_a', mean_abs_x_sq_dtype, 'type_b', abs_m_x_sq_dtype);
[den_unr, den_unr_dtype] = mult_35x25_fi('square', sec_moment_unr, sec_moment_unr, 'type_a', sec_moment_unr_dtype, 'type_b', sec_moment_unr_dtype);
[den_sig, den_dtype] = round_to_wordlength_fi('conv', den_unr, bit_width, den_unr_dtype);
den = den_sig;

end

function num = kurtosis_num(d, f, a, c, b, e, h, varargin)
defaults = {'bit_width', 96, 'conv_latency', 1};
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

% adder layer 1
[e_minus_f, e_minus_f_dtype] = subtract_fi('sub_ef', e, f, 'latency', add_latency, 'full_precision', 1, 'type_a', e_dtype, 'type_b', f_dtype);
[c_minus_b, c_minus_b_dtype] = subtract_fi('sub_cb', c, b, 'latency', add_latency, 'full_precision', 1, 'type_a', c_dtype, 'type_b', b_dtype);
[d_minus_h, d_minus_h_dtype] = subtract_fi('sub_dh', d, h, 'latency', add_latency, 'full_precision', 1, 'type_a', d_dtype, 'type_b', h_dtype);
a_del = delay_srl_fi('a_del', a, add_latency);

% adder layer 2
[a_plus_e_minus_f, a_plus_e_minus_f_dtype] = add_fi('add_aef', e_minus_f, a_del, 'latency', add_latency, 'full_precision', 1, 'type_a', e_minus_f_dtype, 'type_b', a_dtype);
[d_plus_c_minus_b, d_plus_c_minus_b_dtype] = add_fi('add_dcb', d_minus_h, c_minus_b, 'latency', add_latency, 'full_precision', 1, 'type_a', d_minus_h_dtype, 'type_b', c_minus_b_dtype);

[num_unr, num_unr_dtype] = add_fi('add', a_plus_e_minus_f, d_plus_c_minus_b, 'latency', add_latency', 'full_precision', 1, 'type_a', a_plus_e_minus_f_dtype, 'type_b', d_plus_c_minus_b_dtype);
num_sig = round_to_wordlength_fi('conv', num_unr, bit_width, num_unr_dtype);


num = num_sig;

end
