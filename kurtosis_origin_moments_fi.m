function [x_re_fi, x_im_fi, x_sq_re, x_sq_im, abs_x_sq, abs_x_4th, x_3rd_re, x_3rd_im, origin_moment_error] = kurtosis_origin_moments_fi(x_re_fi, x_im_fi, varargin)
type_x_default = fi_dtype(1, 18, 17);
defaults = {'type_x', type_x_default};
type_x = get_var('type_x', 'defaults', defaults, varargin{:});

%-- Calculate origin_moments: 5 reg multipliers, 1 large multiplier, and 2 adders
[x_re_sq, x_re_sq_dtype] = mult_fi('mult1', x_re_fi, x_re_fi, 'type_a', type_x, 'type_b', type_x);
[x_im_sq, x_im_sq_dtype] = mult_fi('mult2', x_im_fi, x_im_fi, 'type_a', type_x, 'type_b', type_x);

[abs_x_sq, abs_x_sq_dtype] = add_fi('add1', x_re_sq, x_im_sq, 'type_a', x_re_sq_dtype, 'type_b', x_im_sq_dtype);
[abs_x_4th, abs_x_4th_dtype] = mult_fi('mult3', abs_x_sq, abs_x_sq, 'type_a', abs_x_sq_dtype, 'type_b', abs_x_sq_dtype);

[x_sq_re, x_sq_re_dtype] = subtract_fi('add2', x_re_sq, x_im_sq, 'type_a', x_re_sq_dtype, 'type_b', x_im_sq_dtype);
[x_sq_im_unscaled, x_sq_im_unscaled_dtype] = mult_fi('mult4', x_re_fi, x_im_fi, 'type_a', type_x, 'type_b', type_x);
[x_sq_im, x_sq_im_dtype] = scale_fi('', x_sq_im_unscaled, 1, 'type_x', x_sq_im_unscaled_dtype);

[x_3rd_re, x_3rd_re_dtype] = mult_fi('mult4', x_re_fi, abs_x_sq, 'type_a', type_x, 'type_b', abs_x_sq_dtype);
[x_3rd_im, x_3rd_im_dtype] = mult_fi('mult5', x_im_fi, abs_x_sq, 'type_a', type_x, 'type_b', abs_x_sq_dtype);

x_re_fl = double(x_re_fi);
x_im_fl = double(x_im_fi);
[x_re_fl, x_im_fl, x_sq_re_fl, x_sq_im_fl, abs_x_sq_fl, abs_x_4th_fl, ...
    x_3rd_re_fl, x_3rd_im_fl] = kurtosis_origin_moments_float(x_re_fl, x_im_fl);

origin_moment_error.x_sq_re = max(max(abs(x_sq_re_fl - double(x_sq_re))));
origin_moment_error.x_sq_im = max(max(abs(x_sq_im_fl - double(x_sq_im))));
origin_moment_error.x_3rd_re = max(max(abs(x_3rd_re_fl - double(x_3rd_re))));
origin_moment_error.x_3rd_im = max(max(abs(x_3rd_im_fl - double(x_3rd_im))));
origin_moment_error.abs_x_sq = max(max(abs(abs_x_sq_fl - double(abs_x_sq))));
origin_moment_error.abs_x_4th = max(max(abs(abs_x_4th_fl - double(abs_x_4th))));

end
