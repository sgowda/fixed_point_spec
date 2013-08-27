function [x_re_fl, x_im_fl, x_sq_re, x_sq_im, abs_x_sq, abs_x_4th, x_3rd_re, x_3rd_im] = kurtosis_origin_moments_float(x_re_fl, x_im_fl, varargin)

%-- Calculate origin_moments: 5 reg multipliers, 1 large multiplier, and 2 adders
x_re_sq = x_re_fl .* x_re_fl;
x_im_sq = x_im_fl .* x_im_fl;

abs_x_sq = x_re_sq + x_im_sq;
abs_x_4th = abs_x_sq .* abs_x_sq;

x_sq_re = x_re_sq - x_im_sq;

x_sq_im_unscaled = x_re_fl .* x_im_fl;
x_sq_im          = 2 * x_sq_im_unscaled;

x_3rd_re = x_re_fl .* abs_x_sq;
x_3rd_im = x_im_fl .* abs_x_sq;

end
