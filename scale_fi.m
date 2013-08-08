function [x_scale] = scale_fi(name, x_fi, scale_factor)

% Fixed-point scaling by power of 2 (shift binary point)
new_fraction_length = x_fi.FractionLength - scale_factor;
x_scale = fi(double(x_fi) * 2^scale_factor, x_fi.Signed, ...
    x_fi.WordLength, new_fraction_length);
