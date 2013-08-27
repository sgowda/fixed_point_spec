function [x_rounded, x_rounded_dtype] = round_to_wordlength_fi(name, x, n_bits, type_x, varargin)
% Rounding to infinity with saturation
% [x_rounded, x_rounded_dtype] = round_to_wordlength_fi(name, x, n_bits, varargin)

n_int_bits = type_x.WordLength - type_x.FractionLength;
bin_pt = n_bits - n_int_bits;
quantizer = fixed.Quantizer('WordLength', n_bits, ...
    'FractionLength', bin_pt, 'RoundingMethod', 'round', 'OverflowAction', 'saturate');
x_rounded = quantize(quantizer, x);
x_rounded_dtype = extract_fi_dtype(x_rounded);
