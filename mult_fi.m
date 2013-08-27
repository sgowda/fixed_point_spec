function [ab, ab_dtype] = mult_fi(name, a, b, varargin)
% [ab] = mult_fi(name, a, b, varargin)

ab = a .* b;
ab_dtype = extract_fi_dtype(ab);