function [ab, ab_dtype] = subtract_fi(name, a, b, varargin)
% subtract(name, a, b, varargin)
ab = a - b;
ab_dtype = extract_fi_dtype(ab);

