function [ab, ab_dtype] = add_fi(name, a, b, varargin)
% add_fi(name, a, b, varargin)
ab = a + b;
ab_dtype = extract_fi_dtype(ab);

