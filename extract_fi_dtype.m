function [dtype] = extract_fi_dtype(x)
% [dtype] = extract_fi_dtype(x)
dtype = fi_dtype(x.Signed, x.WordLength, x.FractionLength);
