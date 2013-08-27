function [x_rounded] = round_inf_and_saturate_fi(name, x, dtype)

quantizer = fixed.Quantizer('WordLength', dtype.WordLength, ...
    'FractionLength', dtype.FractionLength, ...
    'RoundingMethod', 'round', 'OverflowAction', 'saturate');

x_rounded = quantize(quantizer, x);
log_rounding(name, double(x_rounded) - double(x));

end
