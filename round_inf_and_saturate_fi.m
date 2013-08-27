function [x_rounded] = round_inf_and_saturate_fi(name, x, dtype, varargin)
defaults = {'logging', 1};
logging = get_var('logging', 'defaults', defaults, varargin{:});

quantizer = fixed.Quantizer('WordLength', dtype.WordLength, ...
    'FractionLength', dtype.FractionLength, ...
    'RoundingMethod', 'round', 'OverflowAction', 'saturate');

x_rounded = quantize(quantizer, x);
if logging
    orig_bitwidth = n_bits(x);
    new_bitwidth = n_bits(x_rounded);
    log_rounding(name, orig_bitwidth, new_bitwidth, double(x_rounded) - double(x));
end

end
