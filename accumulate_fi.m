function [x_acc] = accumulate_fi(x, dtype)
% Software implemenation of a fixed-point accumulator 
% where the bits do NOT grow as more numbers are added,
% just as in a real hardware implementation

quantizer = fixed.Quantizer('WordLength', dtype.WordLength, ...
    'FractionLength', dtype.FractionLength, ...
    'RoundingMethod', 'round', 'OverflowAction', 'saturate');

x_acc = fi(0, dtype.Signed, dtype.WordLength, dtype.FractionLength);

for k=1:length(x)
    x_acc_unr = x_acc + x(k);
    x_acc = quantize(quantizer, x_acc_unr);
end
