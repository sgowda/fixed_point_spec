function [apbw_re, apbw_im, ambw_re, ambw_im, oflow] = ...
    butterfly_fi(a_re, a_im, b_re, b_im, w_re, w_im, shift)
% radix-2 butterfly
% [apbw_re, apbw_im, ambw_re, ambw_im] = ...
%     butterfly_fi(a_re, a_im, b_re, b_im, w_re, w_im)

if ~isa(a_re, 'embedded.fi') || ~isa(a_im, 'embedded.fi') || ~isa(b_re, 'embedded.fi') || ~isa(b_im, 'embedded.fi') || ~isa(w_re, 'embedded.fi') || ~isa(w_im, 'embedded.fi')
    error('Inputs must be fixed point!')    
end

bw_re_unrounded = b_re*w_re - b_im*w_im;
bw_im_unrounded = b_re*w_im + b_im*w_re;

cmult_quantizer = fixed.Quantizer('WordLength', 22, 'FractionLength', 19, ...
    'RoundingMethod', 'round', 'OverflowAction', 'saturate');
bw_re_rounded = quantize(cmult_quantizer, bw_re_unrounded);
bw_im_rounded = quantize(cmult_quantizer, bw_im_unrounded);

apbw_re_unrounded = bw_re_rounded + a_re;
apbw_im_unrounded = bw_im_rounded + a_im;

ambw_re_unrounded = a_re - bw_re_rounded;
ambw_im_unrounded = a_im - bw_im_rounded;

oflow = 0;
if (apbw_re_unrounded >= 1) || (apbw_re_unrounded < -1)
    oflow = 1;
elseif (apbw_im_unrounded >= 1) || (apbw_im_unrounded < -1)
    oflow = 1;
elseif (ambw_re_unrounded >= 1) || (ambw_re_unrounded < -1)
    oflow = 1;
elseif (ambw_im_unrounded >= 1) || (ambw_im_unrounded < -1)
    oflow = 1;
end
    
out_quantizer = fixed.Quantizer('WordLength', 18, 'FractionLength', 17, ...
    'RoundingMethod', 'round', 'OverflowAction', 'saturate');
apbw_re = quantize(out_quantizer, apbw_re_unrounded);
apbw_im = quantize(out_quantizer, apbw_im_unrounded);
ambw_re = quantize(out_quantizer, ambw_re_unrounded);
ambw_im = quantize(out_quantizer, ambw_im_unrounded);


end
