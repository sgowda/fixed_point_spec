function [ab_re, ab_im, ab_dtype] = cplx_sub_fi(name, a, b, varargin)

if iscell(a) && iscell(b)
    [ab_re, ab_dtype] = subtract_fi([name, '_re'], a{1}, b{1}, varargin{:});    
    ab_im = subtract_fi([name, '_im'], a{2}, b{2}, varargin{:});    
else
    % TODO make cell array out of an array for signals
    NotImplementedError();
end
