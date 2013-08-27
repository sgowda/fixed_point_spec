function [x_rounded] = round_inf_and_saturate_cplx_fi(name, x, dtype, varargin)
defaults = {'logging', 0};
logging = get_var('logging', 'defaults', defaults, varargin{:});
x_real = real(x);
x_imag = imag(x);
x_real_rounded = round_inf_and_saturate_fi([name, '_re'], x_real, dtype, 'logging', logging);
x_imag_rounded = round_inf_and_saturate_fi([name, '_im'], x_imag, dtype, 'logging', logging);
x_rounded = complex(x_real_rounded, x_imag_rounded);

end
