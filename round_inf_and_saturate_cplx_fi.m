function [x_rounded] = round_inf_and_saturate_cplx_fi(name, x, dtype)

x_real = real(x);
x_imag = imag(x);
x_real_rounded = round_inf_and_saturate_fi(name, x_real, dtype);
x_imag_rounded = round_inf_and_saturate_fi(name, x_imag, dtype);
x_rounded = complex(x_real_rounded, x_imag_rounded);

end
