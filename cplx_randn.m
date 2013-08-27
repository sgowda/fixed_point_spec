function [x] = cplx_randn(size, stdev)

x_real = randn(size) * stdev;
x_imag = randn(size) * stdev;

x = complex(x_real, x_imag);

end