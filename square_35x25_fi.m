function [x_sq] = square_35x25_fi(name, x, n_int_bits)

x_35bit = round_inf_and_saturate_fi('', x, fi_dtype(1, 35, 35-n_int_bits));
x_25bit = round_inf_and_saturate_fi('', x, fi_dtype(1, 25, 25-n_int_bits));
x_sq = x_35bit * x_25bit;

end