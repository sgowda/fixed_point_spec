function [single_type_x_mean, double_type_x_mean, triple_type_x_mean, ...
    quad_type_x_mean] = kurtosis_mean_types(type_x, acc_len)

single_type_x = type_x;
double_type_x = type_x^2 + type_x^2;
triple_type_x = type_x^3;
quad_type_x = double_type_x^2;

single_type_x_mean = fi_dtype(1, single_type_x.WordLength+acc_len, single_type_x.FractionLength+acc_len);
double_type_x_mean = fi_dtype(1, double_type_x.WordLength+acc_len, double_type_x.FractionLength+acc_len);
triple_type_x_mean = fi_dtype(1, triple_type_x.WordLength+acc_len, triple_type_x.FractionLength+acc_len);
quad_type_x_mean   = fi_dtype(1, quad_type_x.WordLength+acc_len, quad_type_x.FractionLength+acc_len);
