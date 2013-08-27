function [x_acc, x_acc_dtype] = accumulate_fi(x, dtype)

x_acc = sum(x);
x_acc_dtype = extract_fi_dtype(x_acc);

