
tic;
mdl_name = 'butterfly';
T_sim = 2000

a = (2*rand(1, T_sim) - 1) + 1j*(2*rand(1, T_sim) - 1);
b = (2*rand(1, T_sim) - 1) + 1j*(2*rand(1, T_sim) - 1);
w = (2*rand(1, T_sim) - 1) + 1j*(2*rand(1, T_sim) - 1);

a_re  = timeseries(fi(real(a), 1, 18, 17));
a_im  = timeseries(fi(imag(a), 1, 18, 17));
b_re  = timeseries(fi(real(b), 1, 18, 17));
b_im  = timeseries(fi(imag(b), 1, 18, 17));
w_re  = timeseries(fi(real(w), 1, 18, 17));
w_im  = timeseries(fi(imag(w), 1, 18, 17));

set_param(mdl_name, 'StopTime', num2str(T_sim-1));
sim(mdl_name)
toc

apbw_fl = (a + b .* w).';
ambw_fl = (a - b .* w).';

apbw_fi = apbw_re + 1j*apbw_im;
ambw_fi = ambw_re + 1j*ambw_im;

apbw_error = apbw_fl(1:end-1) - apbw_fi(2:end);
ambw_error = ambw_fl(1:end-1) - ambw_fi(2:end);
max(abs(apbw_error(~logical(oflow(2:end)))))
max(abs(ambw_error(~logical(oflow(2:end)))))
