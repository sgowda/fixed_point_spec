function [kappa_x, kurtosis_num, kurtosis_den] = excess_kurtosis_complex_stream(x)

% excess_kurtosis_complex
%
%   kappa_x = excess_kurtosis_complex(x)
%
%   Calculates the complex random variable version of the excess kurtosis
%   of a set of sample data. The set of sample data can be either real or
%   complex. Here, the biased estimate of the kurtosis is calculated.
%
%   Inputs:
%       x           =   vector of input sample data (real or complex)
%
%   Outputs:
%       kappa_x     =   biased sample excess kurtosis of the vector x
%

% fourth abs moment
acc_len = log2(length(x));

m_x = mean(x);
a = sum(abs(x).^4);
b = 4*real( sum(x .* (abs(x).^2)) .* conj(m_x) );
c = 2*real( sum(x.^2) * conj(m_x^2) );
d = 4*sum( abs(x).^2 * abs(m_x)^2 );
e = abs(m_x)^4;
f = 4*e;
fourth_abs_moment = 2^-acc_len*(a - b + c + d) + (e - f);

% second abs moment
kurtosis_den = (mean(abs(x).^2) - abs(m_x).^2).^2;

% third term
third_term = abs( mean(x.^2) - m_x.^2 ).^2;
kurtosis_num = (fourth_abs_moment - third_term);
kappa_x = kurtosis_num/kurtosis_den - 2;
