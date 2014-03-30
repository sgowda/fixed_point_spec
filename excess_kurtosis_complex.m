function [kappa_x, kurtosis_num, kurtosis_den] = excess_kurtosis_complex(x)

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
%
%--------------------------------------------------------------------------
%
%   Author: Andre Tkacenko
%
%   Revision History:
%       1.0 - 08/06/2013 - Original script started
%             08/06/2013 - Original script completed
%
%--------------------------------------------------------------------------


% Centralize the data vector about its mean
z = x - mean(x);

% Compute the mean of the square of z
E_z_2 = mean(z.^2);
E_z_2_abs_2 = (abs(E_z_2)).^2;

% Compute the mean of the magnitude square of z
E_abs_z_2 = mean((abs(z)).^2);
E_abs_z_2_2 = E_abs_z_2.^2;

% Compute the mean of the magnitude fourth power of z
E_abs_z_4 = mean((abs(z)).^4);

% Calculate the kurtosis using the above moments of z
kurtosis_num = (E_abs_z_4 - E_z_2_abs_2);
kurtosis_den = E_abs_z_2_2;
kappa_x = (E_abs_z_4 - E_z_2_abs_2 - (2.*E_abs_z_2_2))./(E_abs_z_2_2);