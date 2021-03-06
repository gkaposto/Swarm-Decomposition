function [x ,x_c] = FMGenerator(Q, A, F, beta, theta, Fs)

% --------------------------------------------------------------------------------------------------
%   FMGenerator: Frequency modulated trigonometric signal generator
%   args: - Q: signal's length
%         - A: vector that contains the amplitude of every component, the length of the vector indicates the number of generated components
%         - F: vector that contains the frequency of every component,          -"-
%         - beta: vector contains the frequency modulation factor of every component,               -"-
%         - theta: vector contains the phase of every component,               -"-
%         - Fs: sampling frequency
%   returns:
%         - x:   the generated signal
%         - x_c: the individual components
%   
%   This signal model is adopted from the paper:
%   "An iterative approach for decomposition of multi-component non-stationary signals based on eigenvalue decomposition of the Hankel matrix",
%   P.Jain & R.B.Pachori, Journal of the Franklin Institute (2015) 
%
%   Developer: Georgios Apostolidis
% --------------------------------------------------------------------------------------------------

L_A  = length(A);
L_F  = length(F);
L_b  = length(beta);
L_th = length(theta);
 
if (L_A == L_F) && (L_A == L_b) && (L_A == L_th)
    L = L_A;
else
    error('error with parameters lengths');
end

if nargin == 5
    Fs = 1;
end

x_c = zeros(L, Q);
n   = 0:1:Q-1;
f   = F / Fs;

for i = 1:1:L
    x_c(i, :) = A(i) * cos(2 * pi * f(i) * (1 + beta(i)) * n + theta(i));  
end

x = sum(x_c, 1);

end

