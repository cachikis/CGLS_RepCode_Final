%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: solve_rho_impl.m
% Author: Craig A. Chikis
% Note(s):
% MA-MFA
% Date: 10/10/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = solve_rho_impl(rho, lambda, s, kappa)

	

	out = zeros(max(size(s)), 1);


	out = rho - (lambda.^(-s)).*(kappa + rho.^(1 - kappa))./(kappa + rho.^(kappa - 1));


	out = reshape(out, max(size(out)),1);

end
