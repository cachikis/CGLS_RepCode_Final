%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: sales_func.m
% Author: Craig A. Chikis
% Date: 07/18/2022 
% Note(s):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sales_divisor = sales_func(kappa, s_max, nu_s)
	if (kappa >= 9999)
		sales_divisor = [zeros(1,s_max), .5, repelem(1, length(1:s_max))];
	else
		del_s = nu_s.^(1-kappa)./(nu_s.^(1-kappa) + 1);
		del_ms = 1./(nu_s.^(1-kappa)+1);

		sales_divisor = [flip(del_ms), del_s(2:end)];
	end

end