%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: annualize_sim_data.m
% Author: Craig A. Chikis
% Date: 07/02/2020
% Note(s): 
% MA-MFA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function annual = annualize_sim_data(monthly, month_increment, type)

	

	if (strcmp(type, "sum"))

		it = reshape(monthly, size(monthly,1), month_increment, size(monthly,2)/month_increment);
		it = nansum(it, 2);

		annual = reshape(it, size(it,1), size(it,3));


	elseif (strcmp(type, "EoP"))

		it = reshape(monthly, size(monthly,1), month_increment, size(monthly,2)/month_increment);
		it = it(:, size(it,2), :);
		annual = reshape(it, size(it,1), size(it,3));
	elseif (strcmp(type, "BoP"))

		it = reshape(monthly, size(monthly,1), month_increment, size(monthly,2)/month_increment);
		it = it(:, 1, :);
		annual = reshape(it, size(it,1), size(it,3));

	end









end
