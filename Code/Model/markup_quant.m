%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: markup_quant.m
% Author: Craig A. Chikis
% Date: 11/28/2021
% Note(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p50_approx, markup_array, gmarkup_all, markupL, markupF] =  markup_quant(kappa, lambda, s_max, firm_distributions, nu_s, quant)
    	
	if (kappa >= 9999)
		markup_array = [lambda.^(0:s_max); abs(firm_distributions)];
		[~, sort_idx] = sort(markup_array(1,:));

		markup_array = markup_array(:, sort_idx);
		markup_array(3,:) = cumsum(markup_array(2,:));

		markupF = ones(1,s_max+1);
		markupL = lambda.^(0:s_max); 
		
		gmarkup_all = [ones(1, s_max), lambda.^(0:s_max)];
	else
		weightL = nu_s.^(1-kappa)./(nu_s.^(1-kappa)+1);
		weightF = 1./(nu_s.^(1-kappa) + 1);

		markupL = ((kappa+nu_s.^(1-kappa))./(kappa-1));
		markupF = ((kappa*nu_s.^(1-kappa)+1)./((kappa-1)*nu_s.^(1-kappa)));

		gmarkup_all = [flip(markupF), markupL(2:end)];

		markup_array = [weightL.*markupL + weightF.*markupF; abs(firm_distributions)];
		[~, sort_idx] = sort(markup_array(1,:));

		markup_array = markup_array(:, sort_idx);
		markup_array(3,:) = cumsum(markup_array(2,:));

	end

	cutpoints = [quant] - markup_array(3, :);

	p50_loc = [find(cutpoints(1, :) > 0 ,1,'last'), find(cutpoints(1, :) < 0, 1,'first')];
	if (length(p50_loc) < 2)
		p50_loc(2) = p50_loc(1);
	end



	p50_approx = ((markup_array(1, p50_loc(2)) - markup_array(1, p50_loc(1)))/...
					(markup_array(3, p50_loc(2)) - markup_array(3, p50_loc(1))))*(quant - markup_array(3,p50_loc(1))) + ...
				markup_array(1, p50_loc(1));

	
	if (isnan(p50_approx))
		p50_approx = markup_array(1, p50_loc(1));
	end


end
