%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: logtfprdispersion.m
% Author: Craig A. Chikis
% Date: 09/03/2023
% Note(s): 
% MA-MFA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = logtfprdispersion(m)

    [p50_approx, markup_array, gmarkup_all, markupL, markupF] =  markup_quant(m.kappa, m.lambda, m.s_max, m.firm_distributions, m.nu_s, 0.5); 	
 

    vs = sqrt( ((log(markupL) - (log(markupL) + log(markupF))./2).^2  + ...
                (log(markupF) - (log(markupL) + log(markupF))./2).^2) ./ 2 );
    
    m.logtfprdispersion = sum(m.firm_distributions .* vs); 



end