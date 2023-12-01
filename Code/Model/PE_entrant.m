%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: PE_entrant.m
% Author: Craig A. Chikis
% Date: 06/24/2020
% Note(s):
% MA-MFA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = PE_entrant(m)

	probE = m.probE;
	investment_entrant = m.investment_entrant;
	subsidyE = m.subsidyE; 
	B_entrant = m.B_entrant; 
	gamma = m.gamma; 
	s_max = m.s_max; 
	value_functions = m.value_functions;
	firm_distributions = m.firm_distributions; 
	parv_parrho = m.par_Dv_parrho;
	parv_parw = m.par_Dv_parw; 									 
	wage_share = m.wage_share; 
	entry_type = m.entry_type;

	
	probE = flip(probE, 1);

	G = @(x, B) (x/B).^(1/gamma);
	GP = @(x, B) (1/(B*gamma))*(x/B).^((1/gamma)-1);
	GdP = @(x, B) (1/((B^2)*gamma))*((1/gamma)-1)*(x/B).^((1/gamma)-2);
	cs = @(x, B, subsidy) 1./(GdP(x, B)*wage_share.*(1-subsidy));

	A_rho = zeros(s_max+1);
	b_rho = zeros(s_max+1, 1);

	A_w = zeros(s_max+1);
	b_w = zeros(s_max+1, 1);

	if (strcmp(entry_type, "Undirected"))

		other_termw = zeros(1, s_max+1);
		other_termrho = zeros(1, s_max+1);

		for (sigma = 1:length(other_termw))
			other_termw(sigma) = sum(probE(sigma, (s_max+3-sigma):end).*parv_parw((s_max+3-sigma):end));
			other_termrho(sigma) = sum(probE(sigma, (s_max+3-sigma):end).*parv_parrho((s_max+3-sigma):end));
		end

		parve_parw = ...
			 -(1-subsidyE(1))*G(investment_entrant(1), B_entrant) + investment_entrant(1)*sum(firm_distributions.*other_termw);
		parve_parrho = investment_entrant(1)*sum(firm_distributions.*other_termrho);

		parxe_parrho = cs(investment_entrant(1), B_entrant, subsidyE(1))*sum(firm_distributions.*other_termrho);
		parxe_parw = cs(investment_entrant(1), B_entrant, subsidyE(1))*(-GP(investment_entrant(1), B_entrant)*(1-subsidyE(1)) + ...
																		sum(firm_distributions.*other_termw));

		parDve_parrho = sum(firm_distributions.*other_termrho);
		parDve_parw = sum(firm_distributions.*other_termw);


	else

		parve_parw = 0;
		parve_parrho = 0;
		parxe_parrho = 0;
		parxe_parw = 0;
		parDve_parrho = 0;
		parDve_parw = 0;


	end


	m.parve_parw = parve_parw; 
	m.parve_parrho = parve_parrho;
	m.parxe_parrho = parxe_parrho; 
	m.parxe_parw = parxe_parw; 
	m.parDve_parrho = parDve_parrho; 
	m.parDve_parw = parDve_parw;

	

end
