%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: sim_noentry.m
% Author: Craig A. Chikis
% Date: 10/18/2022
% Note(s):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = sim_noentry(m)


	N = m.N;
	T = m.T;
	investment = m.investment;  
	fric_binding = m.fric_binding; 
	investment_entrant = m.investment_entrant; 
	eta = m.eta;
	prob = m.prob; 
	prob_exog = m.prob_exog; 
	probE = m.probE; 
	B = m.B;
	growth_rate = m.growth_rate; 
	firm_distributions = m.firm_distributions; 
	lab_demand = m.lab_demand; 
	value_functions = m.value_functions; 
	kappa = m.kappa; 
	nu_s = m.nu_s; 
	profit = m.profit; 
	subsidy = m.subsidy; 
	gamma = m.gamma; 
	entrance = m.entrance; 
	zeta = m.zeta; 
	lambda = m.lambda; 
	fcitcount = m.fcitcount; 
	numCPC = m.numCPC; 
	nbucket = m.nbucket; 
	nbucket_patent = m.nbucket_patent; 
	qtrtype = m.qtrtype; 
	winsor_vec_drop = m.winsor_vec_drop;  
	prctile_inp = m.prctile_inp; 
	wage_share = m.wage_share; 


    if (entrance)
        error('Wrong code for entrance.')
    end

	s_max = (size(prob, 1) - 1)/2;
	[~, tmp2] = find(prob(1:s_max, :));
	l = max(tmp2) - (s_max+1);
	
    rng(20220830, 'Threefry') 

    s_max = (length(investment) - 1)/2;

    if ((2*max(investment) + max(eta) + max(investment_entrant)) > 1)
		scale_time = 0.99/(2*max(investment) + max(eta) + max(investment_entrant));
	else
		scale_time = 1;
	end

    N = 100*ceil(N/100);

	if (scale_time < 0.1)
		error('scale_time is too small')
	end

    divisor = 4*round((12*1/scale_time)/4);

	T = divisor*ceil(T/scale_time/divisor);

	draws = rand(N, T);
	coinflip = rand(N, T);

    investment = scale_time*investment;
	eta = scale_time*eta;
	investment_entrant = scale_time*investment_entrant;
	B = B*scale_time;
	growth_rate = power(1+growth_rate, 12)-1; 
	growth_rate = power(1+growth_rate, 1/(12/scale_time))-1; 
    value_functions = value_functions/scale_time;

    ageL = Inf(N, T);
	ageF = Inf(N, T);

	leader_innovL = zeros(N, T);
	laggard_innovL = zeros(N,T);
	entrant_innovL = zeros(N,T);
	

	leader_innovF = zeros(N, T);
	laggard_innovF = zeros(N,T);
	entrant_innovF = zeros(N,T);

	leader_innovL2 = zeros(N, T);
	laggard_innovL2 = zeros(N,T);
	entrant_innovL2 = zeros(N,T);
	leader_innovF2 = zeros(N, T);
	laggard_innovF2 = zeros(N,T);
	entrant_innovF2 = zeros(N,T);
	tied_innovL = zeros(N, T);
	tied_innovF = zeros(N, T);
	entrant_tied_innovL = zeros(N, T);
	entrant_tied_innovF = zeros(N, T);

	stateL = zeros(N, T);
	stateF = zeros(N, T);

	prodL = ones(N, T);
	prodF = ones(N, T); 
 
	xiL = zeros(N, T);
	xiF = zeros(N, T);

    lab_demandL = zeros(N,T);
	lab_demandF = zeros(N,T);
    
	init_dist = firm_distributions*N;

	tmp = [];
	for (s = 1:length(init_dist))
		tmp = [tmp, repelem(s-1, round(init_dist(s)))];
	end

	if (length(tmp) < size(stateL,1))
		tmp = [tmp, repelem(find(init_dist == max(init_dist))-1, size(stateL,1)-length(tmp))];
	else
		tmp = tmp(1:size(stateL,1));
	end

	column_index = -s_max:s_max;


    stateL(:, 1) = tmp;
	stateF(:, 1) = -tmp;

    lab_demandL(:, 1) = lab_demand{1}(s_max+1+stateL(:,1)) ./ (1 - subsidy(s_max+1+stateL(:, 1))) + ...
                        lab_demand{2}(s_max+1+stateL(:,1));
    
    lab_demandF(:, 1) = lab_demand{1}(s_max+1+stateF(:,1)) ./ (1 - subsidy(s_max+1+stateF(:, 1))) + ...
                        lab_demand{2}(s_max+1+stateF(:,1));

    nothing_o = 1 - investment((s_max+1):end) - flip(investment(1:(s_max+1))) - eta - investment_entrant;


	L_outcome_o = (investment((s_max+1):end)').*prob((s_max+1):end, :);
	F_outcome_o = flip((investment(1:(s_max+1))').*prob(1:(s_max+1), :));
	E_outcome_o = (investment_entrant').*flip(probE, 1);
	exog_outcome_o = [zeros(1, 2*s_max+1); (eta(2:end)').*flip(prob_exog, 1)];


    for (t = 1:(T-1))
		
		for (n = 1:N)
	
			if (draws(n, t) > nothing_o(abs(stateL(n,t))+1))

				if (stateL(n, t) > 0)

					start_pt = nothing_o(abs(stateL(n,t))+1);
					for (k = 1:(2*s_max+1))
						if (draws(n, t) >= start_pt) && (draws(n,t) < (start_pt + L_outcome_o(abs(stateL(n,t))+1, k)))

							stateL(n, t+1) = column_index(k);
							stateF(n,t+1) = -column_index(k);

							prodL(n, t+1) = (m.lambda^(stateL(n,t+1)-stateL(n,t))) * prodL(n,t);  
							prodF(n, t+1) = prodF(n,t); 


							ageL(n,t+1) = ageL(n,t)+1;
							ageF(n,t+1) = ageF(n,t)+1;

							xiL(n, t+1) = value_functions(find(column_index == stateL(n, t+1))) - ...
										  value_functions(find(column_index == stateL(n, t)));

							leader_innovL(n,t+1) = 1;
							leader_innovL2(n,t+1) = 1;



						elseif (draws(n,t) >= (start_pt + L_outcome_o(abs(stateL(n,t))+1, k))) && (draws(n,t) < (start_pt+L_outcome_o(abs(stateL(n,t))+1, k)+F_outcome_o(abs(stateL(n,t))+1, k)))
							stateF(n,t+1) = column_index(k);
							stateL(n,t+1) = -column_index(k);

							prodF(n, t+1) = (m.lambda^(stateF(n,t+1)-stateF(n,t))) * prodF(n,t); 
							prodL(n,t+1) = prodL(n,t); 

							ageL(n,t+1) = ageL(n,t)+1;
							ageF(n,t+1) = ageF(n,t)+1;

							xiF(n, t+1) = value_functions(find(column_index == stateF(n, t+1))) - ...
										  value_functions(find(column_index == stateF(n, t)));

							laggard_innovF(n,t+1) = 1;
							laggard_innovF2(n,t+1) = 1;


						elseif (draws(n,t) >= (start_pt+L_outcome_o(abs(stateL(n,t))+1, k)+F_outcome_o(abs(stateL(n,t))+1, k))) && ...
								(draws(n,t) < (start_pt+L_outcome_o(abs(stateL(n,t))+1, k)+F_outcome_o(abs(stateL(n,t))+1, k)+E_outcome_o(abs(stateL(n,t))+1, k)))

							stateF(n,t+1) = column_index(k);
							stateL(n,t+1) = -column_index(k);

							prodF(n, t+1) = (m.lambda^(stateF(n,t+1)-stateF(n,t))) * prodF(n,t); 
							prodL(n, t+1) = prodL(n, t); 


							ageF(n,t+1) = 0;
							ageL(n,t+1) = ageL(n,t)+1;


							entrant_innovF(n,t+1) = 1;
							entrant_innovF2(n,t+1) = 1;


							xiF(n,t+1) = nan;


						elseif (draws(n,t) >= (start_pt+L_outcome_o(abs(stateL(n,t))+1, k)+F_outcome_o(abs(stateL(n,t))+1, k)+E_outcome_o(abs(stateL(n,t))+1, k))) && ...
								(draws(n,t) < (start_pt+L_outcome_o(abs(stateL(n,t))+1, k)+F_outcome_o(abs(stateL(n,t))+1, k)+E_outcome_o(abs(stateL(n,t))+1, k)+exog_outcome_o(abs(stateL(n,t))+1, k)))

							stateF(n,t+1) = column_index(k);
							stateL(n,t+1) = -column_index(k);

							prodL(n,t+1) = prodL(n,t);
							prodF(n,t+1) = prodF(n,t); 


							ageL(n,t+1) = ageL(n,t)+1;
							ageF(n,t+1) = ageF(n,t)+1;

							laggard_innovF(n,t+1) = 1;
							

						end

						start_pt = start_pt + L_outcome_o(abs(stateL(n,t))+1, k) + F_outcome_o(abs(stateL(n,t))+1, k) + E_outcome_o(abs(stateL(n,t))+1, k) + exog_outcome_o(abs(stateL(n,t))+1, k);
					end

				elseif (stateL(n, t) < 0)
					start_pt = nothing_o(abs(stateL(n,t))+1);
					for (k = 1:(2*s_max+1))
						if (draws(n, t) >= start_pt) && (draws(n,t) < (start_pt + L_outcome_o(abs(stateL(n,t))+1, k)))

							stateF(n, t+1) = column_index(k);
							stateL(n,t+1) = -column_index(k);

							prodF(n, t+1) = (m.lambda^(stateF(n,t+1)-stateF(n,t))) * prodF(n,t); 
							prodL(n, t+1) = prodL(n, t); 

							ageL(n,t+1) = ageL(n,t)+1;
							ageF(n,t+1) = ageF(n,t)+1;

							xiF(n, t+1) = value_functions(find(column_index == stateF(n, t+1))) - ...
										  value_functions(find(column_index == stateF(n, t)));

							leader_innovF(n,t+1) = 1;
							leader_innovF2(n,t+1) = 1;


						elseif (draws(n,t) >= (start_pt + L_outcome_o(abs(stateL(n,t))+1, k))) && (draws(n,t) < (start_pt+L_outcome_o(abs(stateL(n,t))+1, k)+F_outcome_o(abs(stateL(n,t))+1, k)))
							stateL(n,t+1) = column_index(k);
							stateF(n,t+1) = -column_index(k);

							prodL(n, t+1) = (m.lambda^(stateL(n,t+1)-stateL(n,t))) * prodL(n,t); 
							prodF(n, t+1) = prodF(n, t);


							ageL(n,t+1) = ageL(n,t)+1;
							ageF(n,t+1) = ageF(n,t)+1;

							xiL(n, t+1) = value_functions(find(column_index == stateL(n,t+1))) - ...
										  value_functions(find(column_index == stateL(n,t)));

							laggard_innovL(n,t+1) = 1;
							laggard_innovL2(n,t+1) = 1;



						elseif (draws(n,t) >= (start_pt+L_outcome_o(abs(stateL(n,t))+1, k)+F_outcome_o(abs(stateL(n,t))+1, k))) && ...
								(draws(n,t) < (start_pt+L_outcome_o(abs(stateL(n,t))+1, k)+F_outcome_o(abs(stateL(n,t))+1, k)+E_outcome_o(abs(stateL(n,t))+1, k)))

							stateL(n,t+1) = column_index(k);
							stateF(n,t+1) = -column_index(k);

							prodL(n, t+1) = (m.lambda^(stateL(n,t+1)-stateL(n,t))) * prodL(n,t); 
							prodF(n, t+1) = prodF(n, t); 


							ageF(n,t+1) = ageF(n,t)+1;
							ageL(n,t+1) = 0;


							entrant_innovL(n,t+1) = 1;
							entrant_innovL2(n,t+1) = 1;


							xiL(n,t+1) = nan;

							

						elseif (draws(n,t) >= (start_pt+L_outcome_o(abs(stateL(n,t))+1, k)+F_outcome_o(abs(stateL(n,t))+1, k)+E_outcome_o(abs(stateL(n,t))+1, k))) && ...
								(draws(n,t) < (start_pt+L_outcome_o(abs(stateL(n,t))+1, k)+F_outcome_o(abs(stateL(n,t))+1, k)+E_outcome_o(abs(stateL(n,t))+1, k)+exog_outcome_o(abs(stateL(n,t))+1, k)))

							stateL(n,t+1) = column_index(k);
							stateF(n,t+1) = -column_index(k);

							prodL(n,t+1) = prodL(n, t);
							prodF(n,t+1) = prodF(n, t);

							ageL(n,t+1) = ageL(n,t)+1;
							ageF(n,t+1) = ageF(n,t)+1;

							laggard_innovL(n,t+1) = 1;



						end

						start_pt = start_pt + L_outcome_o(abs(stateL(n,t))+1, k) + F_outcome_o(abs(stateL(n,t))+1, k) + E_outcome_o(abs(stateL(n,t))+1, k) + exog_outcome_o(abs(stateL(n,t))+1, k);
					end



				elseif (stateL(n, t) == 0)

					coinflip_num = true*(coinflip(n,t) > .5) + false*(coinflip(n,t) <= .5);

					start_pt = nothing_o(abs(stateL(n,t))+1);

					for (k = 1:(2*s_max+1))
						if (draws(n, t) >= start_pt) && (draws(n,t) < (start_pt + L_outcome_o(abs(stateL(n,t))+1, k)))
						
							stateF(n, t+1) = column_index(k);
							stateL(n,t+1) = -column_index(k);

							prodF(n, t+1) = (m.lambda^(stateF(n,t+1)-stateF(n,t))) * prodF(n,t); 
							prodL(n, t+1) = prodL(n, t);

							ageL(n,t+1) = ageL(n,t)+1;
							ageF(n,t+1) = ageF(n,t)+1;

							xiF(n,t+1) = value_functions(find(column_index == stateF(n,t+1))) - ...
										 value_functions(find(column_index == stateF(n,t)));

							leader_innovF(n,t+1) = 1;
							tied_innovF(n, t+1) = 1;

						elseif (draws(n,t) >= (start_pt + L_outcome_o(abs(stateL(n,t))+1, k))) && (draws(n,t) < (start_pt+L_outcome_o(abs(stateL(n,t))+1, k)+F_outcome_o(abs(stateL(n,t))+1, k)))
							
							stateL(n,t+1) = column_index(k);
							stateF(n,t+1) = -column_index(k);

							prodL(n, t+1) = (m.lambda^(stateL(n,t+1)-stateL(n,t))) * prodL(n,t); 
							prodF(n, t+1) = prodF(n, t);


							ageL(n,t+1) = ageL(n,t)+1;
							ageF(n,t+1) = ageF(n,t)+1;

							xiL(n,t+1) = value_functions(find(column_index == stateL(n,t+1))) - ...
										 value_functions(find(column_index == stateL(n,t)));

							leader_innovL(n,t+1) = 1;
							tied_innovL(n, t+1) = 1;

						elseif (draws(n,t) >= (start_pt+L_outcome_o(abs(stateL(n,t))+1, k)+F_outcome_o(abs(stateL(n,t))+1, k))) && ...
								(draws(n,t) < (start_pt+L_outcome_o(abs(stateL(n,t))+1, k)+F_outcome_o(abs(stateL(n,t))+1, k)+E_outcome_o(abs(stateL(n,t))+1, k)))
							
							if (coinflip_num)
								stateL(n,t+1) = column_index(k);
								stateF(n,t+1) = -column_index(k);

								prodL(n, t+1) = (m.lambda^(stateL(n,t+1)-stateL(n,t))) * prodL(n,t); 
								prodF(n, t+1) = prodF(n, t); 


								ageF(n,t+1) = ageF(n,t)+1;
								ageL(n,t+1) = 0;


								entrant_innovL(n,t+1) = 1;
								entrant_tied_innovL(n, t+1) = 1;

								xiL(n,t+1) = nan;

							else
								stateF(n,t+1) = column_index(k);
								stateL(n,t+1) = -column_index(k);

								prodF(n, t+1) = (m.lambda^(stateF(n,t+1)-stateF(n,t))) * prodF(n,t);
								prodL(n, t+1) = prodL(n, t); 


								ageF(n,t+1) = 0;
								ageL(n,t+1) = ageL(n,t)+1;


								entrant_innovF(n,t+1) = 1;
								entrant_tied_innovF(n, t+1) = 1;

								xiF(n,t+1) = nan;
							end
						
						elseif (draws(n,t) >= (start_pt+L_outcome_o(abs(stateL(n,t))+1, k)+F_outcome_o(abs(stateL(n,t))+1, k)+E_outcome_o(abs(stateL(n,t))+1, k))) && ...
								(draws(n,t) < (start_pt+L_outcome_o(abs(stateL(n,t))+1, k)+F_outcome_o(abs(stateL(n,t))+1, k)+E_outcome_o(abs(stateL(n,t))+1, k)+exog_outcome_o(abs(stateL(n,t))+1, k)))

							stateL(n,t+1) = column_index(k);
							stateF(n,t+1) = -column_index(k);

							prodL(n, t+1) = prodL(n, t);
							prodF(n, t+1) = prodF(n, t); 

							ageL(n,t+1) = ageL(n,t)+1;
							ageF(n,t+1) = ageF(n,t)+1;



						end

						start_pt = start_pt + L_outcome_o(abs(stateL(n,t))+1, k) + F_outcome_o(abs(stateL(n,t))+1, k) + E_outcome_o(abs(stateL(n,t))+1, k) + exog_outcome_o(abs(stateL(n,t))+1, k);
					end

				end
			




			else
				
				stateL(n, t+1) = stateL(n, t);
				stateF(n, t+1) = stateF(n, t);
				ageL(n,t+1) = ageL(n,t) + 1;
				ageF(n,t+1) = ageF(n,t) + 1;
				prodL(n,t+1) = prodL(n,t);
				prodF(n,t+1) = prodF(n,t); 
			end




		end

		lab_demandL(:, t+1) = lab_demand{1}(s_max+1+stateL(:,t+1)) ./ (1 - subsidy(s_max+1+stateL(:, t+1))) + ...
							  lab_demand{2}(s_max+1+stateL(:,t+1));
		lab_demandF(:, t+1) = lab_demand{1}(s_max+1+stateF(:,t+1)) ./ (1 - subsidy(s_max+1+stateF(:, t+1))) + ...
							  lab_demand{2}(s_max+1+stateF(:,t+1));

	end

    month_increment = 4*round(12*1/scale_time/4);


    innovation_mat = zeros(size(stateL, 1), size(stateL, 2), 10, 2); 
    innovation_mat(:, :, 1, 1) = leader_innovL2;
    innovation_mat(:, :, 2, 1) = laggard_innovL2;
    innovation_mat(:, :, 3, 1) = entrant_innovL2;
    innovation_mat(:, :, 4, 1) = leader_innovF2;
    innovation_mat(:, :, 5, 1) = laggard_innovF2;
    innovation_mat(:, :, 6, 1) = entrant_innovF2;
    innovation_mat(:, :, 7, 1) = tied_innovL; 
    innovation_mat(:, :, 8, 1) = tied_innovF;
    innovation_mat(:, :, 9, 1) = entrant_tied_innovL;
    innovation_mat(:, :, 10, 1) = entrant_tied_innovF;
    innovation_mat(:, :, 1, 2) = stateL;
    innovation_mat(:, :, 2, 2) = stateL;
    innovation_mat(:, :, 3, 2) = stateL;
    innovation_mat(:, :, 4, 2) = stateF;
    innovation_mat(:, :, 5, 2) = stateF;
    innovation_mat(:, :, 6, 2) = stateF;
    innovation_mat(:, :, 7, 2) = stateL;
    innovation_mat(:, :, 8, 2) = stateF;
    innovation_mat(:, :, 9, 2) = stateL;
    innovation_mat(:, :, 10, 2) = stateF;

	countL_innov = sum(innovation_mat(:, :, [1:3, 7], 1), 3);
    countF_innov = sum(innovation_mat(:, :, [4:6, 8], 1), 3);
    count2d = [countL_innov; countF_innov]; 

    sales_divisor = sales_func(kappa, s_max, nu_s);
    state2d = [stateL; stateF]; 
    xi2d = [xiL; xiF];
  
    salesL = sales_divisor(s_max+1+state2d);
	fric_binding = fric_binding(s_max+1+state2d);
    profitL = profit(s_max+1+state2d);
    valueL = value_functions(s_max+1+state2d);

    sales2d = ([salesL]/month_increment).*((1+growth_rate).^(0:(size(salesL, 2)-1)));
    profit2d = ([profitL]/month_increment).*((1+growth_rate).^(0:(size(salesL, 2)-1))); 
    value2d = [valueL];
	output2d = salesL / month_increment; 

    count_t = annualize_sim_data(count2d, month_increment, "sum"); 
	lab_annual_t = annualize_sim_data([lab_demandL; lab_demandF], month_increment, "sum");
    sales_annual_t = annualize_sim_data(sales2d, month_increment, "sum");
    profit_annual_t = annualize_sim_data(profit2d, month_increment, "sum");   
    value_t_eop = annualize_sim_data(value2d, month_increment, "EoP");
	prod_t_eop = annualize_sim_data([prodL; prodF], month_increment, "EoP");
    xi_annual_t = annualize_sim_data(xi2d, month_increment, "sum");
    sales_annual_tm1 = [nan(size(sales_annual_t, 1), 1), sales_annual_t(:, 1:(end-1))];
    profit_annual_tm1 = [nan(size(profit_annual_t, 1), 1), profit_annual_t(:, 1:(end-1))]; 
    sales_annual_t1 = [sales_annual_t(:, 2:end), nan(size(sales_annual_t, 1), 1)];
	prod_annual_t1 = [prod_t_eop(:, 2:end), nan(size(sales_annual_t, 1), 1)]; 
    profit_annual_t1 = [profit_annual_t(:, 2:end), nan(size(profit_annual_t, 1), 1)]; 
	lab_annual_t1 = [lab_annual_t(:, 2:end), nan(size(lab_annual_t, 1), 1)]; 
	output_annual_t = annualize_sim_data(output2d, month_increment, "sum"); 

	% Did the constraint bind for most of the year?
	fric_annual_t = annualize_sim_data(fric_binding, month_increment, "sum");
	fric_annual_t = (fric_annual_t >= 0.5*month_increment); 

    
    timevar = reshape(repmat(1:size(sales_annual_t, 2), [size(sales_annual_t, 1), 1]), [], 1);
    firm = repmat([(1:N)'; ((N+1):(2*N))'], [1, size(sales_annual_t, 2)]);
    industry = repmat([(1:N)'; (1:N)'], [1, size(sales_annual_t, 2)]);

    firm = reshape(firm, [], 1);
    industry = reshape(industry, [], 1);

    sale_t = reshape(sales_annual_t, [], 1);
	fric_t = reshape(fric_annual_t, [], 1); 
	output_t = reshape(output_annual_t, [], 1); 
    sale_tm1 = reshape(sales_annual_tm1, [], 1);
    sale_t1 = reshape(sales_annual_t1, [], 1);
	prod_t1 = reshape(prod_annual_t1, [], 1); 
    profit_t = reshape(profit_annual_t, [], 1);
    profit_tm1 = reshape(profit_annual_tm1, [], 1);
    profit_t1 = reshape(profit_annual_t1, [], 1);
	lab_t1= reshape(lab_annual_t1, [], 1); 
	lab_t = reshape(lab_annual_t, [], 1);
    value_t_eop = reshape(value_t_eop, [], 1); 
	prod_t_eop = reshape(prod_t_eop, [], 1); 
    xi_annual_t = reshape(xi_annual_t, [], 1); 
    count_t_vec = reshape(count_t, [], 1);

    sales_reshape = reshape(sales2d', [1, month_increment/4, 4, size(sales2d, 2)/month_increment, size(sales2d, 1)]);
    profit_reshape = reshape(profit2d', [1, month_increment/4, 4, size(sales2d, 2)/month_increment, size(sales2d, 1)]);
    
	profit_reshape = sum(profit_reshape, 2);
	sales_reshape = sum(sales_reshape, 2);

	profit_reshape = reshape(profit_reshape, size(profit_reshape,3)*size(profit_reshape, 4), size(profit_reshape, 5))';
	profit_reshape_t1 = [profit_reshape(:, 2:end), nan(size(profit_reshape, 1), 1)];
	profit_reshape_tm1 = [nan(size(profit_reshape, 1), 1), profit_reshape(:, 1:(end-1))];
	sales_reshape = reshape(sales_reshape,  size(sales_reshape,3)*size(sales_reshape, 4), size(sales_reshape, 5))';
	qtrmat = repmat(1:4, [size(profit_reshape, 1), size(profit_reshape, 2)/4]);
	yearmat = repmat(min(timevar(:, 1)):max(timevar(:, 1)), [size(profit_reshape, 1), 4]);
	[~, sortidx] = sort(yearmat(1, :));
	yearmat = yearmat(:, sortidx);
	firmmat = repmat(unique(firm), [1, size(qtrmat, 2)]);

	qtrpanel = array2table([reshape(profit_reshape, [], 1), reshape(sales_reshape, [], 1), reshape(qtrmat, [], 1), reshape(yearmat, [], 1), ...
							reshape(firmmat, [], 1), reshape(profit_reshape_t1, [], 1), reshape(profit_reshape_tm1, [], 1)], ...
						   'VariableNames', ["Profit_q", "Sale_q", "Quarter", "Time", "Firm", "Profit_q_t1", "Profit_q_tm1"]);


    G = @(x,B,gamma,subsidy) wage_share*(1-subsidy).*(x./B).^(1./gamma); 
    % xrd = (G(investment(s_max+1+state2d), B, gamma, subsidy(s_max+1+state2d))/month_increment).*((1+growth_rate).^(0:(size(salesL, 2)-1))); 
    xrd = (G(investment(s_max+1+state2d), B, gamma, 0)/month_increment).*((1+growth_rate).^(0:(size(salesL, 2)-1))); 
	xrd_ann = annualize_sim_data(xrd, month_increment, "sum");
	xrd_ann_t1 = [xrd_ann(:, 2:end), nan(size(xrd_ann, 1), 1)]; 
    rdsale = reshape(xrd_ann./sales_annual_t, [], 1); 
    rdsale(isinf(rdsale)) = nan;
    xrd_ann_pivot = reshape(xrd_ann, [], 1);
	xrd_ann_pivot_t1 = reshape(xrd_ann_t1, [], 1); 


    panel = array2table([timevar, firm, industry, fric_t, sale_t, sale_tm1, sale_t1, prod_t1, lab_t1, profit_t, profit_tm1, profit_t1, ...
                         value_t_eop, prod_t_eop, lab_t, xi_annual_t, rdsale, xrd_ann_pivot, xrd_ann_pivot_t1, count_t_vec, output_t], ...
                        'VariableNames', ...
                        {'Time', 'Firm', 'Industry', 'Friction_t', 'Sale_t', 'Sale_tm1', 'Sale_t1', 'Prod_t1', 'Lab_t1', 'Profit_t', 'Profit_tm1', 'Profit_t1', ...
                        'Value_t', 'Prod_t', 'Lab_t', 'xi_t', 'rdsales', 'xrd', 'xrd_t1', 'count_patent', 'Output_t'});
	
	if (strcmp(qtrtype, "threeqtr"))

		threeqtr = groupsummary(qtrpanel(:, ["Firm", "Time", "Profit_q"]), ["Firm", "Time"], @(x) sum(x > 0) >= 3); 
		threeqtr = renamevars(threeqtr, "fun1_Profit_q", "two_qtr_restriction");

		twoqtrjoin = threeqtr;
	elseif (strcmp(qtrtype, "0qtr"))

		threeqtr = groupsummary(qtrpanel(:, ["Firm", "Time", "Profit_q"]), ["Firm", "Time"], @(x) sum(x > 0) >= 0); 
		threeqtr = renamevars(threeqtr, "fun1_Profit_q", "two_qtr_restriction");

		twoqtrjoin = threeqtr;
	end
	panel = outerjoin(panel, twoqtrjoin, 'Type', 'left', 'MergeKeys', true, 'keys', ["Time", "Firm"]);
	

    panel.theta = panel.xi_t./panel.Value_t;
    panel.profitgrowth = panel.Profit_t1./panel.Profit_t - 1;
    panel.profitgrowth(isinf(panel.profitgrowth)) = nan;
	panel.profitgrowth_halti = (panel.Profit_t1 - panel.Profit_t) ./ (0.5*(panel.Profit_t1 + panel.Profit_t)); 
    panel.salegrowth = panel.Sale_t1./panel.Sale_t - 1; 
	panel.salegrowth_halti = (panel.Sale_t1 - panel.Sale_t) ./ (0.5*(panel.Sale_t1 + panel.Sale_t)); 
	panel.rdgrowth = panel.xrd_t1 ./ panel.xrd - 1; 
	panel.rdgrowth_halti = (panel.xrd_t1 - panel.xrd) ./ (0.5*(panel.xrd_t1 + panel.xrd) ); 
	panel.xrdgrowth = panel.xrd_t1./panel.xrd - 1; 
    panel.salegrowth(isinf(panel.salegrowth)) = nan;
	panel.prodgrowth = panel.Prod_t1 ./ panel.Prod_t - 1; 
	panel.empgrowth = (panel.Lab_t1 - panel.Lab_t) ./ (0.5*(panel.Lab_t1 + panel.Lab_t)); 


	panel.salegrowth_halti(panel.rdgrowth < winsor_vec_drop(1) | panel.rdgrowth > winsor_vec_drop(2)) = nan; 
	panel.rdgrowth_halti(panel.rdgrowth < winsor_vec_drop(1) | panel.rdgrowth > winsor_vec_drop(2)) = nan; 
	panel.rdgrowth(panel.rdgrowth < winsor_vec_drop(1) | panel.rdgrowth > winsor_vec_drop(2)) = nan; 
	panel.profitgrowth(panel.profitgrowth < winsor_vec_drop(1) | panel.profitgrowth > winsor_vec_drop(2)) = nan;
	panel.salegrowth(panel.salegrowth < winsor_vec_drop(1) | panel.salegrowth > winsor_vec_drop(2)) = nan;
	panel.rdsales(panel.rdsales < winsor_vec_drop(1) | panel.rdsales > winsor_vec_drop(2)) = nan;
	panel.prodgrowth(panel.prodgrowth < winsor_vec_drop(1) | panel.prodgrowth > winsor_vec_drop(2)) = nan;



	condition_uniform =  ~(panel.xrd > 0 & panel.Sale_t > 0 & panel.Profit_t > 0 & ...
						   ~isnan(panel.rdsales) & ~isnan(panel.profitgrowth) & ~isnan(panel.salegrowth) & ...
						   ~isnan(panel.Profit_t) & ~isnan(panel.Profit_tm1) & ~isnan(panel.Profit_t1));

	condition_uniform = condition_uniform | (~panel.two_qtr_restriction);

	panel_save = panel; 
	panel(condition_uniform, :) = [];

	panel = bucket_quantile(panel, "Profit_t", nbucket);
	panel_2 = bucket_quantile(panel, "Profit_t", 2); 
	panel_patent = bucket_quantile(panel, "Profit_t", nbucket_patent);
	panel_2 = renamevars(panel_2, ["bucket", "bucket_tm1", "bucket_t1"], ["bucket_2", "bucket_tm1_2", "bucket_t1_2"]); 
	panel_patent = renamevars(panel_patent, ["bucket", "bucket_tm1", "bucket_t1"], ["bucket_patent", "bucket_tm1_patent", "bucket_t1_patent"]);

	panel = outerjoin(panel, panel_2(:, ["Time", "Firm", "Industry", "bucket_2", "bucket_tm1_2", "bucket_t1_2"]), ...
			    	  'Type', 'left', 'MergeKeys', true, 'keys', ["Time", "Firm", "Industry"]); 
	panel = outerjoin(panel, panel_patent(:, ["Time", "Firm", "Industry", "bucket_patent", "bucket_tm1_patent", "bucket_t1_patent"]), ...
					  'Type', 'left', 'MergeKeys', true, 'keys', ["Time", "Firm", "Industry"]);

	panel_p90 = groupsummary(panel(:, ["Time", "Profit_t", "Profit_t1"]), "Time", @(x) prctile(x, 90));
	panel_p90 = renamevars(panel_p90, [strcat("fun1_", "Profit_t"), strcat("fun1_", "Profit_t1")], ["p90_t", "p90_t1"]);
	panel_p90 = outerjoin(panel, panel_p90, 'Type', 'left', 'MergeKeys', true, 'keys', ["Time"]);
	panel_p90.top10_t = double(panel_p90.Profit_t > panel_p90.p90_t); 
	panel_p90.top10_t1 = double(panel_p90.Profit_t1 > panel_p90.p90_t1);
	panel_p90 = panel_p90(:, ["Firm", "Time", "Industry", "top10_t", "top10_t1"]);

	panel = outerjoin(panel, panel_p90, 'Type', 'left', 'MergeKeys', true, 'keys', ["Firm", "Time", "Industry"]); 


	panel_p80 = groupsummary(panel(:, ["Time", "Profit_t", "Profit_t1"]), "Time", @(x) prctile(x, 80));
	panel_p80 = renamevars(panel_p80, [strcat("fun1_", "Profit_t"), strcat("fun1_", "Profit_t1")], ["p80_t", "p80_t1"]);
	panel_p80 = outerjoin(panel, panel_p80, 'Type', 'left', 'MergeKeys', true, 'keys', ["Time"]);
	panel_p80.top20_t = double(panel_p80.Profit_t > panel_p80.p80_t); 
	panel_p80.top20_t1 = double(panel_p80.Profit_t1 > panel_p80.p80_t1);
	panel_p80 = panel_p80(:, ["Firm", "Time", "Industry", "top20_t", "top20_t1"]);


	panel = outerjoin(panel, panel_p80, 'Type', 'left', 'MergeKeys', true, 'keys', ["Firm", "Time", "Industry"]); 
	panel = renamevars(panel, 'Prod_t', 'Productivity_index'); 
	panel_save = renamevars(panel_save, 'Prod_t', 'Productivity_index'); 


    m.panel = panel;
    m.panel_save = panel_save; 
    m.month_increment = month_increment;
    m.innovation_mat = innovation_mat; 
    m.scale_time = scale_time;
    m.leader_innovL = leader_innovL;
    m.laggard_innovL = laggard_innovL;
    m.entrant_innovL = entrant_innovL; 
    m.leader_innovF = leader_innovF;
    m.laggard_innovF = laggard_innovF; 
    m.entrant_innovF = entrant_innovF; 
    m.ageL = ageL;
    m.ageF = ageF;


end