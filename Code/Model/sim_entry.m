%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: sim_entry.m
% Author: Craig A. Chikis
% Date: 08/27/2022
% Note(s):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = sim_entry(m)

	N = m.N;
	T = m.T;
	investment = m.investment;  
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



    rng(20220830, 'Threefry') 

    s_max = (length(investment) - 1)/2;

    if ((2*max(investment) + max(eta) + max(investment_entrant)) > 1)
		scale_time = 0.99/(2*max(investment) + max(eta) + max(investmentE));
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

	xiL = zeros(N, T);
	xiF = zeros(N, T);

    lab_demandL = zeros(N,T);
	lab_demandF = zeros(N,T);
    
	init_dist = firm_distributions*N;

	prodL = ones(N, T);
	prodF = ones(N, T); 

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

    lab_demandL(:, 1) = lab_demand{1}(s_max+1+stateL(:,1)) + ...
                        lab_demand{2}(s_max+1+stateF(:,1));
    
    lab_demandF(:, 1) = lab_demand{1}(s_max+1+stateF(:,1)) + ...
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

							ageL(n,t+1) = ageL(n,t)+1;
							ageF(n,t+1) = ageF(n,t)+1;

							xiL(n, t+1) = value_functions(find(column_index == stateL(n, t+1))) - ...
										  value_functions(find(column_index == stateL(n, t)));

							leader_innovL(n,t+1) = 1;
							leader_innovL2(n,t+1) = 1;

							prodL(n,t+1) = prodL(n,t) * m.lambda ^ (stateL(n,t+1) - stateL(n,t));
							prodF(n,t+1) = prodF(n,t);



						elseif (draws(n,t) >= (start_pt + L_outcome_o(abs(stateL(n,t))+1, k))) && (draws(n,t) < (start_pt+L_outcome_o(abs(stateL(n,t))+1, k)+F_outcome_o(abs(stateL(n,t))+1, k)))
							stateF(n,t+1) = column_index(k);
							stateL(n,t+1) = -column_index(k);

							ageL(n,t+1) = ageL(n,t)+1;
							ageF(n,t+1) = ageF(n,t)+1;

							xiF(n, t+1) = value_functions(find(column_index == stateF(n, t+1))) - ...
										  value_functions(find(column_index == stateF(n, t)));

							laggard_innovF(n,t+1) = 1;
							laggard_innovF2(n,t+1) = 1;

							prodL(n,t+1) = prodL(n,t);
							prodF(n,t+1) = prodF(n,t)*m.lambda^(stateF(n,t+1) - stateF(n,t)); 


						elseif (draws(n,t) >= (start_pt+L_outcome_o(abs(stateL(n,t))+1, k)+F_outcome_o(abs(stateL(n,t))+1, k))) && ...
								(draws(n,t) < (start_pt+L_outcome_o(abs(stateL(n,t))+1, k)+F_outcome_o(abs(stateL(n,t))+1, k)+E_outcome_o(abs(stateL(n,t))+1, k)))

							stateF(n,t+1) = column_index(k);
							stateL(n,t+1) = -column_index(k);


							ageF(n,t+1) = 0;
							ageL(n,t+1) = ageL(n,t)+1;


							entrant_innovF(n,t+1) = 1;
							entrant_innovF2(n,t+1) = 1;


							xiF(n,t+1) = nan;

							prodL(n,t+1) = prodL(n,t); 
							prodF(n,t+1) = prodF(n,t)*m.lambda^(stateF(n,t+1) - stateF(n,t)); 


						elseif (draws(n,t) >= (start_pt+L_outcome_o(abs(stateL(n,t))+1, k)+F_outcome_o(abs(stateL(n,t))+1, k)+E_outcome_o(abs(stateL(n,t))+1, k))) && ...
								(draws(n,t) < (start_pt+L_outcome_o(abs(stateL(n,t))+1, k)+F_outcome_o(abs(stateL(n,t))+1, k)+E_outcome_o(abs(stateL(n,t))+1, k)+exog_outcome_o(abs(stateL(n,t))+1, k)))

							stateF(n,t+1) = column_index(k);
							stateL(n,t+1) = -column_index(k);


							ageL(n,t+1) = ageL(n,t)+1;
							ageF(n,t+1) = ageF(n,t)+1;

							laggard_innovF(n,t+1) = 1;

							prodL(n,t+1) = prodL(n,t);
							prodF(n,t+1) = prodF(n,t);
							

						end

						start_pt = start_pt + L_outcome_o(abs(stateL(n,t))+1, k) + F_outcome_o(abs(stateL(n,t))+1, k) + E_outcome_o(abs(stateL(n,t))+1, k) + exog_outcome_o(abs(stateL(n,t))+1, k);
					end

				elseif (stateL(n, t) < 0)
					start_pt = nothing_o(abs(stateL(n,t))+1);
					for (k = 1:(2*s_max+1))
						if (draws(n, t) >= start_pt) && (draws(n,t) < (start_pt + L_outcome_o(abs(stateL(n,t))+1, k)))

							stateF(n, t+1) = column_index(k);
							stateL(n,t+1) = -column_index(k);

							ageL(n,t+1) = ageL(n,t)+1;
							ageF(n,t+1) = ageF(n,t)+1;

							xiF(n, t+1) = value_functions(find(column_index == stateF(n, t+1))) - ...
										  value_functions(find(column_index == stateF(n, t)));

							leader_innovF(n,t+1) = 1;
							leader_innovF2(n,t+1) = 1;

							prodL(n,t+1) = prodL(n,t);
							prodF(n,t+1) = prodF(n,t)*m.lambda^(stateF(n,t+1) - stateF(n,t)); 


						elseif (draws(n,t) >= (start_pt + L_outcome_o(abs(stateL(n,t))+1, k))) && (draws(n,t) < (start_pt+L_outcome_o(abs(stateL(n,t))+1, k)+F_outcome_o(abs(stateL(n,t))+1, k)))
							stateL(n,t+1) = column_index(k);
							stateF(n,t+1) = -column_index(k);

							ageL(n,t+1) = ageL(n,t)+1;
							ageF(n,t+1) = ageF(n,t)+1;

							xiL(n, t+1) = value_functions(find(column_index == stateL(n,t+1))) - ...
										  value_functions(find(column_index == stateL(n,t)));

							laggard_innovL(n,t+1) = 1;
							laggard_innovL2(n,t+1) = 1;

							prodL(n,t+1) = prodL(n,t) * m.lambda^(stateL(n,t+1) - stateL(n,t)); 
							prodF(n,t+1) = prodF(n,t);



						elseif (draws(n,t) >= (start_pt+L_outcome_o(abs(stateL(n,t))+1, k)+F_outcome_o(abs(stateL(n,t))+1, k))) && ...
								(draws(n,t) < (start_pt+L_outcome_o(abs(stateL(n,t))+1, k)+F_outcome_o(abs(stateL(n,t))+1, k)+E_outcome_o(abs(stateL(n,t))+1, k)))

							stateL(n,t+1) = column_index(k);
							stateF(n,t+1) = -column_index(k);

							ageF(n,t+1) = ageF(n,t)+1;
							ageL(n,t+1) = 0;


							entrant_innovL(n,t+1) = 1;
							entrant_innovL2(n,t+1) = 1;


							xiL(n,t+1) = nan;

							prodL(n,t+1) = prodL(n,t)*m.lambda^(stateL(n,t+1) - stateL(n,t));
							prodF(n,t+1) = prodF(n,t); 
							

						elseif (draws(n,t) >= (start_pt+L_outcome_o(abs(stateL(n,t))+1, k)+F_outcome_o(abs(stateL(n,t))+1, k)+E_outcome_o(abs(stateL(n,t))+1, k))) && ...
								(draws(n,t) < (start_pt+L_outcome_o(abs(stateL(n,t))+1, k)+F_outcome_o(abs(stateL(n,t))+1, k)+E_outcome_o(abs(stateL(n,t))+1, k)+exog_outcome_o(abs(stateL(n,t))+1, k)))

							stateL(n,t+1) = column_index(k);
							stateF(n,t+1) = -column_index(k);

							ageL(n,t+1) = ageL(n,t)+1;
							ageF(n,t+1) = ageF(n,t)+1;

							laggard_innovL(n,t+1) = 1;
							
							prodL(n,t+1)=  prodL(n,t);
							prodF(n,t+1) = prodF(n,t); 



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

							ageL(n,t+1) = ageL(n,t)+1;
							ageF(n,t+1) = ageF(n,t)+1;

							xiF(n,t+1) = value_functions(find(column_index == stateF(n,t+1))) - ...
										 value_functions(find(column_index == stateF(n,t)));

							leader_innovF(n,t+1) = 1;
							tied_innovF(n, t+1) = 1;

							prodL(n,t+1) = prodL(n,t);
							prodF(n,t+1) = prodF(n,t)*m.lambda^(stateF(n,t+1) - stateF(n,t)); 

						elseif (draws(n,t) >= (start_pt + L_outcome_o(abs(stateL(n,t))+1, k))) && (draws(n,t) < (start_pt+L_outcome_o(abs(stateL(n,t))+1, k)+F_outcome_o(abs(stateL(n,t))+1, k)))
							
							stateL(n,t+1) = column_index(k);
							stateF(n,t+1) = -column_index(k);

							ageL(n,t+1) = ageL(n,t)+1;
							ageF(n,t+1) = ageF(n,t)+1;

							xiL(n,t+1) = value_functions(find(column_index == stateL(n,t+1))) - ...
										 value_functions(find(column_index == stateL(n,t)));

							leader_innovL(n,t+1) = 1;
							tied_innovL(n, t+1) = 1;

							prodL(n,t+1) = prodL(n,t)*m.lambda^(stateL(n,t+1) - stateL(n,t));
							prodF(n,t+1) = prodF(n,t); 

						elseif (draws(n,t) >= (start_pt+L_outcome_o(abs(stateL(n,t))+1, k)+F_outcome_o(abs(stateL(n,t))+1, k))) && ...
								(draws(n,t) < (start_pt+L_outcome_o(abs(stateL(n,t))+1, k)+F_outcome_o(abs(stateL(n,t))+1, k)+E_outcome_o(abs(stateL(n,t))+1, k)))
							
							if (coinflip_num)
								stateL(n,t+1) = column_index(k);
								stateF(n,t+1) = -column_index(k);

								ageF(n,t+1) = ageF(n,t)+1;
								ageL(n,t+1) = 0;


								entrant_innovL(n,t+1) = 1;
								entrant_tied_innovL(n, t+1) = 1;

								xiL(n,t+1) = nan;

								prodL(n,t+1) = prodL(n,t)*m.lambda^(stateL(n,t+1) - stateL(n,t));
								prodF(n,t+1) = prodF(n,t);

							else
								stateF(n,t+1) = column_index(k);
								stateL(n,t+1) = -column_index(k);

								ageF(n,t+1) = 0;
								ageL(n,t+1) = ageL(n,t)+1;


								entrant_innovF(n,t+1) = 1;
								entrant_tied_innovF(n, t+1) = 1;

								xiF(n,t+1) = nan;

								prodL(n,t+1) = prodL(n,t); 
								prodF(n,t+1) = prodF(n,t)*m.lambda^(stateF(n,t+1) - stateF(n,t)); 
							end
						
						elseif (draws(n,t) >= (start_pt+L_outcome_o(abs(stateL(n,t))+1, k)+F_outcome_o(abs(stateL(n,t))+1, k)+E_outcome_o(abs(stateL(n,t))+1, k))) && ...
								(draws(n,t) < (start_pt+L_outcome_o(abs(stateL(n,t))+1, k)+F_outcome_o(abs(stateL(n,t))+1, k)+E_outcome_o(abs(stateL(n,t))+1, k)+exog_outcome_o(abs(stateL(n,t))+1, k)))

							stateL(n,t+1) = column_index(k);
							stateF(n,t+1) = -column_index(k);

							ageL(n,t+1) = ageL(n,t)+1;
							ageF(n,t+1) = ageF(n,t)+1;

							prodL(n,t+1) = prodL(n,t);
							prodF(n,t+1) = prodF(n,t);

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

		lab_demandL(:, t+1) = lab_demand{1}(s_max+1+stateL(:,t+1)) + ...
							  lab_demand{2}(s_max+1+stateL(:,t+1));
		lab_demandF(:, t+1) = lab_demand{1}(s_max+1+stateF(:,t+1)) + ...
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


    normL = sum(innovation_mat(:, :, [1,2,7], 1), 3);
    normF = sum(innovation_mat(:, :, [4,5,8], 1), 3); 
    entL = sum(innovation_mat(:, :, [3,9], 1), 3);
    entF = sum(innovation_mat(:, :, [6,10], 1), 3);  

    ent2d = reshape([entL; entF], [], 1); 
    norm2d = reshape([normL; normF], [], 1);
    state2d = reshape([stateL; stateF], [], 1); 
	prod2d = reshape([prodL; prodF], [], 1); 
    xi2d = reshape([xiL; xiF], [], 1);
    age2d = reshape([ageL; ageF], [], 1); 

    quarter = reshape(repmat(sort(repmat(1:4, 1, month_increment/4)), size(entL, 1)*2, size(entL, 2)/month_increment), [], 1);
    month = reshape(repmat(1:month_increment, 2*size(entL, 1), size(entL, 2)/month_increment), [], 1); 
    year = reshape(sort(repmat(1:(size(entL,2)/month_increment), 2*size(entL, 1), month_increment), 2), [], 1);
    firm = reshape(repmat([(1:N), ((N+1):(2*N))]', 1, size(entL, 2)), [], 1);
    industry = reshape(repmat([(1:N), (1:N)]', 1, size(entL, 2)), [], 1); 

    panel_ent_month = array2table([ent2d, quarter, month, year, firm, industry],... 
                                  'VariableNames', {'Entry', 'Quarter', 'Month', 'Year', 'Firm', 'Industry'}); 
    panel_norm_month = array2table([norm2d, quarter, month, year, firm, age2d, industry, state2d, ...    
                                    xi2d, prod2d], 'VariableNames', ...
                                   {'Norm', 'Quarter', 'Month', 'Year', 'Firm', 'Age', 'Industry', 'State', 'xi_t', 'Productivity_index'});

    [c1,c2] = meshgrid(unique(panel_norm_month.Year), unique(panel_norm_month.Month));
    cp = array2table([c1(:), c2(:), (0:(size(stateL,2)-1))'], 'VariableNames', ...
                     ["Year", "Month", "factor"]);    

    panel_norm_month = outerjoin(panel_norm_month, cp, 'keys', ["Year", "Month"], 'MergeKeys', true, ...
                                 'Type', 'left');

    panel_norm_month.Profit = (profit(s_max+1+panel_norm_month.State)'/month_increment).*((1 + growth_rate).^...
                                            (panel_norm_month.factor));
    sales_divisor = sales_func(kappa, s_max, nu_s);


	panel_norm_month.Output = (sales_divisor(s_max+1+panel_norm_month.State)' ./ month_increment);

    panel_norm_month.Sale = (sales_divisor(s_max+1+panel_norm_month.State)'/month_increment).*((1 + growth_rate).^...
                                            (panel_norm_month.factor));

   
    G = @(x,B,gamma,subsidy) wage_share*(1-subsidy).*(x./B).^(1./gamma); 
    % panel_norm_month.xrd = ...
    %     (G(investment(s_max+1+panel_norm_month.State)', B, gamma, ...
    %                     subsidy(s_max+1+panel_norm_month.State)')/month_increment).*(...
    %         (1 + growth_rate).^(panel_norm_month.factor));

	panel_norm_month.xrd = ...
	(G(investment(s_max+1+panel_norm_month.State)', B, gamma, ...
					0)/month_increment).*(...
		(1 + growth_rate).^(panel_norm_month.factor));
    panel_norm_month.Value_t = value_functions(s_max+1+panel_norm_month.State)';
    panel_norm_month.emp = lab_demand{1}(s_max+1+panel_norm_month.State)' + ...
                           lab_demand{2}(s_max+1+panel_norm_month.State)'; 

    if (~entrance)
        panel_ent_month(1:size(panel_ent_month, 1), :) = []; 
    end

    replace = outerjoin(panel_norm_month, panel_ent_month, 'keys', ["Firm", "Industry", "Quarter", "Month", "Year"], 'MergeKeys', true, ...
                        'Type', 'left');
    if (~entrance)
        replace.Entry(:) = 0;
    end
	replace.xrd(logical(replace.Entry)) = nan;
	replace.Value_t(logical(replace.Entry)) = nan;
	replace.emp(logical(replace.Entry)) = nan;
	replace.Sale(logical(replace.Entry)) = nan;
	replace.Profit(logical(replace.Entry)) = nan; 
	replace.Output(logical(replace.Entry)) = nan; 

    replace = sortrows(replace, ["Firm", "Year", "Month"]);
    if (entrance)
        uniquefirm = unique(replace.Firm);
        last = uniquefirm(1);
        replace.Firmmodify(:) = nan;
        for (ii = 1:length(uniquefirm))
            replace.Firmmodify(replace.Firm == uniquefirm(ii)) = last + ...
                cumsum(replace.Entry(replace.Firm == uniquefirm(ii)));
        
            last = max(replace.Firmmodify(replace.Firm == uniquefirm(ii))) +1; 
        end
        replace.Firm = replace.Firmmodify;
        replace.Firmmodify = [];
    end


    if (strcmp(qtrtype, "threeqtr"))
        key_stay = groupsummary(replace(:, ["Firm", "Industry", "Year", "Quarter", "Profit"]), ...
                                ["Firm", "Industry", "Year", "Quarter"], ...
                                @(x) nansum(x));

        key_stay = groupsummary(key_stay    (:, ["Firm", "Industry", "Year", "fun1_Profit"]), ["Firm", "Industry", "Year"], ...
                                @(x) nansum(x > 0) >= 3);
        key_stay= renamevars(key_stay, "fun1_fun1_Profit", "two_qtr_restriction");
    end

    panel_annual_sum = groupsummary(replace(:, ["Firm", "Industry", "Year", "Month", "Sale", "Output", "Profit", "emp", "xrd", "xi_t"]), ...
                                   ["Firm", "Industry", "Year"], ...
                                   @(x) nansum(x));

    panel_annual_sum = renamevars(panel_annual_sum, ["fun1_Profit", "fun1_Sale", "fun1_Output", "fun1_emp", "fun1_xrd", "fun1_xi_t"], ...
                                  ["Profit_t", "Sale_t", "Output_t", "emp", "xrd", "xi_t"]); 
    panel_annual_eop = groupsummary(replace(:, ["Firm", "Industry", "Year", "Month", "Age", "State", "Value_t", "Productivity_index"]), ...
                                   ["Firm", "Industry", "Year"], ...
                                   @(x,y) unique(x(y == max(y))), ...
                                   {["Age", "State", "Value_t", "Productivity_index"], ["Month", "Month", "Month", "Month"]});

    panel_annual_eop = renamevars(panel_annual_eop, ["fun1_Age_Month", "fun1_State_Month", "fun1_Value_t_Month", "fun1_Productivity_index_Month"], ...
                                  ["Age", "State_t", "Value_t", "Productivity_index"]); 
    

    panel_annual_tm1 = panel_annual_sum;
    panel_annual_tm1.Year = panel_annual_tm1.Year + 1; 
    panel_annual_tm1 = renamevars(panel_annual_tm1, ["Profit_t", "Sale_t", "Output_t", "xi_t", "emp"], ...
                                  ["Profit_tm1", "Sale_tm1", "Output_tm1", "xi_tm1", "emp_tm1"]);
    panel_annual_t1 = panel_annual_sum;
    panel_annual_t1.Year = panel_annual_t1.Year-1;
    panel_annual_t1 = renamevars(panel_annual_t1, ["Profit_t", "Sale_t", "Output_t", "xi_t", "emp"], ...
                                 ["Profit_t1", "Sale_t1", "Output_t1", "xi_t1", "emp_t1"]); 

    panel_annual_tm1_eop = panel_annual_eop;
    panel_annual_tm1_eop.Year = panel_annual_tm1_eop.Year + 1; 
    panel_annual_tm1_eop = renamevars(panel_annual_tm1_eop, ["State_t", "Age", "Value_t", "Productivity_index"], ...
                                  ["State_tm1", "Age_tm1", "Value_tm1", "Productivity_index_tm1"]);
    panel_annual_t1_eop = panel_annual_eop;
    panel_annual_t1_eop.Year = panel_annual_t1_eop.Year-1;
    panel_annual_t1_eop = renamevars(panel_annual_t1_eop, ["State_t", "Age", "Value_t", "Productivity_index"], ...
                                 ["State_t1", "Age_t1", "Value_t1", "Productivity_index_t1"]);                  

    panel = outerjoin(panel_annual_sum, panel_annual_tm1, 'Type', 'left', ...
                      'MergeKeys', true, 'keys', ["Firm", "Industry", "Year"]);
    panel = outerjoin(panel, panel_annual_t1, 'Type', 'left', ...
                      'MergeKeys', true, 'keys', ["Firm", "Industry", "Year"]);
    panel = outerjoin(panel, panel_annual_eop, 'Type', 'left', ...
                      'MergeKeys', true, 'keys', ["Firm", "Industry", "Year"]);
    panel = outerjoin(panel, panel_annual_tm1_eop, 'Type', 'left', ...
                      'MergeKeys', true, 'keys', ["Firm", "Industry", "Year"]);
    panel = outerjoin(panel, panel_annual_t1_eop, 'Type', 'left', ...
                      'MergeKeys', true, 'keys', ["Firm", "Industry", "Year"]);
    panel = outerjoin(panel, key_stay, 'Type', 'left', ...
                      'MergeKeys', true, 'keys', ["Firm", "Industry", "Year"]);
    panel.Time = panel.Year;
    panel.Age = panel.Age/month_increment;

    panel.Value_t = value_functions(s_max+1+panel.State_t)';
    panel.rdsales = panel.xrd./panel.Sale_t; 

    panel.theta = panel.xi_t./panel.Value_t;
    panel.profitgrowth = panel.Profit_t1./panel.Profit_t - 1;
    panel.profitgrowth(isinf(panel.profitgrowth)) = nan;
    panel.salegrowth = panel.Sale_t1./panel.Sale_t - 1; 
    panel.salegrowth(isinf(panel.salegrowth)) = nan;


	panel.profitgrowth(panel.profitgrowth < winsor_vec_drop(1) | panel.profitgrowth > winsor_vec_drop(2)) = nan;
	panel.salegrowth(panel.salegrowth < winsor_vec_drop(1) | panel.salegrowth > winsor_vec_drop(2)) = nan;
	panel.rdsales(panel.rdsales < winsor_vec_drop(1) | panel.rdsales > winsor_vec_drop(2)) = nan;


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


	
	% tmp = groupsummary(panel_save(:, ["Time", "Productivity_index_t1", "Productivity_index"]), "Time", @(x,y) nanmean(log(y) - log(x)), ... 
	% 			 	   {"Productivity_index", "Productivity_index_t1"}); 








	
end