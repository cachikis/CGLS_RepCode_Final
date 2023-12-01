%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: run_benchmark.m
% Author: Craig A. Chikis
% Date: 09/30/2022
% Note(s): 
% MA-MFA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = run_benchmark(rho, gamma, x, str) 
    

	
    % m = benchmark_current();
    m = benchmark_new(); 
    if (nargin > 0)

        m.rho = rho;
        m.gamma = gamma;
        m.B = x(1);
        m.lambda = x(2);
        m.phi_wt = x(3);

        [prob, ~, ~] = F_mat_main(m.phi_wt, m.phi_wt_tilde, m.phi_wt_exog, m.l, m.l_tilde, m.s_max);

        m.prob = prob; 
    else
        str= "benchmark";
    end

    m = workhorse_nested_robust(m);

    [m.p50_approx, ~, ~] =  markup_quant(m.kappa, m.lambda, m.s_max, m.firm_distributions, m.nu_s, 0.5);
    [m.p90_approx, m.markup_array, ~] =  markup_quant(m.kappa, m.lambda, m.s_max, m.firm_distributions, m.nu_s, 0.9);
    [m.lerner_p50_approx, ~, ~] =  lerner_quant(m.kappa, m.lambda, m.s_max, m.firm_distributions, m.nu_s, 0.5);
    [m.lerner_p90_approx, m.lerner_array, ~] =  lerner_quant(m.kappa, m.lambda, m.s_max, m.firm_distributions, m.nu_s, 0.9);

    lerner_mean = sum(m.lerner_array(1, :).*m.lerner_array(2, :));


    prctile_inp_vec = [90, 80]; 
    m_vec = cell(1, length(prctile_inp_vec));
    for (ii = 1:length(prctile_inp_vec))

        m_tmp = m; 
        m_tmp.prctile_inp = prctile_inp_vec(ii); 
        
        m_tmp = sim_wrap(m_tmp);
        m_tmp = significance_construct_sim(m_tmp); 
        m_tmp = CSTAT_sim(m_tmp); 
        m_tmp = pcm_sim(m_tmp); 

        m_vec{ii} = m_tmp;

    end

                            
    writetable(m_vec{1}.panel_out, strcat("Output/Store_Data/panel_out_", num2str(m.numCPC), "_", str, ".csv"));

    writetable(m_vec{1}.qualityprof_tm1, strcat("Output/Store_Data/qualitysale_tm1_", num2str(m.numCPC), "_", num2str(prctile_inp_vec(1)), "_", str, ".csv")); 
    writetable(m_vec{2}.qualityprof_tm1, strcat("Output/Store_Data/qualitysale_tm1_", num2str(m.numCPC), "_", num2str(prctile_inp_vec(2)), "_", str, ".csv")); 

    writetable(m_vec{1}.citation_track2, strcat("Output/Store_Data/citation_track2_", num2str(m.numCPC), "_", str, ".csv"));


    params = [m.phi_wt, m.B, m.lambda, m.zeta, m.kappa]';
    names = ["phi", "B", "lambda", "zeta", "kappa"]';
    parameterization = array2table(params, 'VariableNames', ["values"]);
    parameterization.parameter = names; 


    writetable(parameterization, strcat("Output/Store_Data/parameters", num2str(m.numCPC), "_", str, ".csv")); 



end
