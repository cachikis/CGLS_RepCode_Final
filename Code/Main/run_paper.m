%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: run_paper.m
% Author: Craig A. Chikis
% Date: 10/19/2023
% Note(s):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = run_paper()
    % Run data if you have it (if not, what you need is in Store_Data/)
    % status code == 0 means that it sucessfully ran; otherwise, it failed
    % [status, ~] = system("Rscript Code/Data/data_filter_v4.R"); 

    % Benchmark
    main_run();
    % Entry 
    entry_new_run();
    % Identification
    moment_sensitivity_new(); 
    % Financial frictions
    finfric_sec(); 
    % Financial frictions scatter
    finfric_scatter(); 
    % Vanishing QCU
    alpha_iterate(); 
    % Transition
    transition_solve(); 
    % Robustness 
    robust_out(); 
    % Error figure
    error_figure(); 
    % Multiplier
    multiplier_new();
    % IA paper on CFG facts
    data_figure_IA();  
    % Monetary policy shock 
    mps(); 

end