%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: transmat_sim.m
% Author: Craig A. Chikis
% Date: 08/27/2022
% Note(s):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = transmat_sim(m)

    
    panel = m.panel_out;


    transmat_out = zeros(5,1); 
    fivediv = zeros(5,1); 
    transmat_out_top10 = zeros(2,2);
    transmat_out_top20 = zeros(2,2);  
    transmat_out_halves = zeros(2,2);
    for (ii = 1:size(panel, 1))
        if (panel.bucket(ii) == 5)
            if (panel.bucket_t1(ii) == 1)
                transmat_out(5,1) = transmat_out(5,1) + 1;
            end
            if (ismember(panel.bucket_t1(ii), [1,2,3,4,5]))
                fivediv(5,1) = fivediv(5,1) + 1; 
            end
        elseif (panel.bucket(ii) == 4)
            if (panel.bucket_t1(ii) == 1)
                transmat_out(4,1) = transmat_out(4,1) + 1;
            end
            if (ismember(panel.bucket_t1(ii), [1,2,3,4,5]))
                fivediv(4,1) = fivediv(4,1) + 1; 
            end
        elseif (panel.bucket(ii) == 3)
            if (panel.bucket_t1(ii) == 1)
                transmat_out(3,1) = transmat_out(3,1) + 1;
            end
            if (ismember(panel.bucket_t1(ii), [1,2,3,4,5]))
                fivediv(3,1) = fivediv(3,1) + 1; 
            end
        elseif (panel.bucket(ii) == 2)
            if (panel.bucket_t1(ii) == 1)
                transmat_out(2,1) = transmat_out(2,1) + 1; 
            end
            if (ismember(panel.bucket_t1(ii), [1,2,3,4,5]))
                fivediv(2,1) = fivediv(2,1) + 1; 
            end
        elseif (panel.bucket(ii) == 1)
            if (panel.bucket_t1(ii) == 1)
                transmat_out(1,1) = transmat_out(1,1) + 1;
            end
            if (ismember(panel.bucket_t1(ii), [1,2,3,4,5]))
                fivediv(1,1) = fivediv(1,1) + 1; 
            end
        end
    end
    transmat_out = transmat_out./fivediv; 


    for (ii = 1:size(panel, 1))
        if (panel.top10_t(ii) == 0)
            if (panel.top10_t1(ii) == 1)
                transmat_out_top10(2,1) = transmat_out_top10(2,1) + 1;
            elseif (panel.top10_t1(ii) == 0)
                transmat_out_top10(2,2) = transmat_out_top10(2,2) + 1;
            end
        elseif (panel.top10_t(ii) == 1)
            if (panel.top10_t1(ii) == 1)
                transmat_out_top10(1,1) = transmat_out_top10(1,1) + 1;
            elseif (panel.top10_t1(ii) == 0)
                transmat_out_top10(1,2) = transmat_out_top10(1,2) + 1;
            end
        end
    end
    transmat_out_top10 = transmat_out_top10./sum(transmat_out_top10, 2); 

    for (ii = 1:size(panel, 1))
        if (panel.top20_t(ii) == 0)
            if (panel.top20_t1(ii) == 1)
                transmat_out_top20(2,1) = transmat_out_top20(2,1) + 1;
            elseif (panel.top20_t1(ii) == 0)
                transmat_out_top20(2,2) = transmat_out_top20(2,2) + 1;
            end
        elseif (panel.top20_t(ii) == 1)
            if (panel.top20_t1(ii) == 1)
                transmat_out_top20(1,1) = transmat_out_top20(1,1) + 1;
            elseif (panel.top20_t1(ii) == 0)
                transmat_out_top20(1,2) = transmat_out_top20(1,2) + 1;
            end
        end
    end
    transmat_out_top20 = transmat_out_top20./sum(transmat_out_top20, 2); 


    

    for (ii = 1:size(panel, 1))
        if (panel.bucket_2(ii) == 2)
            if (panel.bucket_t1_2(ii) == 1)
                transmat_out_halves(2,1) = transmat_out_halves(2,1) + 1;
            elseif (panel.bucket_t1_2(ii) == 2)
                transmat_out_halves(2,2) = transmat_out_halves(2,2) + 1;
            end
        elseif (panel.bucket_2(ii) == 1)
            if (panel.bucket_t1_2(ii) == 1)
                transmat_out_halves(1,1) = transmat_out_halves(1,1) + 1;
            elseif (panel.bucket_t1_2(ii) == 2)
                transmat_out_halves(1,2) = transmat_out_halves(1,2) + 1;
            end
        end
    end
    transmat_out_halves = transmat_out_halves./sum(transmat_out_halves, 2);
    
    

       
    m.transmat_out = transmat_out;
    m.transmat_out_top10 = transmat_out_top10;
    m.transmat_out_top20 = transmat_out_top20;
    m.transmat_out_halves = transmat_out_halves; 
end