function Net_WP = CRsubj_net(mean_data, thresh)
% mean_data is the outputs of mean_data_cross_trials
% select proper subjects for network construction
% if there are multiple levels of temperature in each condition, randomly
% select one.

% input: mean_data: struct varible from mean_data_cross_trials.m
%        thresh: temperature threshold of warm and pain (we choose 45.3¡æ)

Nstudies = length(mean_data.study_name);

% Grouping
%--------------------------------------------------------------------------
k = 1; 

Wdata_for_net = []; Pdata_for_net = [];
Wdata_ori = []; Pdata_ori = []; 

for i = 1:Nstudies
    
    for j = 1:length(mean_data.temp_level{i}) % number of temperature level
        
        temp_tmp = mean_data.temp_level{i}{j};       
        activity_ori_tmp = mean_data.activity_ori_mean_mad{i}{j};
        activity_r_tmp = mean_data.activity_mean_r_mad{i}{j};
        rating_tmp = mean_data.rating{i}{j};
                
        W_idx_tmp = find(temp_tmp < thresh);
        P_idx_tmp = find(temp_tmp > thresh);
               
        if ~isempty(W_idx_tmp) && ~isempty(P_idx_tmp) % if both conditions have more than one level of temperatures
                
            r1 = randperm(length(W_idx_tmp)); r1 = W_idx_tmp(r1(1)); % randomly select 1 warm level
            r2 = randperm(length(P_idx_tmp)); r2 = P_idx_tmp(r2(1)); % randomly select 1 painful level            
            
            % Warm
            A = activity_ori_tmp(r1, :); % randomly select one temperature in warm
            A1 = activity_r_tmp(r1, :);
            temp_warm(k,:) = temp_tmp(r1);
            rating_warm(k,:) = rating_tmp(r1);
            
            % Painful
            B = activity_ori_tmp(r2,:);
            B1 = activity_r_tmp(r2,:);
            temp_pain(k,:) = temp_tmp(r2);
            rating_pain(k,:) = rating_tmp(r2);           
               
         else
                
             continue;
                              
         end
                            
        Wdata_ori = [Wdata_ori; A]; clear A;
        Wdata_for_net = [Wdata_for_net; A1]; clear A1;
        
        Pdata_ori = [Pdata_ori; B]; clear B;
        Pdata_for_net = [Pdata_for_net; B1]; clear B1;
         
        k = k + 1;
        
    end
      
end


%% Network construction (Pearson correlation)
% -------------------------------------------------------------------------
Nnodes = size(Wdata_for_net, 2);
Warm_net = zeros(Nnodes, Nnodes); Pain_net = zeros(Nnodes, Nnodes);
Warm_p = zeros(Nnodes, Nnodes); Pain_p = zeros(Nnodes, Nnodes);

for i = 1:Nnodes
    
    for j = (i+1):Nnodes
        
        [Warm_net(i, j), Warm_p(i, j)] = corr(Wdata_for_net(:, i), Wdata_for_net(:, j));
        [Pain_net(i, j), Pain_p(i, j)] = corr(Pdata_for_net(:, i), Pdata_for_net(:, j));
        
    end
        
end

Warm_net = Warm_net + Warm_net'; Pain_net = Pain_net + Pain_net';
Warm_p = Warm_p + Warm_p'; Pain_p = Pain_p + Pain_p';

Net_WP.Warm.data_for_net = Wdata_for_net;
Net_WP.Warm.data_mad = Wdata_ori;
Net_WP.Warm.r_mat = Warm_net;
Net_WP.Warm.p_mat = Warm_p;
Net_WP.Warm.temp = temp_warm;
Net_WP.Warm.rating = rating_warm;

Net_WP.Pain.data_for_net = Pdata_for_net; 
Net_WP.Pain.data_mad = Pdata_ori;
Net_WP.Pain.r_mat = Pain_net;
Net_WP.Pain.p_mat = Pain_p;
Net_WP.Pain.temp = temp_pain;
Net_WP.Pain.rating = rating_pain;

end
