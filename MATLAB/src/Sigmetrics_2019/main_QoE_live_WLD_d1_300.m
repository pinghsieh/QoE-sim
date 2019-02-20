%% QoE Main Program
clear;
tic;

MIN_TOL = 1e-10;

%% Part 1: Configuration
%network_config = "unreliable_20links_asymmetric_v3";
%network_config = "unreliable_2links_asymmetric";
%network_config = "unreliable_10links_asymmetric_v2";
network_config = "unreliable_10links_asymmetric_heavy_traffic_v3";
%network_config = "unreliable_10links_asymmetric_under_loaded_v3";

%% For policy
policy = "WLD";
%policy = "WRR";
%policy = "WRand";
%policy = "LDF";
%policy = "LDF-no-dummy-packet";
%policy = "EDF";  

% For figure saving
filepath = './figures/';
mode = 'heavy-traffic-d1=300';
%mode = 'under-loaded';

for k=1:5
    dn{1}(k) = 300; % unit: slot
end
for k=6:10
    dn{1}(k) = 600; % unit: slot
end

switch network_config
    case "unreliable_20links_asymmetric"
        config.config_unreliable_20links_asymmetric;
    case "unreliable_20links_asymmetric_v2"
        config.config_unreliable_20links_asymmetric_v2;
    case "unreliable_20links_asymmetric_v3"
        config.config_unreliable_20links_asymmetric_v3;        
    case "unreliable_2links_asymmetric"
        config.config_unreliable_2links_asymmetric; 
    case "unreliable_10links_asymmetric"
        config.config_unreliable_10links_asymmetric;
    case "unreliable_10links_asymmetric_v2"
        config.config_unreliable_10links_asymmetric_v2;  
    case "unreliable_10links_asymmetric_heavy_traffic_v3"
        config.config_unreliable_10links_asymmetric_heavy_traffic_v3;
    case "unreliable_10links_asymmetric_under_loaded_v3"
        config.config_unreliable_10links_asymmetric_under_loaded_v3;           
end

%% Part 2: Initialization
for m=1:N_trials
% For live streaming
Dn_Avg = zeros(N, Ttot);
Un_Avg = zeros(N, Ttot);
Bn_max = Boffset{m};

%% Part 3: Main Operations
for i=1:Run
    InsRate_history = zeros(N,Ttot);
    AvgRate_history = zeros(N,Ttot);  % AvgRate(i,t) denotes average throughput up to time t
    Schedule_history = zeros(N,Ttot);
    Xn_history = zeros(N,Ttot);
    Dn_history = zeros(N,Ttot);
    Bn_history = zeros(N,Ttot);
    Zn_history = zeros(N,Ttot);
    Un_history = zeros(N,Ttot);
    Qn_history = zeros(N,Ttot);
    Xn = zeros(N,1);
    Dn = zeros(N,1);
    Bn = zeros(N,1);
    Zn = zeros(N,1);
    Un = zeros(N,1);
    Qn = (lambda_n.*dn{m}).*ones(N,1);
    
    %% slot-wise update
    for t=1:Ttot
        % Get playback rates: dummy_client == 1 means there is a dummy client
        InsRate = get_channel_rate(N, cdf_pn, channel_rate_vec);   
        
        % Scheduling: 
        % schedule is a N-by-1 boolean vector
        switch policy
            case "WLD"
                schedule = WLD(Zn, beta_n_WLD);
            case "LDF"
                schedule = LDF(Xn);
            case "LDF-no-dummy-packet"
                schedule = LDF_no_dummy(Xn, Qn);
            case "WRand"
                schedule = WRand(lambda_n, pn);            
            case "WRR"
                schedule = WRR(Tn_WRR, t);
            case "EDF"
                schedule = EDF(Qn, lambda_n, frame_interval);                
            case "WPF"
                schedule = WPF(); 
            otherwise
                schedule = WLD(Zn, beta_n_WLD);
        end
        %% Get next-slot state variables
        % Dn = number of outage frames (or equivalently, the number of sropped frames)
        % Qn = AP-side buffer length
        % Bn = client-side buffer length
        % Un = number of delivered dummy packets
        % An = number of delivered valid data packets
        % Xn = An - t*lambda_n
        % Zn = Xn + Un
        switch video_type
            %{
            case "On-demand"
                Xn_next = Xn + schedule.*InsRate; % for on-demand video streaming
                Bn = Bn + (Xn_next - Xn);
                Xn = Xn_next;  
                Zn = Zn + schedule.*InsRate;
            case "Live"
                Xn_next = min(Bn_max, Xn + schedule.*InsRate); % for live streaming
                Un = Un + max(Xn + schedule.*InsRate - Bn_max, 0);
                Bn = Bn + (Xn_next - Xn);
                Xn = Xn_next;
                Zn = Zn + schedule.*InsRate;
            %}
            % for live streaming with packet dropping  
            case "Live-with-drop"
                is_Qn_empty = (Qn < MIN_TOL);
                Xn_next = Xn + schedule.*InsRate.*(is_Qn_empty == 0) - lambda_n;         
                Un_next = Un + schedule.*InsRate.*(is_Qn_empty == 1);
                Qn_next = Qn - schedule.*InsRate.*(is_Qn_empty == 0);
                Bn_next = Bn + schedule.*InsRate.*(is_Qn_empty == 0);
                Zn_next = Zn + schedule.*InsRate - lambda_n;
                Dn_next = Dn;
        end
        if t > 1
            AvgRate_history(:,t) = (AvgRate_history(:,t-1)*(t-1) + InsRate)/t;
        else
            AvgRate_history(:,t) = InsRate;
        end

        if t > 1
            % Playback
            about_to_play = double(rem(t - Dn.*frame_interval, frame_interval) == 1);
            is_Bn_almost_empty = double(Bn < ((lambda_n.*frame_interval) - MIN_TOL));
            about_to_rebuffer = about_to_play.*is_Bn_almost_empty;
            ready_for_playback = about_to_play - about_to_rebuffer;
            Dn_next = Dn + ones(N,1).*about_to_rebuffer;
            Bn_next = max(0, Bn_next - frame_interval.*(lambda_n.*about_to_play)); % If no enough video content, then the buffer would be depleted (packet dropping)
        
            % Get the new frame
            about_to_get_new_frame = double(rem(t, frame_interval) == 1);      
            Qn_next = Qn_next + frame_interval.*(lambda_n.*about_to_get_new_frame); 
        end
        
        % Update history
        Xn_history(:,t) = Xn;
        Dn_history(:,t) = Dn;
        Bn_history(:,t) = Bn;
        Zn_history(:,t) = Zn;
        Un_history(:,t) = Un;
        Qn_history(:,t) = Qn;
        Schedule_history(:,t) = schedule;
        
        % Update state variables
        Xn = Xn_next;
        Dn = Dn_next;
        Bn = Bn_next;
        Zn = Zn_next;
        Un = Un_next;
        Qn = Qn_next;
    end
    
    
    %% Update average history
    Dn_Avg = Dn_Avg + Dn_history;
    Un_Avg = Un_Avg + Un_history;
    
end
Dn_Avg = Dn_Avg./Run;
DAll_Avg = sum(Dn_Avg);
Un_Avg = Un_Avg./Run;
UAll_Avg = sum(Un_Avg);

%% Part 4: Plotting

switch network_config
    case "unreliable_20links_asymmetric"
        step_size = 1;
        tplot_start = 1;
        taxis = tplot_start:step_size:Ttot;
        Dn_Avg_group1 = mean(Dn_Avg(1:10,taxis));
        Dn_Avg_group2 = mean(Dn_Avg(11:20,taxis));
        Dn_Avg_max_group1 = max(Dn_Avg(1:10,taxis));
        Dn_Avg_max_group2 = max(Dn_Avg(11:20,taxis));
        Un_Avg_group1 = mean(Un_Avg(1:10,taxis));
        Un_Avg_group2 = mean(Un_Avg(11:20,taxis));
        Un_Avg_max_group1 = max(Un_Avg(1:10,taxis));
        Un_Avg_max_group2 = max(Un_Avg(11:20,taxis));        
        createfigure_all(taxis, [Dn_Avg_group1; Dn_Avg_group2], policy);
        createfigure_all(taxis, [Dn_Avg_max_group1; Dn_Avg_max_group2], policy);
        createfigure_all_Un(taxis, [Un_Avg_group1; Un_Avg_group2], policy);
        createfigure_all_Un(taxis, [Un_Avg_max_group1; Un_Avg_max_group2], policy);
    case "unreliable_20links_asymmetric_v2"
        step_size = 1;
        tplot_start = 1;
        taxis = tplot_start:step_size:Ttot;
        Dn_Avg_group1 = mean(Dn_Avg(1:10,taxis));
        Dn_Avg_group2 = mean(Dn_Avg(11:20,taxis));
        Dn_Avg_max_group1 = max(Dn_Avg(1:10,taxis));
        Dn_Avg_max_group2 = max(Dn_Avg(11:20,taxis));
        Un_Avg_group1 = mean(Un_Avg(1:10,taxis));
        Un_Avg_group2 = mean(Un_Avg(11:20,taxis));
        Un_Avg_max_group1 = max(Un_Avg(1:10,taxis));
        Un_Avg_max_group2 = max(Un_Avg(11:20,taxis));        
        createfigure_all(taxis, [Dn_Avg_group1; Dn_Avg_group2], policy);
        createfigure_all(taxis, [Dn_Avg_max_group1; Dn_Avg_max_group2], policy);
        createfigure_all_Un(taxis, [Un_Avg_group1; Un_Avg_group2], policy);
        createfigure_all_Un(taxis, [Un_Avg_max_group1; Un_Avg_max_group2], policy);   
   case "unreliable_20links_asymmetric_v3"
        step_size = 1;
        tplot_start = 1;
        taxis = tplot_start:step_size:Ttot;
        Dn_Avg_group1 = mean(Dn_Avg(1:10,taxis));
        Dn_Avg_group2 = mean(Dn_Avg(11:20,taxis));
        Dn_Avg_max_group1 = max(Dn_Avg(1:10,taxis));
        Dn_Avg_max_group2 = max(Dn_Avg(11:20,taxis));
        Un_Avg_group1 = mean(Un_Avg(1:10,taxis));
        Un_Avg_group2 = mean(Un_Avg(11:20,taxis));
        Un_Avg_max_group1 = max(Un_Avg(1:10,taxis));
        Un_Avg_max_group2 = max(Un_Avg(11:20,taxis));        
        createfigure_all(taxis, [Dn_Avg_group1; Dn_Avg_group2], policy);
        createfigure_all(taxis, [Dn_Avg_max_group1; Dn_Avg_max_group2], policy);
        createfigure_all_Un(taxis, [Un_Avg_group1; Un_Avg_group2], policy);
        createfigure_all_Un(taxis, [Un_Avg_max_group1; Un_Avg_max_group2], policy);          
    case "unreliable_2links_asymmetric"
        step_size = 1;
        tplot_start = 1;
        taxis = tplot_start:step_size:Ttot;
        createfigure_all(taxis, Dn_Avg(:,taxis), policy);
        createfigure_all_Un(taxis, Un_Avg(:,taxis), policy);
    case "unreliable_10links_asymmetric"
        step_size = 10;
        tplot_start = 1;
        taxis = tplot_start:step_size:Ttot;
        Dn_Avg_group1 = mean(Dn_Avg(1:5,taxis));
        Dn_Avg_group2 = mean(Dn_Avg(6:10,taxis));
        Dn_Avg_max_group1 = max(Dn_Avg(1:5,taxis));
        Dn_Avg_max_group2 = max(Dn_Avg(6:10,taxis));
        Un_Avg_group1 = mean(Un_Avg(1:5,taxis));
        Un_Avg_group2 = mean(Un_Avg(6:10,taxis));
        Un_Avg_max_group1 = max(Un_Avg(1:5,taxis));
        Un_Avg_max_group2 = max(Un_Avg(6:10,taxis));        
        createfigure_all(taxis, [Dn_Avg_group1; Dn_Avg_group2], policy);
        createfigure_all(taxis, [Dn_Avg_max_group1; Dn_Avg_max_group2], policy);
        createfigure_all_Un(taxis, [Un_Avg_group1; Un_Avg_group2], policy);
        createfigure_all_Un(taxis, [Un_Avg_max_group1; Un_Avg_max_group2], policy);
    case "unreliable_10links_asymmetric_v2"
        step_size = 20;
        tplot_start = 1;
        taxis = tplot_start:step_size:Ttot;
        Dn_Avg_group1 = mean(Dn_Avg(1:5,taxis));
        Dn_Avg_group2 = mean(Dn_Avg(6:10,taxis));
        Dn_Avg_max_group1 = max(Dn_Avg(1:5,taxis));
        Dn_Avg_max_group2 = max(Dn_Avg(6:10,taxis));
        Un_Avg_group1 = mean(Un_Avg(1:5,taxis));
        Un_Avg_group2 = mean(Un_Avg(6:10,taxis));
        Un_Avg_max_group1 = max(Un_Avg(1:5,taxis));
        Un_Avg_max_group2 = max(Un_Avg(6:10,taxis));        
        createfigure_all(taxis, [Dn_Avg_group1; Dn_Avg_group2], policy);
        createfigure_all(taxis, [Dn_Avg_max_group1; Dn_Avg_max_group2], policy);
        createfigure_all_Un(taxis, [Un_Avg_group1; Un_Avg_group2], policy);
        createfigure_all_Un(taxis, [Un_Avg_max_group1; Un_Avg_max_group2], policy);      
    case "unreliable_10links_asymmetric_heavy_traffic_v3"
        step_size = 8000;
        tplot_start = step_size;
        taxis = [1, tplot_start:step_size:Ttot];
        taxis_second = taxis/Nslots_per_second;
        Dn_Avg_group1 = mean(Dn_Avg(1:5,taxis));
        Dn_Avg_group2 = mean(Dn_Avg(6:10,taxis));
        Dn_Avg_all = mean(Dn_Avg(:,taxis));
        Dn_Avg_max_group1 = max(Dn_Avg(1:5,taxis));
        Dn_Avg_max_group2 = max(Dn_Avg(6:10,taxis));
        Dn_Avg_max = max(Dn_Avg(:,taxis));
        Un_Avg_group1 = mean(Un_Avg(1:5,taxis));
        Un_Avg_group2 = mean(Un_Avg(6:10,taxis));
        Un_Avg_all = mean(Un_Avg(:,taxis));
        Un_Avg_max_group1 = max(Un_Avg(1:5,taxis));
        Un_Avg_max_group2 = max(Un_Avg(6:10,taxis));  
        Un_Avg_max = max(Un_Avg(:,taxis));
        createfigure_all(taxis_second, Dn_Avg_group1, policy, mode, filepath, 'group1_Dn_avg');
        createfigure_all(taxis_second, Dn_Avg_group2, policy, mode, filepath, 'group2_Dn_avg');
        createfigure_all(taxis_second, Dn_Avg_all, policy, mode, filepath, 'all_Dn_avg');
        createfigure_all(taxis_second, Dn_Avg_max_group1, policy, mode, filepath, 'group1_Dn_max');
        createfigure_all(taxis_second, Dn_Avg_max_group2, policy, mode, filepath, 'group2_Dn_max');
        createfigure_all(taxis_second, Dn_Avg_max, policy, mode, filepath, 'all_Dn_max');
        createfigure_all_Un(taxis_second, Un_Avg_group1, policy, mode, filepath, 'group1_Un_avg');
        createfigure_all_Un(taxis_second, Un_Avg_group2, policy, mode, filepath, 'group2_Un_avg');
        createfigure_all_Un(taxis_second, Un_Avg_all, policy, mode, filepath, 'all_Un_avg');
        createfigure_all_Un(taxis_second, Un_Avg_max_group1, policy, mode, filepath, 'group1_Un_max');
        createfigure_all_Un(taxis_second, Un_Avg_max_group2, policy, mode, filepath, 'group2_Un_max');
        createfigure_all_Un(taxis_second, Un_Avg_max, policy, mode, filepath, 'all_Un_max');  
    case "unreliable_10links_asymmetric_under_loaded_v3"
        step_size = 8000;
        tplot_start = step_size;
        taxis = [1, tplot_start:step_size:Ttot];
        taxis_second = taxis/Nslots_per_second;
        Dn_Avg_group1 = mean(Dn_Avg(1:5,taxis));
        Dn_Avg_group2 = mean(Dn_Avg(6:10,taxis));
        Dn_Avg_all = mean(Dn_Avg(:,taxis));
        Dn_Avg_max_group1 = max(Dn_Avg(1:5,taxis));
        Dn_Avg_max_group2 = max(Dn_Avg(6:10,taxis));
        Dn_Avg_max = max(Dn_Avg(:,taxis));
        Un_Avg_group1 = mean(Un_Avg(1:5,taxis));
        Un_Avg_group2 = mean(Un_Avg(6:10,taxis));
        Un_Avg_all = mean(Un_Avg(:,taxis));
        Un_Avg_max_group1 = max(Un_Avg(1:5,taxis));
        Un_Avg_max_group2 = max(Un_Avg(6:10,taxis));  
        Un_Avg_max = max(Un_Avg(:,taxis));
        createfigure_all(taxis_second, Dn_Avg_group1, policy, mode, filepath, 'group1_Dn_avg');
        createfigure_all(taxis_second, Dn_Avg_group2, policy, mode, filepath, 'group2_Dn_avg');
        createfigure_all(taxis_second, Dn_Avg_all, policy, mode, filepath, 'all_Dn_avg');
        createfigure_all(taxis_second, Dn_Avg_max_group1, policy, mode, filepath, 'group1_Dn_max');
        createfigure_all(taxis_second, Dn_Avg_max_group2, policy, mode, filepath, 'group2_Dn_max');
        createfigure_all(taxis_second, Dn_Avg_max, policy, mode, filepath, 'all_Dn_max');
        createfigure_all_Un(taxis_second, Un_Avg_group1, policy, mode, filepath, 'group1_Un_avg');
        createfigure_all_Un(taxis_second, Un_Avg_group2, policy, mode, filepath, 'group2_Un_avg');
        createfigure_all_Un(taxis_second, Un_Avg_all, policy, mode, filepath, 'all_Un_avg');
        createfigure_all_Un(taxis_second, Un_Avg_max_group1, policy, mode, filepath, 'group1_Un_max');
        createfigure_all_Un(taxis_second, Un_Avg_max_group2, policy, mode, filepath, 'group2_Un_max');
        createfigure_all_Un(taxis_second, Un_Avg_max, policy, mode, filepath, 'all_Un_max');         
end

end

toc;