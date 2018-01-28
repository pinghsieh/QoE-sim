%% QoE Main Program
clear;
tic;

%% Part 1: Configuration
config.config_ONOFF_10links;


%% Part 2: Initialization
% For live streaming
Dn_Avg = zeros(N, Ttot);

%% Part 3: Main Operations
for i=1:Run
    %[D1, D2, D3, D4, D5, Z_MW, Z_PF, Z_HDR, Z_NOVA, Z_MW_r, W_NOVA, X_PF, DALL_MW, InsRate, TP1, TP2, TP3, TP4, TP5] = qoefade(N, Ttot, Pfade, qnN, wn, Roots, alpha, beta);
    InsRate_history = zeros(N,Ttot);
    AvgRate_history = zeros(N,Ttot);  % AvgRate(i,t) denotes average throughput up to time t
    Schedule_history = zeros(N,Ttot);
    Xn_history = zeros(N,Ttot);
    Dn_history = zeros(N,Ttot);
    Bn_history = zeros(N,Ttot);
    Zn_history = zeros(N,Ttot);
    Xn = zeros(N,1);
    Dn = zeros(N,1);
    Bn = zeros(N,1);
    Zn = zeros(N,1);
    
    %% slot-wise update
    for t=1:Ttot
        % Get playback rates
        InsRate = get_channel_rate(N, cdf_fade, channel_rate_vec);
        
        % Scheduling: 
        % schedule is a N-by-1 boolean vector
        switch policy
            case "HDR"
                schedule = HDR(Zn, wn_HDR, InsRate);
                
            case "MW"
                schedule = MW(Xn, wn_MW, InsRate);
                
            case "WPF"
                
                
            case "MW-r"
            
                
            case "NOVA"
                
                
        end
        % Update state variables
        switch video_type
            case "On-demand"
                Xn_next = Xn + schedule.*InsRate; % for on-demand video streaming
                Bn = Bn + (Xn_next - Xn);
                Xn = Xn_next;  
                Zn = Zn + schedule.*InsRate;
            case "Live"
                Xn_next = min(Toffset*qn, Xn + schedule.*InsRate); % for live streaming
                Bn = Bn + (Xn_next - Xn);
                Xn = Xn_next;
                Zn = Zn + schedule.*InsRate;
            case "Live-with-drop"
                Xn_next = min(Toffset*qn, Xn + Dn.*qn + schedule.*InsRate) - Dn.*qn; % for live streaming with packet dropping             
                Bn = Bn + (Xn_next - Xn);
                Xn = Xn_next;
                Zn = Zn + schedule.*InsRate;
        end
        if t > 1
            AvgRate_history(:,t) = (AvgRate_history(:,t-1)*(t-1) + InsRate)/t;
        else
            AvgRate_history(:,t) = InsRate;
        end
        about_to_play = double(rem(t - Dn, frame_interval) == 0);
        buffer_empty = double(Bn < (qn*frame_interval));
        about_to_rebuffer = about_to_play.*buffer_empty;
        ready_for_playback = about_to_play - about_to_rebuffer;
        Dn = Dn + ones(N,1).*about_to_rebuffer;
        % Playback
        Bn = Bn - frame_interval*(qn.*ready_for_playback);
        Xn = Xn - qn;
        Zn = Zn - qn;
        % Update history
        Xn_history(:,t) = Xn;
        Dn_history(:,t) = Dn;
        Bn_history(:,t) = Bn;
        Zn_history(:,t) = Zn;
        Schedule_history(:,t) = schedule;
    end
    
    
    %% Update average history
    Dn_Avg = Dn_Avg + Dn_history;
    
end
Dn_Avg = Dn_Avg./Run;
DAll_Avg = sum(Dn_Avg);

%% Part 4: Plotting
step_size = 1;
tplot_start = 1;
taxis = tplot_start:step_size:Ttot;
createfigure_all(taxis, Dn_Avg(:,taxis));
createfigure_all(taxis, DAll_Avg(taxis));

toc;