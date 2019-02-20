%% Configuration
% Unreliable channel, asymmetric case
% 10 links
% tau_n = tau = 200
% 
%% Common settings
N = 10;
Run = 1;
Ttot = 360000;
Taxis = 1:1:Ttot;
Nslots_per_second = 6000;
N_trials = 1;

%% Channel and playback
channel_rate_vec = cell(N, 1);
pn = zeros(N, 1);
cdf_pn = cell(N, 1);
dn{1} = zeros(N, 1);
frame_interval = zeros(N, 1);
lambda_n_max = zeros(N, 1);
beta_n_WLD = zeros(N, 1);
Tn_WRR = zeros(N, 1);
for k=1:5
    frame_interval(k) = 50;    
    channel_rate_vec{k} = [0, 1];
    pn(k) = 0.8;  % Channel reliability
    cdf_pn{k} = [0.2 1];
    lambda_n_max(k) = 0.08;
    dn{1}(k) = 100; % unit: slot
    beta_n_WLD(k) = 2;
    Tn_WRR(k) = 1;
end
for k=6:10
    frame_interval(k) = 100;
    channel_rate_vec{k} = [0, 1];
    pn(k) = 0.6;  % Channel reliability
    cdf_pn{k} = [0.4 1];
    lambda_n_max(k) = 0.06;
    dn{1}(k) = 200; % unit: slot
    beta_n_WLD(k) = 3;
    Tn_WRR(k) = 1;
end
video_type = "Live-with-drop";

% offset time for live streaming
rho = 1;
lambda_n = rho*lambda_n_max;    % Consumption rate
%dn{1} = [1000; 1000];          % unit: slot
%dn{2} = [200; 100];
%dn{3} = [50; 50];
%dn{4} = [30; 30];
%dn{5} = [10; 10];
Boffset = cell(N_trials, 1);
for j=1:N_trials
    Boffset{j} = dn{j}.*lambda_n;
end
       



