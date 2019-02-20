%% Configuration
% Unreliable channel, symmetric case
% 2 links
% tau_n = tau = 200
% 
%% Common settings
N = 2;
Run = 10;
Ttot = 300000;
Taxis = 1:1:Ttot;
frame_interval = 200;
Nslots_per_second = 6000;
N_trials = 1;

%% Channel and playback
channel_rate_vec{1} = [0, 1];
channel_rate_vec{2} = [0, 1];
pn = [0.8; 0.6]; % Channel reliability
cdf_pn{1} = [0.2 1];
cdf_pn{2} = [0.4 1];

video_type = "Live-with-drop";

% offset time for live streaming
rho = 1;
lambda_n_max = [0.6; 0.15];
lambda_n = rho*lambda_n_max;    % Consumption rate
dn{1} = [1000; 1000];          % unit: slot
%dn{2} = [200; 100];
%dn{3} = [50; 50];
%dn{4} = [30; 30];
%dn{5} = [10; 10];
Boffset = cell(N_trials, 1);
for j=1:N_trials
    Boffset{j} = dn{j}.*lambda_n;
end

%% For policy
%policy = "WLD";
%policy = "WRR";
%policy = "WRand";
policy = "LDF";
%policy = "EDF";

% Weights for WLD 
beta_n_WLD = [4; 1];
% Parameters for WRR
Tn_WRR = [3; 1];           



