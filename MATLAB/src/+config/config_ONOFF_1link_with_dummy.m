%% Configuration
% ON/OFF channel
% 1 link + 1 fake link
%% Common settings
N = 2;
Run = 10;
Ttot = 500000;
Taxis = 1:1:Ttot;
Roots = [2 10];
frame_interval = 1;
N_trials = 1;
dummy = 1;

%% Channel and playback
channel_rate_vec{1} = [0, 1];

% Channel probability for independent case
cdf_fade{1} = [0.5 1];

%video_type = "On-demand";
%video_type = "Live";
video_type = "Live-with-drop";

% offset time for live streaming
rho = 0.8;
qn_max = 1/2;
qn = [rho*qn_max; (1-rho)*1/2]; % Consumption rate
Toffset{1} = [10; 1e+10];
Boffset = cell(N_trials, 1);
for j=1:N_trials
    Boffset{j} = Toffset{j}.*qn;
end
%% For policy
policy = "HDR";
%policy = "HDRbuff";
%policy = "MW";
% Weight for HDR 
wn_HDR = [1; 0];
% Weight for MW
wn_MW = 1;           



