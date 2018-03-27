%% Configuration
% ON/OFF channel
% 2 links + 1 fake link
%% Common settings
N = 3;
Run = 10;
Ttot = 500000;
Taxis = 1:1:Ttot;
Roots = [2 10];
frame_interval = 1;
N_trials = 1;
dummy = 1;

%% Channel and playback
channel_rate_vec{1} = [0, 1];
channel_rate_vec{2} = [0, 1];

% Channel probability for independent case
cdf_fade{1} = [0.5 1];
cdf_fade{2} = [0.5 1];

%video_type = "On-demand";
%video_type = "Live";
video_type = "Live-with-drop";

% offset time for live streaming
rho = 0.8;
qn_max = [3/8; 3/8];
qn = [rho*qn_max; (1-rho)*3/4];    % Consumption rate
Toffset{1} = [10; 10; 1e+10];
%Toffset{2} = [70; 70];
%Toffset{3} = [50; 50];
%Toffset{4} = [30; 30];
%Toffset{5} = [10; 10];
Boffset = cell(N_trials, 1);
for j=1:N_trials
    Boffset{j} = Toffset{j}.*qn;
end

%% For policy
policy = "HDR";
%policy = "HDRbuff";
%policy = "MW";
% Weight for HDR 
wn_HDR = [1; 1; 0];
% Weight for MW
wn_MW = [1; 1];           



