%% Configuration
% ON/OFF channel
% 2 links
%% Common settings
N = 2;
Run = 1;
Ttot = 10000;
Taxis = 1:1:Ttot;
frame_interval = 8;
N_trials = 1;
dummy_client = 0;

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
rho = 4/3;
qn_max = [3/8; 3/8];
qn = rho*qn_max;    % Consumption rate
Toffset{1} = [40; 40];
%Toffset{1} = [200; 100];
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
wn_HDR = [1; 1];
% Weight for MW
wn_MW = [1; 1];           



