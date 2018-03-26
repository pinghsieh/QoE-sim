%% Configuration
% ON/OFF channel
% 3 links
%% Common settings
N = 3;
Run = 5;
Ttot = 100000;
Taxis = 1:1:Ttot;
Roots = [2 10];
frame_interval = 1;

% offset time for live streaming
Toffset = [20; 20; 40];

%% Channel and playback
channel_rate_vec{1} = [0, 1];
channel_rate_vec{2} = [0, 1];
channel_rate_vec{3} = [0, 1];

% Channel probability for independent case
cdf_fade{1} = [0.5 1];
cdf_fade{2} = [0.5 1];
cdf_fade{3} = [0.5 1];

%video_type = "On-demand";
%video_type = "Live";
video_type = "Live-with-drop";

rho = 1;
qn = rho*[7/24; 7/24; 7/24];    % Consumption rate

%% For policy
policy = "HDR";
%policy = "HDRbuff";
%policy = "MW";
% Weight for HDR
%wn_HDR = ones(N, 1);  
wn_HDR = [2; 2; 1];
% Weight for MW
wn_MW = ones(N, 1);           



