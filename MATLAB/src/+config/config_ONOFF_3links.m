%% Configuration
% ON/OFF channel
% 3 links
%% Common settings
N = 3;
Run = 50;
Ttot = 10000;
Taxis = 1:1:Ttot;
Roots = [2 10];
frame_interval = 40;

% offset time for live streaming
Toffset = 40;

%% Channel and playback
channel_rate_vec{1} = [0, 1];
channel_rate_vec{2} = [0, 1];
channel_rate_vec{3} = [0, 1];

% Channel probability for independent case
cdf_fade{1} = [0.5 1];
cdf_fade{2} = [0.5 1];
cdf_fade{3} = [0.5 1];

video_type = "On-demand";

rho = 1;
qn = rho*[7/24; 7/24; 7/24];    % Consumption rate

%% For policy
policy = "HDR";
% Weight for HDR
wn_HDR = [1; 1; 1];            



