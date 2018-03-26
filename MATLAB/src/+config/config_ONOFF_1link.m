%% Configuration
% ON/OFF channel
% 1 link
%% Common settings
N = 1;
Run = 10;
Ttot = 50000;
Taxis = 1:1:Ttot;
Roots = [2 10];
frame_interval = 1;

%% Channel and playback
channel_rate_vec{1} = [0, 1];

% Channel probability for independent case
cdf_fade{1} = [0.5 1];

%video_type = "On-demand";
%video_type = "Live";
video_type = "Live-with-drop";

% offset time for live streaming
rho = 1;
qn_max = 1/2;
qn = rho*qn_max;    % Consumption rate
Toffset = 10;
Boffset = Toffset.*qn;

%% For policy
policy = "HDR";
%policy = "HDRbuff";
%policy = "MW";
% Weight for HDR 
wn_HDR = 1;
% Weight for MW
wn_MW = 1;           



