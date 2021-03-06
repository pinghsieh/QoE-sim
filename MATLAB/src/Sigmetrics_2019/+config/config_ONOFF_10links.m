%% Configuration
% ON/OFF channel
% 10 links
%% Common settings
N = 10;
Run = 5;
Ttot = 30000;
Taxis = 1:1:Ttot;
Roots = [2 10];
frame_interval = 5;
N_trials = 1;
dummy = 0;

rho = 1;
qn = rho*(1023/1024)/10*ones(N,1);    % Consumption rate

% offset time for live streaming
%Toffset{1} = [20; 20; 20; 20; 20; 40; 40; 40; 40; 40];
%Toffset{1} = [40; 40; 40; 40; 40; 20; 20; 20; 20; 20];
Toffset{1} = [10; 10; 10; 10; 10; 10; 10; 10; 10; 10];

Boffset = cell(N_trials, 1);
for j=1:N_trials
    Boffset{j} = Toffset{j}.*qn;
end

%% Channel and playback
channel_rate_vec{1} = [0, 1];
channel_rate_vec{2} = [0, 1];
channel_rate_vec{3} = [0, 1];
channel_rate_vec{4} = [0, 1];
channel_rate_vec{5} = [0, 1];
channel_rate_vec{6} = [0, 1];
channel_rate_vec{7} = [0, 1];
channel_rate_vec{8} = [0, 1];
channel_rate_vec{9} = [0, 1];
channel_rate_vec{10} = [0, 1];

% Channel probability for independent case
cdf_fade{1} = [0.5 1];
cdf_fade{2} = [0.5 1];
cdf_fade{3} = [0.5 1];
cdf_fade{4} = [0.5 1];
cdf_fade{5} = [0.5 1];
cdf_fade{6} = [0.5 1];
cdf_fade{7} = [0.5 1];
cdf_fade{8} = [0.5 1];
cdf_fade{9} = [0.5 1];
cdf_fade{10} = [0.5 1];

%video_type = "On-demand";
%video_type = "Live";
video_type = "Live-with-drop";


%% For policy
policy = "HDR";
%policy = "MW";
% Weight for HDR
%wn_HDR = ones(N, 1); 
wn_HDR = [2; 2; 2; 2; 2; 1; 1; 1; 1; 1];
% Weight for MW
wn_MW = ones(N, 1);   


