%qoe main program
%#####Setting of Simple ON/OFF channel
%{
N = 3;
%Pin=[1/4 1/5 1/6];  %Channel probability for independent case
Pin=[1/2 1/2 1/2];  %Channel probability for independent case
%Pin=[1/4 1/3 1/2];  %Channel probability for independent case
TSetFlag=[0 3];
%qn=[0.25 0.14 0.11];    %Consumption rate
%cn=[0.25 0.14 0.11];    %Weight for PF
qn=[0.49 0.21 0.175];    %Consumption rate
cn=[0.49 0.21 0.175];    %Weight for PF
%qn=[0.25 0.2 0.3];    %Consumption rate
%cn=[0.25 0.2 0.3];    %Weight for PF
wn=[1 1 1];            %Weight for MW
Pmc(:,:,1)=[0.2 0.8; 0.2 0.8];  %Channel probability for Markov model
Pmc(:,:,2)=[0.5 0.5; 0.5 0.5];  %Channel probability for Markov model
%}

%#####Setting of Fading channel

N = 3;
%N = 2;
%Pfade=[0.5 0.5 0; 0 0.5 0.5];
%Pfade=[0.25 0.5 0.25; 0.5 0.25 0.25];
%Pfade=[0.5 0 0.5; 0.5 0 0.5];
%Pfade=[1/3 1/3 1/3; 1/3 1/3 1/3];  %Channel probability for independent case
Pfade=[1/3 1/3 1/3; 1/3 1/3 1/3; 1/3 1/3 1/3];  %Channel probability for independent case
%qnN=[16/32*2 5/32*2];    %Consumption rate for fading channels
%qnN=[0.9 0.6];    %Consumption rate for fading channels
%qnN=[3/8 9/8];    %Consumption rate
qnN=[15/27 15/27 15/27];    %Consumption rate
%qnN=[19/27 19/27 7/27];    %Consumption rate
%qnN=[23/27 11/27 11/27];    %Consumption rate
%wn=[1000 1];            %Weight for MW
wn=[1 1 1];            %Weight for MW
%alpha = [0.2 0.2 0.2];
%alpha = [0.145 0.145 0.145];
alpha = [0 0 0];
%beta = [0.725 0.725 0.725];
beta = [0 0 0];
%#####Common Settings
Run = 5;
Ttot = 30000;
Taxis = 0:1:Ttot;
Taxis2 = 1:1:Ttot;
Roots = [2 10];
sizeRoots = length(Roots);
D_MW_Avg = zeros(N, Ttot+1);
DALL_MW_Avg = zeros(1, Ttot+1);
D_PF_Avg = zeros(N, Ttot+1);
D_TMW_Avg = zeros(N, Ttot+1);
D_HDR_Avg = zeros(N, Ttot+1);
D_NOVA_Avg = zeros(N, Ttot+1);
D_MW_r_Avg = zeros(sizeRoots, N, Ttot+1);
%D_MW_Log_Avg = zeros(N, Ttot+1);
Y_PF_Avg = zeros(N, Ttot+1);
Y_MW_Avg = zeros(N, Ttot+1);
Y_TMW_Avg = zeros(N, Ttot+1);
Z_MW_Avg = zeros(1, Ttot+1);
Z_PF_Avg = zeros(1, Ttot+1);
Z_HDR_Avg = zeros(1, Ttot+1);
Z_NOVA_Avg = zeros(1, Ttot+1);
Z_MW_r_Avg = zeros(sizeRoots, 1, Ttot+1);
%Z_MW_Log_Avg = zeros(1, Ttot+1);
W_NOVA_Avg = zeros(N, Ttot+1);
X_PF_Avg = zeros(N, Ttot+1);

TP_MW_Avg = zeros(N, Ttot);
TP_PF_Avg = zeros(N, Ttot);
TP_TMW_Avg = zeros(N, Ttot);
TP_HDR_Avg = zeros(N, Ttot);
TP_NOVA_Avg = zeros(N, Ttot);
TP_MW_r_Avg = zeros(sizeRoots, N, Ttot);

for i=1:Run
    [D1, D2, D3, D4, D5, Z_MW, Z_PF, Z_HDR, Z_NOVA, Z_MW_r, W_NOVA, X_PF, DALL_MW, InsRate, TP1, TP2, TP3, TP4, TP5] = qoefade(N, Ttot, Pfade, qnN, wn, Roots, alpha, beta);
    D_MW_Avg = D_MW_Avg + D1;
    D_PF_Avg = D_PF_Avg + D2;
    D_HDR_Avg = D_HDR_Avg + D3;
    D_NOVA_Avg = D_NOVA_Avg + D4;
    for r=1:sizeRoots
        D_MW_r_Avg(r,:,:) = D_MW_r_Avg(r,:,:) + D5(r,:,:);
    end
    
    TP_MW_Avg = TP_MW_Avg + TP1;
    TP_PF_Avg = TP_PF_Avg + TP2;
    TP_HDR_Avg = TP_HDR_Avg + TP3;
    TP_NOVA_Avg = TP_NOVA_Avg + TP4;
    for r=1:sizeRoots
        TP_MW_r_Avg(r,:,:) = TP_MW_r_Avg(r,:,:) + TP5(r,:,:);
    end
    %D_MW_Log_Avg = D_MW_Log_Avg + D6;
    %D_TMW_Avg = D_TMW_Avg + D3;
    %Y_PF_Avg = Y_PF_Avg + Y_PF;
    %Y_MW_Avg = Y_MW_Avg + Y_MW;
    %Y_TMW_Avg = Y_TMW_Avg + Y_TMW;
    Z_MW_Avg = Z_MW_Avg + Z_MW;
    Z_PF_Avg = Z_PF_Avg + Z_PF;
    Z_HDR_Avg = Z_HDR_Avg + Z_HDR;
    Z_NOVA_Avg = Z_NOVA_Avg + Z_NOVA;
    for r=1:sizeRoots
        Z_MW_r_Avg(r,:,:) = Z_MW_r_Avg(r,:,:) + Z_MW_r(r,:,:);
    end
    %Z_MW_Log_Avg = Z_MW_Log_Avg + Z_MW_Log;
    W_NOVA_Avg = W_NOVA_Avg + W_NOVA;
    %X_PF_Avg = X_PF_Avg + X_PF;
    %DALL_MW_Avg = DALL_MW_Avg + DALL_MW;
end
D_MW_Avg = D_MW_Avg./Run;
D_PF_Avg = D_PF_Avg./Run;
D_HDR_Avg = D_HDR_Avg./Run;
D_NOVA_Avg = D_NOVA_Avg./Run;
D_MW_r_Avg = D_MW_r_Avg./Run;
%D_MW_Log_Avg = D_MW_Log_Avg./Run;
%D_TMW_Avg = D_TMW_Avg./Run;
%Y_PF_Avg = Y_PF_Avg./Run;
%Y_MW_Avg = Y_MW_Avg./Run;
%Y_TMW_Avg = Y_TMW_Avg./Run;
Z_MW_Avg = Z_MW_Avg./Run;
Z_PF_Avg = Z_PF_Avg./Run;
Z_HDR_Avg = Z_HDR_Avg./Run;
Z_NOVA_Avg = Z_NOVA_Avg./Run;
Z_MW_r_Avg = Z_MW_r_Avg./Run;
%Z_MW_Log_Avg = Z_MW_Log_Avg./Run;
W_NOVA_Avg = W_NOVA_Avg./Run;
%DALL_MW_Avg = DALL_MW_Avg./Run;
TP_MW_Avg = TP_MW_Avg./Run;
TP_PF_Avg = TP_PF_Avg./Run;
TP_HDR_Avg = TP_HDR_Avg./Run;
TP_NOVA_Avg = TP_NOVA_Avg./Run;
TP_MW_r_Avg = TP_MW_r_Avg./Run;
%#####Plot Setting for simple ON/OFF channel

figure;
plot(Taxis,D_MW_Avg(1,:),'-*r');
hold on;
plot(Taxis,D_PF_Avg(1,:),'-og');
plot(Taxis,D_HDR_Avg(1,:),'-ob');
plot(Taxis,D_NOVA_Avg(1,:),'-om');
%plot(Taxis,D_MW_Log_Avg(1,:),'-or');
plot(Taxis,squeeze(D_MW_r_Avg(1,1,:)),'-*k');
plot(Taxis,squeeze(D_MW_r_Avg(2,1,:)),'-*c');
%plot(Taxis,squeeze(D_MW_r_Avg(3,1,:)),'-*y');
%plot(Taxis,squeeze(D_MW_r_Avg(4,1,:)),'-*g');
%plot(Taxis,squeeze(D_MW_r_Avg(5,1,:)),'-*m');
%plot(Taxis,D_TMW_Avg(1,:),'-ob');
figure;
plot(Taxis,D_MW_Avg(2,:),'-*r');
hold on;
plot(Taxis,D_PF_Avg(2,:),'-og');
plot(Taxis,D_HDR_Avg(2,:),'-ob');
plot(Taxis,D_NOVA_Avg(2,:),'-om');
%plot(Taxis,D_MW_Log_Avg(2,:),'-or');
plot(Taxis,squeeze(D_MW_r_Avg(1,2,:)),'-*k');
plot(Taxis,squeeze(D_MW_r_Avg(2,2,:)),'-*c');
%plot(Taxis,squeeze(D_MW_r_Avg(3,2,:)),'-*y');
%plot(Taxis,squeeze(D_MW_r_Avg(4,2,:)),'-*g');
%plot(Taxis,squeeze(D_MW_r_Avg(5,2,:)),'-*m');
%plot(Taxis,D_TMW_Avg(2,:),'-ob');


figure;
plot(Taxis,D_MW_Avg(3,:),'-*r');
hold on;
plot(Taxis,D_PF_Avg(3,:),'-og');
plot(Taxis,D_HDR_Avg(3,:),'-ob');
plot(Taxis,D_NOVA_Avg(3,:),'-om');
%plot(Taxis,D_MW_Log_Avg(3,:),'-or');
plot(Taxis,squeeze(D_MW_r_Avg(1,3,:)),'-*k');
plot(Taxis,squeeze(D_MW_r_Avg(2,3,:)),'-*c');
%plot(Taxis,squeeze(D_MW_r_Avg(3,3,:)),'-*y');
%plot(Taxis,squeeze(D_MW_r_Avg(4,3,:)),'-*g');
%plot(Taxis,squeeze(D_MW_r_Avg(5,3,:)),'-*m');
%plot(Taxis,D_TMW_Avg(3,:),'-ob');


figure;
plot(Taxis,Z_MW_Avg(1,:),'-*r');
hold on;
plot(Taxis,Z_PF_Avg(1,:),'-*g');
plot(Taxis,Z_HDR_Avg(1,:),'-*b');
plot(Taxis,Z_NOVA_Avg(1,:),'-*m');
%plot(Taxis,Z_MW_Log_Avg(1,:),'-*r');
plot(Taxis,squeeze(Z_MW_r_Avg(1,1,:)),'-*k');
plot(Taxis,squeeze(Z_MW_r_Avg(2,1,:)),'-*c');
%plot(Taxis,squeeze(Z_MW_r_Avg(3,1,:)),'-*y');
%plot(Taxis,squeeze(Z_MW_r_Avg(4,1,:)),'-*g');
%plot(Taxis,squeeze(Z_MW_r_Avg(5,1,:)),'-*m');

figure;
plot(Taxis2,TP_MW_Avg(1,:),'-*r');
hold on;
plot(Taxis2,TP_PF_Avg(1,:),'-*g');
plot(Taxis2,TP_HDR_Avg(1,:),'-*b');
plot(Taxis2,TP_NOVA_Avg(1,:),'-*m');
%plot(Taxis,Z_MW_Log_Avg(1,:),'-*r');
plot(Taxis2,squeeze(TP_MW_r_Avg(1,1,:)),'-*k');
plot(Taxis2,squeeze(TP_MW_r_Avg(2,1,:)),'-*c');

%{
figure;
plot(Taxis,W_NOVA_Avg(1,:),'-*r');
hold on;
plot(Taxis,W_NOVA_Avg(2,:),'-*g');
plot(Taxis,W_NOVA_Avg(3,:),'-*b');
%}
%figure;
%plot(Taxis,DALL_MW_Avg(1,:),'-*b');
%}
%{
figure;
plot(Taxis,X_PF_Avg(1,:),'-*r');
hold on;
plot(Taxis,X_PF_Avg(2,:),'-*g');
plot(Taxis,X_PF_Avg(3,:),'-*b');
%}
%{
figure;
plot(Taxis,Y_PF_Avg(1,:),'-og');
hold on;
plot(Taxis,Y_MW_Avg(1,:),'-or');
plot(Taxis,Y_TMW_Avg(1,:),'-ob');
%}
%#####Plot Setting for Fading ON/OFF channel
%{
figure;
plot(Taxis,D_MW_Avg(1,:),'-*r');
hold on;
plot(Taxis,D_PF_Avg(1,:),'-og');
figure;
plot(Taxis,D_MW_Avg(2,:),'-*r');
hold on;
plot(Taxis,D_PF_Avg(2,:),'-og');
figure;
plot(Taxis,D_MW_Avg(1,:),'-*r');
hold on;
plot(Taxis,D_MW_Avg(2,:),'-*b');
%}
%{
figure;
plot(Taxis,D_MW_Avg(3,:),'-*m');
hold on;
plot(Taxis,D_PF_Avg(3,:),'-ok');
%}