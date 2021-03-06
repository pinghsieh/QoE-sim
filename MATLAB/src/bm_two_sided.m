%% Simulate the CDF of D(t)
clear;
tic;

%% 1. Initialization
nRun = 10000;
T = 100000;
step = 1;
upper = 25;
lower = 0;
Dfinal = zeros(1,nRun);
Ufinal = zeros(1,nRun);
D = zeros(1,T);
U = zeros(1,T);
dX = zeros(1,T);

for i=1:nRun
%% 2. Random numbers
D = zeros(1,T);
U = zeros(1,T);
%X = upper/2;
X = 0;
dX = (2*randi([0 1],1,T) - 1)*step;

%% 3. Derive D(t) and U(t)
for t=2:T
    if X + dX(t) > upper
        U(t) = U(t-1) + 1;
        D(t) = D(t-1);
        X = min(X + dX(t), upper);
    else
        if X + dX(t) < lower
            D(t) = D(t-1) + 1;
            U(t) = U(t-1);
            X = max(X + dX(t), lower);           
        else
            X = X + dX(t);
            U(t) = U(t-1);
            D(t) = D(t-1);
        end
    end
end

Dfinal(i) = D(T);
Ufinal(i) = U(T);

end

%% 4. Plot CDF
taxis = 1:1:T;
%figure;
%plot(taxis, D, '-^r');
%figure;
%plot(taxis, U, '-ob');
figure;
h1 = cdfplot(Dfinal);
figure;
h2 = cdfplot(Ufinal);
set(h1, 'LineStyle', '-', 'Color', 'r');
set(h2, 'LineStyle', '-', 'Color', 'b');

toc;