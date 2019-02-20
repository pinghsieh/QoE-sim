function Schedule = WLD(Zn, beta_n)
% Weighted Largest Deficit policy
% Schedule a client with the smallest Zn/beta_n
    rand_scale = 1e-2;
    N = length(Zn);
    Schedule = zeros(N, 1);
    [val, id] = min((Zn./(beta_n + rand_scale*rand(N,1))));
    Schedule(id) = 1;
end