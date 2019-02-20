function Schedule = LDF(Xn)
% Largest Deficit First policy
% Schedule a client with the smallest Xn
    MIN_TOL = 1e-10;
    N = length(Xn);
    Schedule = zeros(N, 1);
    [val, id] = min(Xn);
    Schedule(id) = 1;
end