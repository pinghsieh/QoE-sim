function Schedule = MW(Xn, wn, InsRate)
% Highest data rate policy
% Schedule a client with the highest instantaneous rate, and break ties by
% using wn*Xn

    MIN_TOL = 1e-6;
    N = length(Xn);
    Schedule = zeros(N, 1);
    [val, id] = max(InsRate);
    if val > MIN_TOL
        [val_MW, nid] = min((wn.*Xn).*InsRate);
        Schedule(nid) = 1;
    end 
    
end