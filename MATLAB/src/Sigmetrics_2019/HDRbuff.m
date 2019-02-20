function Schedule = HDRbuff(Xn, wn, InsRate, Bn, Bn_max)
% Highest data rate policy with buffer awareness
% Schedule a client with the highest instantaneous rate, and break ties by
% using wn*Xn

    MIN_TOL = 1e-6;
    N = length(Xn);
    Schedule = zeros(N, 1);
    [val, id] = max(InsRate.*(Bn < Bn_max));
    if val > MIN_TOL
        id_max_rate = find((InsRate.*(Bn < Bn_max) == val));
        [val_Xn, nid] = min(wn(id_max_rate).*Xn(id_max_rate));
        Schedule(id_max_rate(nid)) = 1;
    end 
    
end