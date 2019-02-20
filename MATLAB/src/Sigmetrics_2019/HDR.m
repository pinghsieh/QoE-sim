function Schedule = HDR(Xn, wn, InsRate)
% Highest data rate policy
% Schedule a client with the highest instantaneous rate, and break ties by
% using wn*Xn

    MIN_TOL = 1e-6;
    N = length(Xn);
    Schedule = zeros(N, 1);
    [val, id] = max(InsRate > MIN_TOL);
    if val > MIN_TOL
        id_max_rate = find(InsRate == val);
        %[val_Xn, nid] = min(wn(id_max_rate).*Xn(id_max_rate));
        nids = find(wn(id_max_rate).*Xn(id_max_rate) == min(wn(id_max_rate).*Xn(id_max_rate)));
        Schedule(id_max_rate(datasample(nids,1))) = 1;
        
    end 
    
end