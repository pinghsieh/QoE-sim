function Schedule = WRR(Tn,currentT)
% Weighted Round-Robin policy
% Schedule a client with the smallest Zn/beta_n
    N = length(Tn);
    Schedule = zeros(N, 1);
    Tn_cumulated = cumsum(Tn) - 1;
    currentT_modulus = mod(currentT-1, sum(Tn));
    index = find(currentT_modulus <= Tn_cumulated);
    Schedule(index(1)) = 1;
end