function Schedule = WRand(lambda_n, pn)
% Weighted Round-Robin policy
% Schedule a client with the smallest Zn/beta_n
    N = length(lambda_n);
    Schedule = zeros(N, 1);
    workload = lambda_n./pn;
    randP = workload/sum(workload);
    randP_cumulated = cumsum(randP);
    index = find(rand(1) <= randP_cumulated);
    Schedule(index(1)) = 1;
end