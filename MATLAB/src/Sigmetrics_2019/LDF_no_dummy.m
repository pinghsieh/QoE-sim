function Schedule = LDF_no_dummy(Xn, Qn)
% Largest Deficit First policy
% Schedule a client with the smallest Xn
    MIN_TOL = 1e-10;
    N = length(Xn);
    index_Qn_nonempty = find((Qn > MIN_TOL));
    Schedule = zeros(N, 1);
    if isempty(index_Qn_nonempty) == 0
        [val, id] = min(Xn(index_Qn_nonempty));
        sid = index_Qn_nonempty(id);
    else
        sid = randi(N);
    end
    Schedule(sid) = 1;
end