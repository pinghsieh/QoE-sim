function InsRate = get_channel_rate(N, cdf_fade, channel_rate_vec) 
    InsRate = zeros(N, 1);
    for i=1:N
        temp = rand(1); 
        comp = cdf_fade{i} >= temp; 
        [val, id] = max(comp);
        InsRate(i) = channel_rate_vec{i}(id);
    end
end