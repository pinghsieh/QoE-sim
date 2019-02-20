function Schedule = EDF(Qn, lambda_n, frame_interval)
    N = length(Qn);
    Schedule = zeros(N, 1);
    head_of_queue_deadline = frame_interval.*(ceil(Qn./(lambda_n.*frame_interval)));
    [val, id] = max(head_of_queue_deadline);
    Schedule(id) = 1;
end