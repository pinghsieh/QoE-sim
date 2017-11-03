function [expval_all, expval_subset, expval_subset_QoS, hist_, variance_] = admctrl2(Pr_, rate_, groupA, groupB)
% Probability matrix
%    r1    r2    r3    ...    rn  
%1  0.2   0.3  0.1  ...   
%2  0.1   0.5  0.1   ...
%3  0.3   0.3   0.1 ...
%
% Usage: 
% Pr_       :  describes the distribution of rates
% rate_    :  displays the possible rates in an descending order 
% groupA :  specifies the target group of clients
% groupB :  specifies the rest of the clients
[rdim, cdim] = size(Pr_);
sizeA = length(groupA);
sizeB = length(groupB);
% Check if all elements are nonnegative
for i=1:rdim
    for j=1:cdim
        if Pr_(i,j) < 0
            disp('All elements of the probability matrix should be nonnegative!!')
            return
        end
    end
end
%{
% Check row sum
for i=1:rdim
    rsum = 0;
    for j=1:cdim
        rsum = rsum + Pr_(i,j);
    end
    if ~(rsum == 1)
        disp('Row sum not equal to 1 ! Not a probability matrix !')
        return
    end
end
%}

hist_ = zeros(rdim,cdim);
% Get the histogram for prabability
for i=1:rdim
    for j=1:cdim
        if j == 1
            hist_(i,j) = 1 - Pr_(i,j);
        else
            hist_(i,j) = hist_(i,j-1) - Pr_(i,j);
        end
    end
end
expval_all = 0;
moment_2nd = 0;
% For HDR policy:
for i=1:cdim
    prob_p = 1;
    prob_n = 1;
    if i == 1
        for j=1:rdim         
            prob_n = prob_n*hist_(j,i); 
        end
        expval_all = expval_all + (1-prob_n)*rate_(1,i);
        moment_2nd = moment_2nd + (1-prob_n)*(rate_(1,i)^2);  
    else
        for j=1:rdim
            prob_n = prob_n*hist_(j,i); 
        end
        for k=1:rdim
            prob_p = prob_p*hist_(k,i-1); 
        end
        expval_all = expval_all + (prob_p-prob_n)*rate_(1,i);     
        moment_2nd = moment_2nd + (prob_p-prob_n)*(rate_(1,i)^2);  
    end
end
variance_ = moment_2nd - (expval_all^2);

% For HDR policy:
expval_subset = 0;
for i=1:cdim
    prob_p = 1;
    prob_n = 1;
    prob_o = 1;
    if i == 1
        for j=1:sizeA   
            prob_n = prob_n*hist_(groupA(1,j),i); 
        end
        expval_subset = expval_subset + (1-prob_n)*rate_(1,i);
    else
        for j=1:sizeA
            prob_n = prob_n*hist_(groupA(1,j),i); 
        end
        for k=1:sizeA
            prob_p = prob_p*hist_(groupA(1,k),i-1); 
        end
        for l=1:sizeB
            prob_o = prob_o*hist_(groupB(1,l),i-1);
        end
        expval_subset = expval_subset + prob_o*(prob_p-prob_n)*rate_(1,i);     
    end
end

% For QoS
expval_subset_QoS = 0;
for i=1:cdim
    prob_p = 1;
    prob_n = 1;
    %prob_o = 1;
    if i == 1
        for j=1:sizeA   
            prob_n = prob_n*hist_(groupA(1,j),i); 
        end
        expval_subset_QoS = expval_subset_QoS + (1-prob_n)*rate_(1,i);
    else
        for j=1:sizeA
            prob_n = prob_n*hist_(groupA(1,j),i); 
        end
        for k=1:sizeA
            prob_p = prob_p*hist_(groupA(1,k),i-1); 
        end
        %for l=1:sizeB
        %    prob_o = prob_o*hist_(groupB(1,l),i-1);
        %end
        expval_subset_QoS = expval_subset_QoS + (prob_p-prob_n)*rate_(1,i);     
    end
end


end