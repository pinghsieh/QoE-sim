function [D_MW, D_PF, D_HDR, D_NOVA, D_MW_r, Z_MW, Z_PF, Z_HDR, Z_NOVA, Z_MW_r, W_NOVA, X_PF, DALL_MW, InsRate, AvgRate_MW, AvgRate_PF, AvgRate_HDR, AvgRate_NOVA, AvgRate_MW_r] = qoefade10(N, Ttot, Pch, qn, wn, Roots, alpha, beta, rate_, offset_)
%**********QoE-driven simulation**********
%N=2;                %Number of clients
%Ttot=200000;           %Total simulation time
%Pch=[0.5 0.6];  %Channel probability for independent case
%qn=[0.4 0.3];    %Consumption rate
%Taxis=1:1:Ttot;
%Taxis2=0:1:Ttot;
InsRate=zeros(N,Ttot);
[rdim_, cdim_] = size(Pch);
%*****Max-Weight Policy
AvgRate_MW=zeros(N,Ttot);  %AvgRate(i,t) denotes average throughput up to time t
Schedule_MW=zeros(N,Ttot);
Receive_MW=zeros(N,Ttot);
VQueue_MW=zeros(N,Ttot+1);
X_MW=zeros(N,Ttot+1);
D_MW=zeros(N,Ttot+1);
Z_MW=zeros(1,Ttot+1);
XALL_MW=zeros(1,Ttot+1);
DALL_MW=zeros(1,Ttot+1);
sizeRoots = length(Roots);
%***Determine channel states
for t=1:Ttot
    for i=1:N
        y=rand;
        cumul_ = 0;
        flag_ = 0;
        for k=1:cdim_
            cumul_ = cumul_ + Pch(i,k);
            if y<=cumul_ && flag_ == 0
                InsRate(i,t) = rate_(1,k);        
                flag_ = 1;
            end
        end
    end
end

%MW: Slot-wise update
for t=1:Ttot
    head = 0;   %Pick the one with ON channel and smallest Xn(t)
    Dbmax = -1e10;   %Track minimum Xn(t) in time slot t
    Ratemax = 0;
    best = 0;
    for i=1:N
        if InsRate(i,t) > Ratemax
            best = i;
            Ratemax = InsRate(i,t);
        end
    end
    
    for i=1:N
        X_MW(i,t+1)=X_MW(i,t)-qn(1,i);
        if t == 1
            if InsRate(i,t) > 0 && InsRate(i,t) > Dbmax
                head = i;
                Dbmax = InsRate(i,t);
            end    
        else
            VQueue_MW(i,t) = InsRate(i,t)*(subplus(-X_MW(i,t)))*wn(1,i);
%            VQueue_MW(i,t) = InsRate(i,t)*(-X_MW(i,t))*wn(1,i);
            if (InsRate(i,t) > 0) && (VQueue_MW(i,t) > Dbmax)
                Dbmax = VQueue_MW(i,t);
                head = i;
            end
        end
    end
%    z = rand;
    if head > 0 
        Schedule_MW(head,t)=1;
        %if z <= InsRate(head,t)
        X_MW(head,t+1)=X_MW(head,t+1)+InsRate(head,t);
        Receive_MW(head,t)=InsRate(head,t);
        %end
    end
    
    if head > 0 && best > 0 && InsRate(head,t) < InsRate(best,t)
        Z_MW(1,t+1) = Z_MW(1,t) + 1;
    else
        Z_MW(1,t+1) = Z_MW(1,t);
    end
    
    for i=1:N
        if i == head && Receive_MW(i,t) > 0
            if t == 1
                AvgRate_MW(i,t)=Receive_MW(i,t);
            else
                %AvgRate_MW(i,t)=AvgRate_MW(i,t-1)/t*(t-1)+1/t;
                AvgRate_MW(i,t)=AvgRate_MW(i,t-1)+Receive_MW(i,t);
            end
        else
            if t > 1
                %AvgRate_MW(i,t)=AvgRate_MW(i,t-1)/t*(t-1);
                AvgRate_MW(i,t)=AvgRate_MW(i,t-1);
            end
        end
        if X_MW(i,t+1)<0 && floor(-X_MW(i,t+1)/qn(1,i))>D_MW(i,t)
            D_MW(i,t+1) = floor(-X_MW(i,t+1)/qn(1,i));
        else
            D_MW(i,t+1) = D_MW(i,t);
        end            
    end
%Find X(t) and D(t)
    for i=1:N
        XALL_MW(1,t+1) = XALL_MW(1,t+1) + X_MW(i,t+1);
    end
    if XALL_MW(1,t+1)<0 && floor(-XALL_MW(1,t+1))>DALL_MW(1,t)
        DALL_MW(1,t+1) = floor(-XALL_MW(1,t+1));
    else
        DALL_MW(1,t+1) = DALL_MW(1,t);
    end     
end

%*****General Max-Weight Policy to the rth root
AvgRate_MW_r=zeros(sizeRoots, N,Ttot);  %AvgRate(i,t) denotes average throughput up to time t
Schedule_MW_r=zeros(sizeRoots, N,Ttot);
Receive_MW_r=zeros(sizeRoots, N,Ttot);
VQueue_MW_r=zeros(sizeRoots, N,Ttot+1);
X_MW_r=zeros(sizeRoots, N,Ttot+1);
D_MW_r=zeros(sizeRoots, N,Ttot+1);
Z_MW_r=zeros(sizeRoots, 1,Ttot+1);
XALL_MW_r=zeros(sizeRoots, 1,Ttot+1);
DALL_MW_r=zeros(sizeRoots, 1,Ttot+1);

%MW_r: Slot-wise update
for r=1:sizeRoots
for t=1:Ttot
    head = 0;   %Pick the one with ON channel and smallest Xn(t)
    Dbmax = -1e10;   %Track minimum Xn(t) in time slot t
    Ratemax = 0;
    best = 0;
    for i=1:N
        if InsRate(i,t) > Ratemax
            best = i;
            Ratemax = InsRate(i,t);
        end
    end
    
    for i=1:N
        X_MW_r(r,i,t+1)=X_MW_r(r,i,t)-qn(1,i);
        if t == 1
            if InsRate(i,t) > 0 && InsRate(i,t) > Dbmax
                head = i;
                Dbmax = InsRate(i,t);
            end    
        else
            VQueue_MW_r(r,i,t) = InsRate(i,t)*nthroot(subplus((-X_MW_r(r,i,t))*wn(1,i)),Roots(1,r));
%            VQueue_MW_r(r,i,t) = InsRate(i,t)*nthroot((-X_MW_r(r,i,t))*wn(1,i),Roots(1,r));
            if (InsRate(i,t) > 0) && (VQueue_MW_r(r,i,t) > Dbmax)
                Dbmax = VQueue_MW_r(r,i,t);
                head = i;
            end
        end
    end
%    z = rand;
    if head > 0 
        Schedule_MW_r(r,head,t)=1;
        %if z <= InsRate(head,t)
        X_MW_r(r,head,t+1)=X_MW_r(r,head,t+1)+InsRate(head,t);
        Receive_MW_r(r,head,t)=InsRate(head,t);
        %end
    end
    
    if head > 0 && best > 0 && InsRate(head,t) < InsRate(best,t)
        Z_MW_r(r,1,t+1) = Z_MW_r(r,1,t) + 1;
    else
        Z_MW_r(r,1,t+1) = Z_MW_r(r,1,t);
    end
    
    for i=1:N
        if i == head && Receive_MW_r(r,i,t) > 0
            if t == 1
                AvgRate_MW_r(r,i,t)=Receive_MW_r(r,i,t);
            else
                %AvgRate_MW_r(i,t)=AvgRate_MW_r(i,t-1)/t*(t-1)+1/t;
                AvgRate_MW_r(r,i,t)=AvgRate_MW_r(r,i,t-1)+Receive_MW_r(r,i,t);
            end
        else
            if t > 1
                %AvgRate_MW_r(i,t)=AvgRate_MW_r(i,t-1)/t*(t-1);
                AvgRate_MW_r(r,i,t)=AvgRate_MW_r(r,i,t-1);
            end
        end
        if X_MW_r(r,i,t+1)<0 && floor(-X_MW_r(r,i,t+1)/qn(1,i))>D_MW_r(r,i,t)
            D_MW_r(r,i,t+1) = floor(-X_MW_r(r,i,t+1)/qn(1,i));
        else
            D_MW_r(r,i,t+1) = D_MW_r(r,i,t);
        end            
    end
%Find X(t) and D(t)
    for i=1:N
        XALL_MW_r(r,1,t+1) = XALL_MW_r(r,1,t+1) + X_MW_r(r,i,t+1);
    end
    if XALL_MW_r(r,1,t+1)<0 && floor(-XALL_MW_r(r,1,t+1))>DALL_MW_r(r,1,t)
        DALL_MW_r(r,1,t+1) = floor(-XALL_MW_r(r,1,t+1));
    else
        DALL_MW_r(r,1,t+1) = DALL_MW_r(r,1,t);
    end     
end
end
%{
%*****Generalized Max-Weight Policy with Logarithm
AvgRate_MW_Log=zeros(N,Ttot);  %AvgRate(i,t) denotes average throughput up to time t
Schedule_MW_Log=zeros(N,Ttot);
Receive_MW_Log=zeros(N,Ttot);
VQueue_MW_Log=zeros(N,Ttot+1);
X_MW_Log=zeros(N,Ttot+1);
D_MW_Log=zeros(N,Ttot+1);
Z_MW_Log=zeros(1,Ttot+1);
XALL_MW_Log=zeros(1,Ttot+1);
DALL_MW_Log=zeros(1,Ttot+1);

%MW: Slot-wise update
for t=1:Ttot
    head = 0;   %Pick the one with ON channel and smallest Xn(t)
    Dbmax = -1e10;   %Track minimum Xn(t) in time slot t
    Ratemax = 0;
    best = 0;
    for i=1:N
        if InsRate(i,t) > Ratemax
            best = i;
            Ratemax = InsRate(i,t);
        end
    end
    
    for i=1:N
        X_MW_Log(i,t+1)=X_MW_Log(i,t)-qn(1,i);
        if t == 1
            if InsRate(i,t) > 0 && InsRate(i,t) > Dbmax
                head = i;
                Dbmax = InsRate(i,t);
            end    
        else
            VQueue_MW_Log(i,t) = InsRate(i,t)*(log(1+subplus(-X_MW_Log(i,t))))*wn(1,i);
%            VQueue_MW_Log(i,t) = InsRate(i,t)*(-X_MW_Log(i,t))*wn(1,i);
            if (InsRate(i,t) > 0) && (VQueue_MW_Log(i,t) > Dbmax)
                Dbmax = VQueue_MW_Log(i,t);
                head = i;
            end
        end
    end
%    z = rand;
    if head > 0 
        Schedule_MW_Log(head,t)=1;
        %if z <= InsRate(head,t)
        X_MW_Log(head,t+1)=X_MW_Log(head,t+1)+InsRate(head,t);
        Receive_MW_Log(head,t)=InsRate(head,t);
        %end
    end
    
    if head > 0 && best > 0 && InsRate(head,t) < InsRate(best,t)
        Z_MW_Log(1,t+1) = Z_MW_Log(1,t) + 1;
    else
        Z_MW_Log(1,t+1) = Z_MW_Log(1,t);
    end
    
    for i=1:N
        if i == head && Receive_MW_Log(i,t) > 0
            if t == 1
                AvgRate_MW_Log(i,t)=Receive_MW_Log(i,t);
            else
                %AvgRate_MW(i,t)=AvgRate_MW(i,t-1)/t*(t-1)+1/t;
                AvgRate_MW_Log(i,t)=AvgRate_MW_Log(i,t-1)+Receive_MW_Log(i,t);
            end
        else
            if t > 1
                %AvgRate_MW(i,t)=AvgRate_MW(i,t-1)/t*(t-1);
                AvgRate_MW_Log(i,t)=AvgRate_MW_Log(i,t-1);
            end
        end
        if X_MW_Log(i,t+1)<0 && floor(-X_MW_Log(i,t+1)/qn(1,i))>D_MW_Log(i,t)
            D_MW_Log(i,t+1) = floor(-X_MW_Log(i,t+1)/qn(1,i));
        else
            D_MW_Log(i,t+1) = D_MW_Log(i,t);
        end            
    end
%Find X(t) and D(t)
    for i=1:N
        XALL_MW_Log(1,t+1) = XALL_MW_Log(1,t+1) + X_MW_Log(i,t+1);
    end
    if XALL_MW_Log(1,t+1)<0 && floor(-XALL_MW_Log(1,t+1))>DALL_MW_Log(1,t)
        DALL_MW_Log(1,t+1) = floor(-XALL_MW_Log(1,t+1));
    else
        DALL_MW_Log(1,t+1) = DALL_MW_Log(1,t);
    end     
end
%}
%*****Proportional Fair Policy
AvgRate_PF=zeros(N,Ttot);  %AvgRate(i,t) denotes average throughput up to time t
Schedule_PF=zeros(N,Ttot);
Receive_PF=zeros(N,Ttot);
VQueue_PF=zeros(N,Ttot+1);
X_PF=zeros(N,Ttot+1);
D_PF=zeros(N,Ttot+1);
Z_PF=zeros(1,Ttot+1);

%***PF: Slot-wise update
for t=1:Ttot
    choice = 0; %Pick the one with ON channel and largest InsRate(t)/AvgRate(t)
    Rmax = -1; %Track maximum InsRate(t)/AvgRate(t) in time slot t
    best = 0;
    Ratemax = 0;
    for i=1:N
        if InsRate(i,t) > Ratemax
            best = i;
            Ratemax = InsRate(i,t);
        end
    end
    for i=1:N
        X_PF(i,t+1)=X_PF(i,t)-qn(1,i);
        if t == 1
            if InsRate(i,t) > 0 && InsRate(i,t) > Rmax
                choice = i;
            end
        else
            VQueue_PF(i,t) = qn(1,i)*InsRate(i,t)/(AvgRate_PF(i,t-1));
            if InsRate(i,t) > 0 && (qn(1,i)*InsRate(i,t)/(AvgRate_PF(i,t-1))) > Rmax
                Rmax = qn(1,i)*InsRate(i,t)/(AvgRate_PF(i,t-1));
                choice = i;
            end
        end
    end
    %w = rand;
    if choice > 0 
        Schedule_PF(choice,t)=1;
        %if w <= InsRate(choice,t)
            X_PF(choice,t+1) = X_PF(choice,t+1)+InsRate(choice,t);
            Receive_PF(choice,t) = InsRate(choice,t);
        %end
    end
    
    if choice > 0 && best > 0 && InsRate(choice,t) < InsRate(best,t)
        Z_PF(1,t+1) = Z_PF(1,t) + 1;
    else
        Z_PF(1,t+1) = Z_PF(1,t);
    end
    
    for i=1:N
        if i == choice && Receive_PF(i,t) > 0
            if t == 1
                AvgRate_PF(i,t)=Receive_PF(i,t);
            else
                %AvgRate_PF(i,t)=AvgRate_PF(i,t-1)/t*(t-1)+1/t;
                AvgRate_PF(i,t)=AvgRate_PF(i,t-1)+Receive_PF(i,t);
            end
        else
            if t > 1
                %AvgRate_PF(i,t)=AvgRate_PF(i,t-1)/t*(t-1);
                AvgRate_PF(i,t)=AvgRate_PF(i,t-1);
            end
        end
        if X_PF(i,t+1)<0 && floor(-X_PF(i,t+1)/qn(1,i))>D_PF(i,t)
            D_PF(i,t+1) = floor(-X_PF(i,t+1)/qn(1,i));
        else
            D_PF(i,t+1) = D_PF(i,t);
        end            
    end
end
 

%*****Highest Data Rate First Policy
AvgRate_HDR=zeros(N,Ttot);  %AvgRate(i,t) denotes average throughput up to time t
Schedule_HDR=zeros(N,Ttot);
Receive_HDR=zeros(N,Ttot);
VQueue_HDR=zeros(N,Ttot+1);
X_HDR=zeros(N,Ttot+1);
X_HDRsh = zeros(N, Ttot+1);
D_HDR=zeros(N,Ttot+1);
Z_HDR=zeros(1,Ttot+1);
XALL_HDR=zeros(1,Ttot+1);
DALL_HDR=zeros(1,Ttot+1);
%HDR: Slot-wise update
for t=1:Ttot
    head = 0;   %Pick the one with ON channel and smallest Xn(t)
    Dbmax = -1e10;   %Track minimum Xn(t) in time slot t
    Ratemax = 0;
    best = 0;
    condition = 0;
    ratetrack = 0;
    for i=1:N
        if InsRate(i,t) > Ratemax
            best = i;
            Ratemax = InsRate(i,t);
        end
        if  InsRate(i,t) > 0 %&& X_HDR(i,t) < 0
            condition = 1;
        end
    end
    
    for i=1:N
        X_HDR(i,t+1)=X_HDR(i,t)-qn(1,i);
        if t == 1
            if InsRate(i,t) > 0 && InsRate(i,t) > Dbmax
                head = i;
                Dbmax = InsRate(i,t);
            end    
        else
            VQueue_HDR(i,t) = InsRate(i,t)*(-X_HDR(i,t))*wn(1,i) ;
            if condition > 0 % At least one client with rate>0 and Xn(t) < 0
                %if ((InsRate(i,t) > ratetrack) || (InsRate(i,t) == ratetrack && Dbmax < VQueue_HDR(i,t)))
                if (InsRate(i,t)*max(1, X_HDRsh(i,t)) > ratetrack) || ((InsRate(i,t)*max(1, X_HDRsh(i,t)) == ratetrack) && Dbmax < VQueue_HDR(i,t))
                    Dbmax = InsRate(i,t)*(-X_HDR(i,t))*wn(1,i);
                    %ratetrack = InsRate(i,t);
                    ratetrack = InsRate(i,t)*max(1, X_HDRsh(i,t));
                    head = i;
                end
            else
                if (InsRate(i,t) > 0) && (InsRate(i,t)*(-X_HDR(i,t))*wn(1,i) > Dbmax)
                    Dbmax = InsRate(i,t)*(-X_HDR(i,t))*wn(1,i);
                    head = i;
                end
            end
        end
    end
%    z = rand;
    if head > 0 
        Schedule_HDR(head,t)=1;
        %if z <= InsRate(head,t)
        X_HDR(head,t+1)=X_HDR(head,t+1)+InsRate(head,t);
        Receive_HDR(head,t)=InsRate(head,t);
        %end
    end
    for i=1:N
        %X_HDRsh(i, t+1) =  -X_HDR(i, t+1)-alpha(1,i)*nthroot(t^3, 4);
        X_HDRsh(i, t+1) =  nthroot(subplus(-X_HDR(i, t+1)-alpha(1,i)*nthroot(t^3, 4)-beta(1,i)*nthroot(t, 2)-offset_), 10);
        %X_HDRsh(i, t+1) =  nthroot(subplus(-X_HDR(i, t+1)-qn(1,i)*alpha(1,i)*nthroot(t^3, 4)-beta(1,i)*nthroot(t, 2)), 10);
    end
    if head > 0 && best > 0 && InsRate(head,t) < InsRate(best,t)
        Z_HDR(1,t+1) = Z_HDR(1,t) + 1;
    else
        Z_HDR(1,t+1) = Z_HDR(1,t);
    end
    
    for i=1:N
        if i == head && Receive_HDR(i,t) > 0
            if t == 1
                AvgRate_HDR(i,t)=Receive_HDR(i,t);
            else
                AvgRate_HDR(i,t)=AvgRate_HDR(i,t-1)+Receive_HDR(i,t);
            end
        else
            if t > 1
                AvgRate_HDR(i,t)=AvgRate_HDR(i,t-1);
            end
        end
        if X_HDR(i,t+1)<0 && floor(-X_HDR(i,t+1)/qn(1,i))>D_HDR(i,t)
            D_HDR(i,t+1) = floor(-X_HDR(i,t+1)/qn(1,i));
        else
            D_HDR(i,t+1) = D_HDR(i,t);
        end            
    end
%Find X(t) and D(t)
    for i=1:N
        XALL_HDR(1,t+1) = XALL_HDR(1,t+1) + X_HDR(i,t+1);
    end
    if XALL_HDR(1,t+1)<0 && floor(-XALL_HDR(1,t+1))>DALL_HDR(1,t)
        DALL_HDR(1,t+1) = floor(-XALL_HDR(1,t+1));
    else
        DALL_HDR(1,t+1) = DALL_HDR(1,t);
    end     
end


%*****NOVA Policy
AvgRate_NOVA=zeros(N,Ttot);  %AvgRate(i,t) denotes average throughput up to time t
Schedule_NOVA=zeros(N,Ttot);
Receive_NOVA=zeros(N,Ttot);
VQueue_NOVA=zeros(N,Ttot+1);
X_NOVA=zeros(N,Ttot+1);
D_NOVA=zeros(N,Ttot+1);
Z_NOVA=zeros(1,Ttot+1);
W_NOVA=zeros(N,Ttot+1);
XALL_NOVA=zeros(1,Ttot+1);
DALL_NOVA=zeros(1,Ttot+1);
eps = 0.05;
tslot = 0.01;
%NOVA: Slot-wise update
for t=1:Ttot
    head = 0;   %Pick the one with ON channel and smallest Xn(t)
    Dbmax = -1e10;   %Track minimum Xn(t) in time slot t
    Ratemax = 0;
    best = 0;
    condition = 0;
    ratetrack = 0;
    for i=1:N
        if InsRate(i,t) > Ratemax
            best = i;
            Ratemax = InsRate(i,t);
        end
        if  InsRate(i,t) > 0
            condition = 1;
        end
    end
    
    for i=1:N
        X_NOVA(i,t+1)=X_NOVA(i,t)-qn(1,i);
        W_NOVA(i,t+1)=W_NOVA(i,t)+(eps*tslot);
        %W_NOVA(i,t+1)=W_NOVA(i,t)+(qn(1,i)*eps*tslot);
        if t == 1
            if InsRate(i,t) > 0 && InsRate(i,t) > Dbmax
                head = i;
                Dbmax = InsRate(i,t);
            end    
        else
            VQueue_NOVA(i,t) = InsRate(i,t)*(hNOVA(W_NOVA(i,t)))*wn(1,i) ;
            if condition > 0 % At least one client with rate>0
                if ((Dbmax < VQueue_NOVA(i,t)) || (InsRate(i,t) > ratetrack && Dbmax == VQueue_NOVA(i,t)))
                    Dbmax = VQueue_NOVA(i,t);
                    ratetrack = InsRate(i,t);
                    head = i;
                end
            end
        end
    end
%    z = rand;
    if head > 0 
        Schedule_NOVA(head,t)=1;
        %if z <= InsRate(head,t)
        X_NOVA(head,t+1)=X_NOVA(head,t+1)+InsRate(head,t);
        W_NOVA(head,t+1)=W_NOVA(head,t+1)-(eps*tslot*InsRate(head,t)/qn(1,head));
        %W_NOVA(head,t+1)=W_NOVA(head,t+1)-(eps*tslot*InsRate(head,t));
        Receive_NOVA(head,t)=InsRate(head,t);
        %end
    end
    
    if head > 0 && best > 0 && InsRate(head,t) < InsRate(best,t)
        Z_NOVA(1,t+1) = Z_NOVA(1,t) + 1;
    else
        Z_NOVA(1,t+1) = Z_NOVA(1,t);
    end
    
    for i=1:N
        if i == head && Receive_NOVA(i,t) > 0
            if t == 1
                AvgRate_NOVA(i,t)=Receive_NOVA(i,t);
            else
                AvgRate_NOVA(i,t)=AvgRate_NOVA(i,t-1)+Receive_NOVA(i,t);
            end
        else
            if t > 1
                AvgRate_NOVA(i,t)=AvgRate_NOVA(i,t-1);
            end
        end
        if X_NOVA(i,t+1)<0 && floor(-X_NOVA(i,t+1)/qn(1,i))>D_NOVA(i,t)
            D_NOVA(i,t+1) = floor(-X_NOVA(i,t+1)/qn(1,i));
        else
            D_NOVA(i,t+1) = D_NOVA(i,t);
        end            
    end
%Find X(t) and D(t)
    for i=1:N
        XALL_NOVA(1,t+1) = XALL_NOVA(1,t+1) + X_NOVA(i,t+1);
    end
    if XALL_NOVA(1,t+1)<0 && floor(-XALL_NOVA(1,t+1))>DALL_NOVA(1,t)
        DALL_NOVA(1,t+1) = floor(-XALL_NOVA(1,t+1));
    else
        DALL_NOVA(1,t+1) = DALL_NOVA(1,t);
    end     
end

%{
plot(Taxis,InsRate(1,:),'-*r');
hold on;
plot(Taxis,InsRate(2,:),'-ob');
figure;

plot(Taxis2,D_MW(1,:),'-*r');
hold on;
plot(Taxis2,D_PF(1,:),'-ob');
figure;
plot(Taxis2,D_MW(2,:),'-*r');
hold on;
plot(Taxis2,D_PF(2,:),'-ob');
%}
%*****Naive Random Policy
end

