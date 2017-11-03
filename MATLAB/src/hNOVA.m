function [hout] = hNOVA(b)
%Function for assessing rebuffering
hout = 0;
b0 = 0;
bi = b + b0;
    if bi > 20
        hout = 0.005*(20*bi + ((bi-20)*20)^2); 
    end
    if bi <= 20 && bi >= 0
        hout = 0.1*bi;
    end
    if bi < 0
        hout = 0;
    end
end
