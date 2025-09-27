function [FER] = RS_FER(Pu, Umax, M, n, T, coeff_2dv1, coeff_2dv2)
 %Now consider the number of byte errors from one RS cw given the number of bit errors u in
    %a BCH cw. Assuming that in each BCH cw, one RS cw constribute 8 or 9
    %symbols (bytes)

    %If contribute 8 bytes
    Pv = zeros(1,9);
    for u=0:Umax
        for v=0:8
             Pv (v+1) = Pv (v+1)  + Pu(u+1)*coeff_2dv1(u+1,v+1) / nck(n,u);
        end
    end


    Pv2 = zeros(1,10);
    for u=0:Umax
        for v=0:9
            Pv2 (v+1) = Pv2 (v+1)  + Pu(u+1)*coeff_2dv2(u+1,v+1) / nck(n,u);
        end
    end



    Pe = Pv;
    for i=1:31
        Pe = conv(Pe,Pv);

    end


    for i=1:32
        Pe = conv(Pe,Pv2);
    
    end


    FER_single_frame = sum(Pe(T+2:3*T+2));
    FER = 1 - (1-FER_single_frame)^M;
end

