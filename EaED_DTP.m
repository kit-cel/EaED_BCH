function [PEaEDsucc,PEaEDfail,PEaEDmc] = EaED_DTP(n_primitive,k,t,extended, Umax, Emax)
dmin = 2*t+1;
n = n_primitive;
filename = sprintf('./BCH_weight_enumerators/A%d_%d_%d.mat', n_primitive, k, dmin);

if exist(filename, 'file')
    Au = load(filename).weights;
else
    Au = zeros(1,n_primitive+1);
    b = log2(n_primitive+1);
    Au(1) = 1;
    for w = dmin:100
        Au(w+1) = 2^(-1*b*t)*nck(n_primitive,w);
    end
end

if extended==1
    n = n+1;
    dmin = dmin + 1;
    Au = [Au,0];
    for i=1:n_primitive/2
        Au(2*i+1) = Au(2*i+1) + Au(2*i);
        Au(2*i) = 0;
    end
end



%row is input error number, column is output error number
PBDDsucc = zeros(Umax+Emax+1,Umax+Emax+t+1);
PBDDfail = zeros(Umax+Emax+1,Umax+Emax+t+1);
PBDDmc = zeros(Umax+Emax+1,Umax+Emax+t+1);

for U=0:t
    PBDDsucc(U+1,1) = 1;
    PBDDfail(U+1,1) = 0;
    PBDDmc(U+1,1) = 0;
end


for U=t+1:Umax+Emax
    Pmc = 0;
    for r=max(0,U-t):U+t
        if Au(r+1) == 0
            continue;
        end  
        L=0;
        for a=0:t
            for b=0:t-a
                if U+a-b == r
                    L = L + Au(r+1)*nck(r,a)*nck(n-r,b);
                end
            end
        end
        PBDDmc(U+1,r+1) = L/nck(n,U);
        Pmc = Pmc + PBDDmc(U+1,r+1);
    end
    PBDDfail(U+1,U+1) = 1-Pmc;
end




PEaEDsucc = zeros(Umax+1,Emax+1, Umax+Emax+t+1);
%here, fail is recorded differently as U+E/2 is not alwasys integer
%so we only record the probability of decoding failure 
PEaEDfail = zeros(Umax+1,Emax+1); 
PEaEDmc = zeros(Umax+1,Emax+1, Umax+Emax+t+1);

%case for no erasure, E=0, EaED reduces to BDD
for U=0:Umax
    PEaEDfail(U+1,1) = PBDDfail(U+1,U+1);
    for R = 0:U+t
        PEaEDsucc(U+1,1,R+1) = PBDDsucc(U+1,R+1);    
        PEaEDmc(U+1,1,R+1) = PBDDmc(U+1,R+1);
    end
end

for U=0:t
    for E=1:dmin-1
        if 2*U + E < dmin
            PEaEDsucc(U+1,E+1,1) = 1;
        end
    end
end
     


for U=0:Umax
    for E=max(1,dmin-2*U):Emax
        %for cases in \mathcal{L}
        for E1 = 0:t-U
            Theta = nck(n,U)*nck(n-U,E) * nck(E,E1);       
            La = 0;
            Lb = 0;
            Lc = 0;
            for a=0:t
                for b=0:t-a
                    r = U+E-E1+a-b;
                    Lar = 0;
                    Lbr = 0;
                    Lcr = 0;
                    if r<dmin
                        continue;
                    end
                    term1 = Au(r+1) * nck(r,a)*nck(n-r,b);
                    for lambda = max(0, E1-a):min(E1,n-r-b)
                        for gamma = 0:min(E-E1,b)
                            if E-E1-gamma<=r-a
                                term2 = nck(n-r-b,lambda) * nck(b,gamma) * nck(r-a,E-E1-gamma)*nck(a,E1-lambda);
                                if lambda-gamma < U+E1-a-b
                                    Lar = Lar + term1 * term2;
                                elseif lambda-gamma > U+E1-a-b
                                    Lbr = Lbr + term1 * term2;
                                else
                                    Lcr = Lcr + term1 * term2;
                                end
                            end
                        end
                    end

                    PEaEDmc(U+1,E+1,r+1) = PEaEDmc(U+1,E+1,r+1) + 2*((Lar + Lcr/2)/Theta)*nck(E,E1)/(2^E);
                    PEaEDsucc(U+1,E+1,1) = PEaEDsucc(U+1,E+1,1) + 2*((Lbr + Lcr/2)/Theta)*nck(E,E1)/(2^E);
                    La = La + Lar;
                    Lb = Lb + Lbr;
                    Lc = Lc + Lcr;
                end
            end
            
            L1 = Theta - (La + Lb + Lc);
            
            PEaEDsucc(U+1,E+1,1) = PEaEDsucc(U+1,E+1,1) + 2*(L1/Theta)*nck(E,E1)/(2^E);
        end
        
        
        
        
        %for cases in \mathcal{M}
        for E1 = max(0,t-U+1) : min(U+E-t-1,E)
            PEaEDfail(U+1,E+1) = PEaEDfail(U+1,E+1) + PBDDfail(U+E1+1,U+E1+1) * PBDDfail(U+E-E1+1,U+E-E1+1) *nck(E,E1)/(2^E);
            for R=0:max(U+E1+t,U+E-E1+t)
                pmc = PBDDfail(U+E1+1,U+E1+1)*PBDDmc(U+E-E1+1,R+1) + PBDDmc(U+E1+1,R+1)*PBDDfail(U+E-E1+1,U+E-E1+1);               
                PEaEDmc(U+1,E+1,R+1) = PEaEDmc(U+1,E+1,R+1) + pmc *nck(E,E1)/(2^E);
            end
            
            for R1=U+E1-t:U+E1+t
                for R2=U+E-E1-t:U+E-E1+t
                    PEaEDmc(U+1,E+1,R1+1) = PEaEDmc(U+1,E+1,R1+1) + PBDDmc(U+E1+1,R1+1)*PBDDmc(U+E-E1+1,R2+1)/2*nck(E,E1)/(2^E);
                    PEaEDmc(U+1,E+1,R2+1) = PEaEDmc(U+1,E+1,R2+1) + PBDDmc(U+E1+1,R1+1)*PBDDmc(U+E-E1+1,R2+1)/2*nck(E,E1)/(2^E);
                   
                end
            end
            
            
        end
    end
end

end

    

