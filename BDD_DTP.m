function [PBDDsucc,PBDDfail,PBDDmc] = BDD_DTP(n_primitive,k,t,extended, Umax)
dmin = 2*t+1;
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
    n_primitive = n_primitive+1;
    dmin = dmin + 1;
    Au = [Au,0];
    for i=1:n_primitive/2
        Au(2*i+1) = Au(2*i+1) + Au(2*i);
        Au(2*i) = 0;
    end
end



%row is input error number, column is output error number
PBDDsucc = zeros(Umax+1,Umax+t+1);
PBDDfail = zeros(Umax+1,Umax+t+1);
PBDDmc = zeros(Umax+1,Umax+t+1);

for U=0:t
    PBDDsucc(U+1,1) = 1;
    PBDDfail(U+1,1) = 0;
    PBDDmc(U+1,1) = 0;
end


for U=t+1:Umax
    Pmc = 0;
    for r=max(0,U-t):U+t
        if Au(r+1) == 0
            continue;
        end  
        L=0;
        for a=0:t
            for b=0:t-a
                if U+a-b == r
                    L = L + Au(r+1)*nck(r,a)*nck(n_primitive-r,b);
                end
            end
        end
        PBDDmc(U+1,r+1) = L/nck(n_primitive,U);
        Pmc = Pmc + PBDDmc(U+1,r+1);
    end
    PBDDfail(U+1,U+1) = 1-Pmc;
end





end

