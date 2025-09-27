function [ber, fer] = BDD_BER(DTP, n, p)
    ber = 0;
    fer = 0;
    for U=0:size(DTP,1)-1
        
        for R=0:size(DTP,1)-1
             ber = ber + R*nck(n,U)*(p^U)*((1-p)^(n-U))*DTP(U+1,R+1);
             
        end
        fer = fer + nck(n,U)*(p^U)*((1-p)^(n-U))*(1-DTP(U+1,1));
    end
    ber = ber/n;
end

