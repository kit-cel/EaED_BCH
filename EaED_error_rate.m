function [ber,fer] = EaED_BER(DTP, Pfail, n, delta, ep)
    ber = 0;
    fer = 0;
    for U=0:size(DTP,1)-1
        for E=0:size(DTP,2)-1
            pe = nck(n,U)*nck(n-U,E)*(delta^U)*(ep^E)*((1-delta-ep)^(n-U-E));
%             ber = ber + (U+E/2)*pe;
            for R=0:size(DTP,3)-1
                 ber = ber + R*pe*DTP(U+1,E+1,R+1);                
            end
            ber = ber + (U+E/2)*pe*Pfail(U+1,E+1);
            fer = fer + pe * (1-DTP(U+1,E+1,1));
        end
    end
    ber = ber / n;

end

