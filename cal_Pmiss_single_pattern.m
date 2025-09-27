function [Pfail, Pmc] = cal_Pmiss_single_pattern(n,t,Au, Pca, Pwa, Umax, Emax, shorten)
%consider the case where D+e1>t, and we know for sure that BDD will fail or
%miscorrect, calcualte the miscorrection probability of a single pattern
Ul = [0:Umax];
El = [1:Emax];
e1l = [0:Emax];
Pfail = zeros(Umax+1, Emax+1, Emax+1);
Pmc = zeros(Umax+1, Emax+1, Emax+1, Umax+Emax+t+1);
for U=Ul
    for E=El
        Pcas = (shorten)/(n - U -E);
        for e1 = 0:E
            Theta = nck(n,U)*nck(n-U,E)*nck(E,e1);
            if U+e1<=t
                continue;
            end
                
            for a=0:t
                for b=0:t-a
                    wc = U+e1+a-b;
                    if wc<t
                        continue;
                    end
                    if (Au(wc+1)==0)
                        continue;
                    end
                    term1 = Au(wc+1)* nck(wc,a)*nck(n-wc,b);
                    
                    for nwa = 0:b
                        for nca = 0:a
                            c = e1-(b-nwa);
                            d = E-e1-(a-nca);
                            if c>=0 && c<=(wc-a) && d>=0 && d<=n-wc-b
                                term2 = nck(b,nwa)*nck(a,nca)*nck(wc-a,c)*nck(n-wc-b,d);
                                term2 = term2*((1-Pca)^nca)*((1-Pcas)^nca)*((1-Pwa)^nwa);
                                Lr = term1*term2;
                            end
                            
                        end
                    end
                    Pmc(U+1,E+1,e1+1,wc+1) = Pmc(U+1,E+1,e1+1,wc+1) + Lr / Theta;
                end
            end
            Pfail(U+1,E+1, e1+1) = 1 - sum(Pmc(U+1,E+1,e1+1,:),'all');
        end
        
    end
end

end

