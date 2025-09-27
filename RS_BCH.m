clear all;
clc;


M = 8;
m = 64;
B = 10;
b = 10;
n_primitive = 2^b-1;
N = 544;
K = 514;
T = 15;
n_shortened = 700;
shorten = n_primitive - n_shortened;
k = 680;
t = 2;
dmin = 2*t+1;

n = n_primitive;
t=2;
dmin = 2*t+1;
k = n - log2(n_primitive+1)*t;
extended = 0;

Umax = 25;
Emax = 15;
if extended==1
    n = n_primitive+1;
    dmin = dmin + 1;
end

rate = B*K*M/(n_shortened*m);


% syms x y
% f = (1+y*((1+x)^10-1))^8*(1+x)^(n-80);
% expanded_f = expand(f);
% [coeff_list, term_list] = coeffs(f);
% 
% coeff_2dv1 = zeros(41,9);
% for u=0:40
%     for v=0:8
%         target = x^u*y^v;
%         idx = find(term_list == target);
% 
%         if ~isempty(idx)
%             c = coeff_list(idx);
%         else
%             fprintf("didn't find this term");
%             c = 0; % term not found
%         end
%         disp(c);
%         coeff_2dv1(u+1,v+1) = c;
%         Pv (v+1) = Pv (v+1)  + Pu(u+1)*c / nck(n,u);
%     end
% end
% 
% 
% Pv2 = zeros(1,10);
% syms x y
% f = (1+y*((1+x)^10-1))^9*(1+x)^(n-90);
% expanded_f = expand(f);
% [coeff_list, term_list] = coeffs(f);
% coeff_2dv2 = zeros(41,9);
% for u=0:40
%     for v=0:9
%         target = x^u*y^v;
%         idx = find(term_list == target);
% 
%         if ~isempty(idx)
%             c = coeff_list(idx);
%         else
%             fprintf("didn't find this term");
%             c = 0; % term not found
%         end
%         disp(c);
%         Pv2 (v+1) = Pv2 (v+1)  + Pu(u+1)*c / nck(n,u);
%         coeff_2dv2(u+1,v+1) = c;
%     end
% end
% 
% save("coeff_2dv1");
% save("coeff_2dv2");

coeff_2dv1 = load("coeff_2dv1").coeff_2dv1;
coeff_2dv2 = load("coeff_2dv2").coeff_2dv2;



Te=0.06;
Ta =0.56;
EbNo_list = [6:0.1:8.01];
BDD_ber_list = [];
uncoded_ber_list = [];
EaEDa_ber_list = [];
FER_list = [];
minFER = 100;
% for Te = [0.01:0.01:0.18]
%     for Ta = [Te+0.1:0.01:1.4]
for EbNo_idx = 1:length(EbNo_list)
    EbNo = EbNo_list(EbNo_idx);
    EbNo_val = 10^(EbNo/10);
    EsNo_val = EbNo_val * rate;
    variance = 0.5 / EsNo_val;
    sigma = sqrt(variance);
    delta = 1-qfunc((-Te-1)/sigma);
    ep = 1-qfunc((Te-1)/sigma) - delta;

    cross_over_prob = 1-qfunc((-1)/sigma);    
    uncoded_ber_list = [uncoded_ber_list, cross_over_prob];

    Pca = (qfunc((Ta-1)/sigma)) / (1-ep-delta);
    Pwa = (1-qfunc((-Ta-1)/sigma)) / (delta);
    
    
    [PEaEDasucc,PEaEDafail,PEaEDamc] = EaED_w_anchor_DTP(n_primitive,n,k,t,extended, Umax, Emax, Pca, Pwa, shorten);
    DTP = PEaEDasucc+PEaEDamc;
%     format longE
%      sum(PEaEDasucc+PEaEDamc,3)+PEaEDafail
    [EaEDa_ber, EaEDa_fer] = EaED_error_rate(PEaEDasucc+PEaEDamc, PEaEDafail, n, delta, ep);
    EaEDa_ber_list = [EaEDa_ber_list EaEDa_ber];
    
    Pu = zeros(1,Umax+Emax+t+1);
    totalP = 0;
    for U=0:Umax
        for E=0:Emax
            pe = nck(n_shortened,U)*nck(n_shortened-U,E)*(delta^U)*(ep^E)*((1-delta-ep)^(n_shortened-U-E));
            totalP = totalP + pe;
            for R=0:Umax+Emax+t
                 Pu(R+1) = Pu(R+1)+ pe*DTP(U+1,E+1,R+1);
            end
            %erasures are randomly assigned hard-decision value
            for e = 0:E
                Pu(U+e+1) = Pu(U+e+1) + pe*PEaEDafail(U+1,E+1)*nck(E,e)/(2^E);
            end
        end
    end

    
    
    FER = RS_FER(Pu, Umax+Emax, M, n_shortened, T, coeff_2dv1, coeff_2dv2);
    
%     if FER<minFER
%         minFER = FER;
%         fprintf('Te = %f, Ta = %f, FER = %.8e\n', Te, Ta, FER);
%     end
   
    FER_list = [FER_list, FER];
    fprintf('%f, %e\n', EbNo, FER);
end


%     end
% end
% semilogy(EbNo_list,uncoded_ber_list, 'k', 'LineWidth', 1.5);

semilogy(uncoded_ber_list, FER_list, 'rx-', 'LineWidth', 1.5);
hold on;
grid on;

%simulated
Pu_matrix = [0.65434 0.00224 0.02357 0.08894 0.10723 0.07127 0.033 0.01314 0.00437 0.00136 0.00043 9e-05 1e-05 1e-05 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 %EbNo6
0.70946 0.00217 0.02109 0.08055 0.09287 0.05681 0.0244 0.009 0.00259 0.00078 0.00024 2e-05 2e-05 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 %EbNo6.1
0.76007 0.00207 0.01885 0.07111 0.07772 0.04332 0.01834 0.00596 0.00202 0.00042 8e-05 3e-05 0 1e-05 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 %EbNo6.2
0.80703 0.00175 0.01627 0.06195 0.06353 0.03313 0.01131 0.00371 0.00103 0.00027 2e-05 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 %EbNo6.3
0.84613 0.00176 0.01425 0.05172 0.05068 0.02438 0.00821 0.00219 0.00055 9e-05 4e-05 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 %EbNo6.4
0.88144 0.00122 0.01127 0.04236 0.03871 0.01769 0.00549 0.00141 0.00031 0.0001 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 %EbNo6.5
0.90991 0.00102 0.00929 0.03374 0.02964 0.01199 0.0033 0.00091 0.00016 3e-05 1e-05 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 %EbNo6.6
0.93303 0.00086 0.00734 0.02713 0.02068 0.0084 0.00205 0.00045 6e-05 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 %EbNo6.7
0.94963 0.00058 0.00571 0.02151 0.0155 0.0055 0.0013 0.00024 2e-05 1e-05 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 %EbNo6.8
0.96389 0.00054 0.0039 0.01619 0.01085 0.0037 0.00075 0.00016 2e-05 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 %EbNo6.9
0.97322 0.00033 0.00305 0.01226 0.00789 0.00263 0.00056 5e-05 0 1e-05 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 %EbNo7
0.98146 0.00025 0.0024 0.00861 0.00543 0.0015 0.0003 5e-05 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 %EbNo7.1
0.98666 0.0003 0.00157 0.00644 0.00395 0.00092 0.00015 1e-05 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 %EbNo7.2
0.99134 0.00022 0.00116 0.00433 0.00234 0.00057 2e-05 1e-05 1e-05 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 %EbNo7.3
0.99417 0.00015 0.00066 0.003 0.00169 0.00031 2e-05 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 %EbNo7.4
0.99638 0.00013 0.00056 0.00183 0.00086 0.00021 3e-05 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; %EbNo7.5
FER_list = [];
for EbNo_idx = 1:16
    EbNo = EbNo_list(EbNo_idx);
    EbNo_val = 10^(EbNo/10);
    EsNo_val = EbNo_val * rate;
    variance = 0.5 / EsNo_val;
    sigma = sqrt(variance);
    delta = 1-qfunc((-Te-1)/sigma);
    ep = 1-qfunc((Te-1)/sigma) - delta;
    cross_over_prob = 1-qfunc((-1)/sigma);    
    uncoded_ber_list = [uncoded_ber_list, cross_over_prob];
    Pu = Pu_matrix(EbNo_idx,:); 
    FER = RS_FER(Pu, 20, M, n_shortened, T, coeff_2dv1, coeff_2dv2);
    FER_list = [FER_list, FER];
    fprintf('%f, %e\n', EbNo, FER);
end
semilogy(uncoded_ber_list(1:16), FER_list, 'kx-', 'LineWidth', 1.5);

