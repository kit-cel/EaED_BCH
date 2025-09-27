clear all;
clc;
n_primitive= 255;
n = n_primitive;
t=2;
dmin = 2*t+1;
k = n - log2(n+1)*t;
extended = 1;
Umax = dmin;
Emax = dmin*2;
if extended==1
    n = n_primitive+1;
    dmin = dmin + 1;
end
rate = k/n;


[PBDDsucc,PBDDfail,PBDDmc] = BDD_DTP(n_primitive,k,t,extended, Umax);
BDDDTP = PBDDsucc+PBDDfail+PBDDmc;
[PEaEDsucc,PEaEDfail,PEaEDmc] = EaED_DTP(n_primitive,k,t,extended, Umax, Emax);

format shortE %check if the probabilities sum up to one
EaEDDTP = PEaEDsucc+PEaEDmc;
sum(EaEDDTP,3)+PEaEDfail


EbNo = 9;
EbNo_val = 10^(EbNo/10);
EsNo_val = EbNo_val * rate;
variance = 0.5 / EsNo_val;
sigma = sqrt(variance);
format short


%for EaED, grid serach for optimal T, fix EbNo = 9dB
T_list = [0.01:0.01:0.2];
minBER = 100;
Topt = 0;
Taopt = 0;
for T = T_list
    delta = 1-qfunc((-T-1)/sigma);
    ep = 1-qfunc((T-1)/sigma) - delta;        
    [EaED_ber, EaED_fer] = EaED_error_rate(EaEDDTP, PEaEDfail, n, delta, ep);
    if EaED_ber < minBER
        minBER = EaED_ber;
        Topt = T
    end

end
T = Topt;
labels = {'uncoded', 'BDD'};
labels = [labels, sprintf('EaED (T = %.2f)', T)];

EbNo_list = [6:0.2:10.01];
BDD_ber_list = [];
uncoded_ber_list = [];
EaED_ber_list = [];
for EbNo_idx = 1:length(EbNo_list)
    EbNo = EbNo_list(EbNo_idx);
    EbNo_val = 10^(EbNo/10);
    EsNo_val = EbNo_val * rate;
    variance = 0.5 / EsNo_val;
    sigma = sqrt(variance);
    delta = 1-qfunc((-T-1)/sigma);
    ep = 1-qfunc((T-1)/sigma) - delta;

    cross_over_prob = 1-qfunc((-1)/sigma);    
%     [BDD_ber, BDD_fer] = BDD_error_rate(BDDDTP, n, cross_over_prob);
%     [EaED_ber, EaED_fer] = EaED_error_rate(EaEDDTP, PEaEDfail, n, delta, ep);
     [BDD_ber, BDD_fer] = BDD_error_rate(PBDDmc, n, cross_over_prob);
     [EaED_ber, EaED_fer] = EaED_error_rate(PEaEDmc, zeros(size(PEaEDfail)), n, delta, ep);
%     [BDD_ber, BDD_fer] = BDD_error_rate(PBDDfail, n, cross_over_prob);
%     [EaED_ber, EaED_fer] = EaED_error_rate(zeros(size(PEaEDmc)), PEaEDfail, n, delta, ep);
    BDD_ber_list = [BDD_ber_list, BDD_ber];
    EaED_ber_list = [EaED_ber_list, EaED_ber];
    uncoded_ber_list = [uncoded_ber_list, cross_over_prob];
end




semilogy(EbNo_list,uncoded_ber_list, 'k', 'LineWidth', 1.5);
hold on;
grid on;
semilogy(EbNo_list, BDD_ber_list, 'bo', 'LineWidth', 1.5);
semilogy(EbNo_list, EaED_ber_list, 'rx-', 'LineWidth', 1.5);

%for EaEDa, grid search for the optimal T and Ta, fix EbNo = 9dB

T_list = [0.08:0.01:0.16];
Ta_list = [0.1:0.01:1.8];
minBER = 100;
Topt = 0;
Taopt = 0;
EbNo = 9;
EbNo_val = 10^(EbNo/10);
EsNo_val = EbNo_val * rate;
variance = 0.5 / EsNo_val;
sigma = sqrt(variance);

for T = T_list
    delta = 1-qfunc((-T-1)/sigma);
    ep = 1-qfunc((T-1)/sigma) - delta;
    for Ta = Ta_list
        if Ta<=T
            continue;
        end
        Pca = (qfunc((Ta-1)/sigma)) / (1-ep-delta);
        Pwa = (1-qfunc((-Ta-1)/sigma)) / (delta);
        [PEaEDasucc,PEaEDafail,PEaEDamc] = EaED_w_anchor_DTP(n_primitive,n,k,t,extended, Umax, Emax, Pca, Pwa, 0);
        [EaEDa_ber, EaEDa_fer] = EaED_error_rate(PEaEDasucc+PEaEDamc, PEaEDafail, n, delta, ep);
        if EaEDa_ber < minBER
            minBER = EaEDa_ber;
            Topt = T
            Taopt = Ta
        end
    end
end

T = Topt;
Ta = Taopt;
labels = [labels, sprintf('EaEDa (T = %.2f, Ta = %.2f)', T, Ta)];

format short
EaEDa_ber_list = zeros(length(EbNo_list));
for EbNo_idx = 1:length(EbNo_list)
    EbNo = EbNo_list(EbNo_idx);
    EbNo_val = 10^(EbNo/10);
    EsNo_val = EbNo_val * rate;
    variance = 0.5 / EsNo_val;
    sigma = sqrt(variance);
    delta = 1-qfunc((-T-1)/sigma);
    ep = 1-qfunc((T-1)/sigma) - delta;
    Pca = (qfunc((Ta-1)/sigma)) / (1-ep-delta);
    Pwa = (1-qfunc((-Ta-1)/sigma)) / (delta);
    [PEaEDasucc,PEaEDafail,PEaEDamc] = EaED_w_anchor_DTP(n_primitive,n,k,t,extended, Umax, Emax, Pca, Pwa, 0);
%     sum(PEaEDasucc+PEaEDamc,3)+PEaEDafail
%     [EaEDa_ber, EaEDa_fer] = EaED_error_rate(PEaEDasucc+PEaEDamc, PEaEDafail, n, delta, ep);
     [EaEDa_ber, EaEDa_fer] = EaED_error_rate(PEaEDamc, zeros(size(PEaEDafail)), n, delta, ep);
%     [EaEDa_ber, EaEDa_fer] = EaED_error_rate(zeros(size(PEaEDamc)), PEaEDafail, n, delta, ep);
    EaEDa_ber_list(EbNo_idx) = EaEDa_ber;

end

semilogy(EbNo_list, EaEDa_ber_list, 'mo-', 'LineWidth', 1.5);

legend(labels);


for EbNo_idx = 1:length(EbNo_list)
    EbNo = EbNo_list(EbNo_idx);
    BER = EaEDa_ber_list(EbNo_idx);
    fprintf('%f %e\n', EbNo, BER);
end