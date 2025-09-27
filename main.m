%verifying the BER reusults of EaEDa with different Ta with the
%simulation results
clear all;
clc;

n_primitive= 255;
n = n_primitive;
t=2;
dmin = 2*t+1;
k = n - log2(n+1)*t;
extended = 0;
Umax = dmin;
Emax = dmin*2;
if extended==1
    n = n_primitive+1;
    dmin = dmin + 1;
end
rate = k/n;

%make DTP
[PBDDsucc,PBDDfail,PBDDmc] = BDD_DTP(n_primitive,k,t,extended, Umax);
BDDDTP = PBDDsucc+PBDDfail+PBDDmc;
[PEaEDsucc,PEaEDfail,PEaEDmc] = EaED_DTP(n_primitive,k,t,extended, Umax, Emax);

format short
EaEDDTP = PEaEDsucc+PEaEDmc;

% sum(EaEDDTP,3)+PEaEDfail


%verifying cases with some random T and Ta values
T=0.12;
Ta_list = [0.2:0.1:1.0];

EbNo_list = [4:0.2:10.01];
BDD_ber_list = [];
uncoded_ber_list = [];
EaED_ber_list = [];
EaEDa_ber_list = zeros(length(Ta_list),length(EbNo_list));
for EbNo_idx = 1:length(EbNo_list)
    EbNo = EbNo_list(EbNo_idx);
    EbNo_val = 10^(EbNo/10);
    EsNo_val = EbNo_val * rate;
    variance = 0.5 / EsNo_val;
    sigma = sqrt(variance);
    delta = 1-qfunc((-T-1)/sigma);
    ep = 1-qfunc((T-1)/sigma) - delta;

    cross_over_prob = 1-qfunc((-1)/sigma);    
    [BDD_ber, BDD_fer] = BDD_error_rate(BDDDTP, n, cross_over_prob);
    [EaED_ber, EaED_fer] = EaED_error_rate(EaEDDTP, PEaEDfail, n, delta, ep);
    BDD_ber_list = [BDD_ber_list, BDD_ber];
    EaED_ber_list = [EaED_ber_list, EaED_ber];
    uncoded_ber_list = [uncoded_ber_list, cross_over_prob];
    
    for Ta_idx = 1:length(Ta_list)
        Ta = Ta_list(Ta_idx);
        Pca = (qfunc((Ta-1)/sigma)) / (1-ep-delta);
        Pwa = (1-qfunc((-Ta-1)/sigma)) / (delta);
        [PEaEDasucc,PEaEDafail,PEaEDamc] = EaED_w_anchor_DTP(n_primitive,n,k,t,extended, Umax, Emax, Pca, Pwa, 0);
        sum(PEaEDasucc+PEaEDamc,3)+PEaEDafail
        [EaEDa_ber, EaEDa_fer] = EaED_error_rate(PEaEDasucc+PEaEDamc, PEaEDafail, n, delta, ep);
        EaEDa_ber_list(Ta_idx, EbNo_idx) = EaEDa_ber;
    end

end

semilogy(EbNo_list,uncoded_ber_list, 'k', 'LineWidth', 1.5);
hold on;
grid on;
semilogy(EbNo_list, BDD_ber_list, 'bo', 'LineWidth', 1.5);
semilogy(EbNo_list, EaED_ber_list, 'rx-', 'LineWidth', 1.5);
% semilogy(EbNo_list, EaEDa_ber_list, 'cx-', 'LineWidth', 1.5);

% labels = {'uncoded', 'BDD', 'EaED','EaEDa','ref T_a=0.2'};
labels = cell(8,1);
labels{1} = 'uncoded';
labels{2} = 'BDD';
labels{3} = 'EaED';


markers = {'o', 's', 'd', 'x', '+', '^','o', 's', 'd', 'x', '+', '^'};
myColors = {'r','k','m','b','g','y','c'};
%recorded from simulation=============================================
Ta=0.2;
idx = 1;
labels{idx+3} = sprintf('Ta = %.2f', Ta);
ber=[4,0.0164733
4.5,0.0106477
5,0.00597613
5.5,0.00288062
6,0.00122197
6.5,0.000482324
7,0.000182984
7.5,6.5002e-05
8,2.17227e-05
8.5,7.17187e-06
9,2.02734e-06
9.5,4.58984e-07
10,1.17188e-07
];
semilogy(ber(:,1),ber(:,2),['--' markers{idx} myColors{idx}], 'LineWidth', 2);

idx=idx+1;
Ta=0.3;
labels{idx+3} = sprintf('Ta = %.2f', Ta);
ber=[4,0.0159195
4.5,0.00980752
5,0.00501149
5.5,0.00207729
6,0.000705016
6.5,0.000208584
7,5.99141e-05
7.5,1.73047e-05
8,4.60938e-06
8.5,1.39062e-06
9,2.94922e-07
9.5,4.29687e-08
];
semilogy(ber(:,1),ber(:,2),['--' markers{idx} myColors{idx}], 'LineWidth', 2);
idx=idx+1;
Ta=0.4;
labels{idx+3} = sprintf('Ta = %.2f', Ta);
ber=[4,0.0155867
4.5,0.00933485
5,0.00452753
5.5,0.00171685
6,0.000496684
6.5,0.00011459
7,2.24375e-05
7.5,4.44141e-06
8,9.90234e-07
8.5,2.51953e-07
];
semilogy(ber(:,1),ber(:,2),['--' markers{idx} myColors{idx}], 'LineWidth', 2);

idx=idx+1;
Ta=0.5;
labels{idx+3} = sprintf('Ta = %.2f', Ta);
ber=[4,0.0154387
4.5,0.00911249
5,0.00431571
5.5,0.00156825
6,0.000422906
6.5,8.45605e-05
7,1.30645e-05
7.5,1.89063e-06
8,1.34766e-07
];
semilogy(ber(:,1),ber(:,2),['--' markers{idx} myColors{idx}], 'LineWidth', 2);
idx=idx+1;
Ta=0.6;
labels{idx+3} = sprintf('Ta = %.2f', Ta);
ber=[4,0.0153979
4.5,0.0090443
5,0.00424688
5.5,0.00152161
6,0.000402361
6.5,7.66523e-05
7,9.93945e-06
7.5,9.375e-07
];
semilogy(ber(:,1),ber(:,2),['--' markers{idx} myColors{idx}], 'LineWidth', 2);




for idx=1:5
    semilogy(EbNo_list, EaEDa_ber_list(idx,:),['-' markers{idx} myColors{idx}], 'LineWidth', 1);
end



legend(labels);