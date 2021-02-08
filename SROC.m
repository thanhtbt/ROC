clc;clear;close all;

%% Load cost/benefit points (Sensitivity and Specificity)
load SEN
load SPE

%% Probit
SEN_log = [];
a = 3.8;
b = 0.24;
N = length(SEN);
x = linspace(0,1,N);

SEN_pro = [];
SEN_pro_AUC = [];
for ii = 1 : N
    r = 1/2 + 1/2*erf( (1+b)/(1-b)*erfinv(1-2*SPE(ii)) + a*sqrt(pi)/(4-4*b) );
    SEN_pro = [SEN_pro; r];
    
    r = 1/2 + 1/2*erf( (1+b)/(1-b)*erfinv(1-2*x(ii)) + a*sqrt(pi)/(4-4*b) );
    SEN_pro_AUC = [SEN_pro_AUC; r];
end


%%  Logit
a = 3.9;
b = 0.1;
SEN_log_AUC = [];
delta = 0.0001;
for ii = 1 : N
    TS      = exp(a/(1-b)) * ( (1-SPE(ii))/SPE(ii) )^( (1+b)/(1-b) );
    MS      = 1 + exp(a/(1-b)) * ( (1-SPE(ii))/SPE(ii) )^( (1+b)/(1-b) );
    SEN_log = [SEN_log; TS/MS];
    
    TS      = exp(a/(1-b)) * ( (1-x(ii))/(x(ii)+delta) )^( (1+b)/(1-b) );
    MS      = 1 + exp(a/(1-b)) * ( (1-x(ii))/(x(ii)+delta) )^( (1+b)/(1-b) );    
    SEN_log_AUC = [SEN_log_AUC; TS/MS];
    
end

%% ERROR
error_log = abs(SEN -SEN_log);
error_log = sum(error_log)/N

error_pro = abs(SEN -SEN_pro);
error_pro = sum(error_pro)/N

%% AUC
S_log = 0;
for ii = 1 : N-1
    S_log = S_log + 1/2*(SEN_log_AUC(ii) + SEN_log_AUC(ii+1))*(x(ii+1)-x(ii));
end
S_log

S_pro = 0;
for ii = 1 : N-1
    S_pro = S_pro + 1/2*(SEN_pro_AUC(ii) + SEN_pro_AUC(ii+1))*(x(ii+1)-x(ii));
end
S_pro

%% PLOT
figure; 
plot(1 - SPE, SEN,'o');
hold on; plot(1 - SPE',SEN_pro','b');
hold on; plot(1 - SPE',SEN_log','r');

legend('Cost/Benefit Point','Gaussian-based ROC','Logistic-based ROC');
xlabel('False Positive Rate');
ylabel('True Positive Rate');
grid on;


