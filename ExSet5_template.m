%% Import and select data & Graph Cumulative Returns to Assets
clear;clc;
cd ''; 
%Read Portfolio Asset Classes
tblAssets = readtable('Econ 4360 DataSets.xlsx','sheet','Your Sheet'); 
date = datetime(tblAssets.Date); 
M = [tblAssets.yourIndex .... ]; 
T = size(M,1);
mean_return = 
AssetNames = {'Asset class names in single quotes'}; 
% Plot Cumulative Factor Returns over the sample
temp =  cumulative returns calc here
figure(1);
plot(date,temp); title('Cumulative Value of $1 Investment'); legend('indexes in quotes'); grid on;
%% 2.a.	Set Portfolio Optimization Constraints & Options 
nAssets = numel(mean_return); r = 'target return';     % number of assets and desired return
Aeq = ones(1,nAssets); beq = 1;              % equality Aeq*x = beq
Aineq = -mean_return; bineq = -r;            % inequality Aineq*x <= bineq
lb = zeros(nAssets,1); ub = ones(nAssets,1); % bounds lb <= x <= ub
c = zeros(nAssets,1);                        % objective has no linear term; set it to zero
% Select Options for Quadprog
options = optimset('Algorithm','interior-point-convex');
options = optimset(options,'Display','off','TolFun',1e-10);

%% 2.b.	Estimate Vcov each month (begin at 36 months) and update monthly

for t = 36:T
    V = cov(M(t-35:t,:)); 
    Vcov(:,:,t) = V; Vcond(t) = cond(V); 
end

%% Optimize three base and three stress cases
% Identify the Stress Risk Matrix (Vcov) for July 2010
dateStr = datestr(date); 
n = strmatch('30-Jul-2010',dateStr,'exact'); 
% Base case:  
Vbase = cov(M);  % Average through whole sample
wmvBase = quadprog(Vbase,c,[],[],Aeq,beq,[],[],[],options); 
wBase = quadprog(Vbase,c,Aineq,bineq,Aeq,beq,lb,ub,[],options);
wrpBase = OptRP(Vbase); 
% Stress case:  
Vstress = Vcov(:,:,n);  % n = location for '30-Jul-2010'
wmvStress = quadprog(Vstress,c,[],[],Aeq,beq,[],[],[],options); 
wStress = quadprog(Vstress,c,Aineq,bineq,Aeq,beq,lb,ub,[],options);
wrpStress = OptRP(Vstress); 
% Expected returns
Rmv_base =  
R_base = 
R_rp_base = 
Rmv_base_stress =
R_base_stress =
R_rp_stress = 

% 2.c. Table results
PortWts = table(wmvBase,wmvStress,wBase,wStress,wrpBase,wrpStress,'RowNames',{'index names in single quotes'})
PortRets = table([Rmv_base;R_base;R_rp_base],[Rmv_base_stress;R_base_stress;R_rp_stress],'RowNames',{'mv','targRet','riskparity'},'VariableNames',{'Base','Stress'})

% Risk Attribution Cases
riskBase = round(sqrt(wBase'*(12*Vbase)*wBase),3);
riskStress = 
riskParityBase = 
riskParityStress = 
riskScenarios = table(riskBase,riskStress,riskParityBase,riskParityStress, ...
    'VariableNames',{'MaxSharpeBase','MaxSharpeStress','RiskParityBase','RiskParityStress'})
% look at how well Risk Parity holds up under stress
MCR = table(round(12*Vbase*wBase/riskBase,3),...),...
    'VariableNames',{'MCRBase','MCRStress','MCRRPBase','MCRRPStress'},'RowNames',{'indexes in single quotes'})
RiskBudgets =  table(round((12*Vbase*wBase/riskBase).*wBase,3),...),...
    'VariableNames',{'RskBudgBase','RskBudgStress','RskBudgRPBase','RskBudgRPStress'},'RowNames',{'indexes in single quotes'})

%% 3. 99% VaR
VaR_01 = ((1+[R_base, R_base_stress, R_rp_base, R_rp_stress]).^12 -1) + norminv('arguments')*[table2array(riskScenarios)]
%% 4. Create Portfolio Returns

% R_bt = returns for base case targeted return;  
% R_st = returns for stress case targeted return; 
% R_brp = returns for base risk parity; 
% R_srp = returns for stress risk parity; 

R_bt = M*wBase; R_st = M*wStress; 
R_brp = M*wrpBase; R_srp = M*wrpStress; 

%% the base (stress) portfolio will lose at least 9.25% (15.30%) of its value 1% of the time

%choose a portfolio 
Y = R_st; 

%% Read in Factors
Factors = readtable('Econ 4360 DataSets.xlsx','sheet','BLK_AQR'); 
Factors = Factors(start:end,:); % note the START date!!! 
dateFactors = datetime(Factors.Date);
X = [Factors.BAB, Factors.MKT-Factors.RF, Factors.SMB, Factors.HML, Factors.UMD]; % Switch to BlackRock Factors
[T,k]= size(X); 
FactorNames = {'BAB','MKT','SMB','HML','MOM'}; % Switch to BlackRock Factors

%% Estimate Betas in 36-month windows
b = zeros(T, ); % create an empty matrix to store betas N rows x cols=(intercept and k factors). 
A = zeros(T, ); % create an empty matrix to store factor attributionex N rows x  cols
% Note: attribution tracks the residual (idiosyncratic), the intercept, and
% the factors.  
for t = 36:size(X,1) 
    LM = fitlm(X(t-35:t,:),Y(t-35:t)); 
    b(t,:) = (LM.Coefficients.Estimate)';
    A(t,:) = [Y(t)-b(t,1)-X(t,:)*b(t,2:end)', b(t,1), b(t,2:end).*X(t,:)]; 
end 
% You should verify the attribution (should sum to Y) 
% Replace the AQR factors names below with BlackRock
Attribution = table(dateFactors,A(:,1), A(:,2),A(:,3),A(:,4),A(:,5),A(:,6),A(:,7),...
    'VariableNames',{'Date','Idiosyn','Alpha','BAB','MKT','SMB','HML','MOM'}); 
%% Plots (Edit all below to accomodate BlackRock Factors)
% Beta exposures first 
figure(2);
subplot(3,2,1); plot(dateFactors(36:end),b(36:end,2)); recessionplot; title('BAB beta'); grid on; 
subplot(3,2,2); plot(dateFactors(36:end),b(36:end,3)); recessionplot; title('MKT beta'); grid on; 
subplot(3,2,3); plot(dateFactors(36:end),b(36:end,4)); recessionplot; title('SMB beta'); grid on; 
subplot(3,2,4); plot(dateFactors(36:end),b(36:end,5)); recessionplot; title('HML beta'); grid on; 
subplot(3,2,5); plot(dateFactors(36:end),b(36:end,6)); recessionplot; title('MOM beta'); grid on; 

% Return Attribution 
figure(3);
subplot(4,2,1); plot(dateFactors(36:end),A(36:end,1)); recessionplot; title('Idiosyn Return Attr'); grid on; 
subplot(4,2,2); plot(dateFactors(36:end),A(36:end,2)); recessionplot; title('Alpha Attr'); grid on; 
subplot(4,2,3); plot(dateFactors(36:end),A(36:end,3)); recessionplot; title('BAB Attr'); grid on; 
subplot(4,2,4); plot(dateFactors(36:end),A(36:end,4)); recessionplot; title('MKT  Attr'); grid on; 
subplot(4,2,5); plot(dateFactors(36:end),A(36:end,5)); recessionplot; title('SMB Attr'); grid on; 
subplot(4,2,6); plot(dateFactors(36:end),A(36:end,6)); recessionplot; title('HML Attr'); grid on; 
subplot(4,2,7); plot(dateFactors(36:end),A(36:end,7)); recessionplot; title('MOM Attr'); grid on; 

%% 5-6 Remarks and Interpretation

% Open a blank Excel Sheet and copy all the tables and plots to Excel.  Reformat and edit ...
%     Change format in Excel, Word to Garamond 11.  Number the tables and add titles ...
%     Copy plots and tables to Word and provide commentary as requested in
%     parts 5-6.  When finished, save as PDF and post the PDF to Canvas.  
