%%
clc; clear; close all;
load('PGEDataMod.mat')
%%
NumUsers = size(EnerCon,1);
UserType = UserInfo.Type;
MonthType = ['W';repmat('S',4,1);repmat('W',8,1);repmat('S',4,1);repmat('W',8,1);repmat('S',2,1)];
%% User tariffs
PiMaxUser = zeros(NumUsers,3);
PiMaxUser(UserType == 2,1) = 22.16;
PiMaxUser(UserType == 3,1) = 39.22;
%%
% For summer
% 12:00 am (01) to 02:00 pm (56): Off-peak
% 02:00 pm (57) to 04:00 pm (64): Partial-Peak
% 04:00 pm (65) to 09:00 pm (84): Peak
% 09:00 pm (85) to 11:00 pm (92): Partial-peak
% 11:00 pm (93) to 12:00 am (96): Off-peak

% For winter
% 12:00 am (01) to 09:00 am (36): Off-peak
% 09:00 am (37) to 02:00 pm (56): Super off-peak
% 02:00 pm (57) to 04:00 pm (64): Off-Peak
% 04:00 pm (65) to 09:00 pm (84): Peak
% 09:00 pm (85) to 12:00 am (96): Off-peak

% for summer
% 1: peak, 2: partial-peak, 0: for rest
SummerPeakIndentifier = zeros(96,1);
SummerPeakIndentifier(65:84) = 1;
SummerPeakIndentifier([57:64,85:92]) = 2;
% WinterPeakIndentifier = zeros(96,1);
% WinterPeakIndentifier(65:84) = 1;

%% BEV tariff
% 12:00 am (01) to 09:00 am (36): Off-peak
% 09:00 am (37) to 02:00 pm (56): Super off-peak
% 02:00 pm (57) to 04:00 pm (64): Off-peak
% 04:00 pm (65) to 09:00 pm (84): peak
% 09:00 pm (85) to 11:45 pm (96): Off-peak
PiT = zeros(96,1);
PiT(1:36) = 0.189;
PiT(37:56) = 0.162;
PiT(57:64) = 0.189;
PiT(65:84) = 0.38;
PiT(85:96) = 0.189;
PiMaxEV = 1.24;
%%
Delta = 0.25;
XDes = UserInfo.xRel;
YDes = UserInfo.yRel;

PiMaxEVPrime = PiMaxEV/Delta;
SOCMax = 30;
SOCMin = 10;
EnerPerMile = 40/150;
DepCostPerkWh = (2420/12330)*(1/EnerPerMile);
eta_rt = 0.87;
eta_ch = sqrt(eta_rt);
eta_dis = eta_ch;
%%
UserCharLim = (15/4)*ones(NumUsers,1);  
UserCharLim(UserType == 3,1) = (30/4);
EVCharLim = 30/4;

for NumDrivers = 2:5                 % 1:5
InitSOC = 20*ones(NumDrivers,1);

for month = 1:2                 % 1:27
T = 96*DaysInMonth(month);
TimesTermSOC = 96:96:T;
PiTEV = repmat(PiT,DaysInMonth(month),1);
PeakHours = find(repmat(SummerPeakIndentifier,DaysInMonth(month),1)==1);
if MonthType(month) == 'S'
    PartialPeakHours = find(repmat(SummerPeakIndentifier,DaysInMonth(month),1)==2);
    PiMaxUser(UserType == 3,2) = 54.17;
    PiMaxUser(UserType == 3,3) = 11.75;
else
    PartialPeakHours = 1;
    PiMaxUser(UserType == 3,2) = 3.2;
    PiMaxUser(UserType == 3,3) = 0;
end
PiMaxUserPrime = PiMaxUser/Delta;
SelEnerCon = EnerCon(:,MonthStartMet(month):MonthEndMet(month));
PeakEnerCon = max(SelEnerCon,[],2);
%%
logfilename = ['Log_Mixed_NumDrivers',num2str(NumDrivers), ...
               '_Month',num2str(month), '.txt'];

diary(logfilename)
diary on

clear e_ch e_ch_total e_dis e_tr SOC a1 a2 M u x y m Profit
cvx_begin
cvx_solver Gurobi_2
cvx_solver_settings('TimeLimit',5*3600)
variables ModEnerCon(NumUsers,T) e_ch(NumDrivers,T) e_ch_total(1,T) e_dis(NumDrivers,T) e_tr(NumDrivers,T) SOC(NumDrivers,T) a1(NumUsers,NumDrivers,T) a2(NumUsers,NumDrivers,T);
variable M(NumUsers,NumDrivers,T) binary;
variable u(NumUsers,T) binary; 

OrgCost = PiMaxUserPrime(:,1)'*max(SelEnerCon,[],2) + PiMaxUserPrime(:,2)'*max(SelEnerCon(:,PeakHours),[],2) + PiMaxUserPrime(:,3)'*max(SelEnerCon(:,PartialPeakHours),[],2);
CICost = PiMaxUserPrime(:,1)'*max(ModEnerCon,[],2) + PiMaxUserPrime(:,2)'*max(ModEnerCon(:,PeakHours),[],2) + PiMaxUserPrime(:,3)'*max(ModEnerCon(:,PartialPeakHours),[],2);
EVCost = PiMaxEVPrime*max(e_ch_total) + e_ch_total*PiTEV + DepCostPerkWh*sum(sum(e_dis+e_tr));
NewCost = CICost + EVCost;
minimize NewCost - OrgCost
subject to
ModEnerCon == SelEnerCon-(UserCharLim*ones(1,T).*u);
M(:,:,1) == 0;                             % not required but adding for faster convergence
M(:,:,[TimesTermSOC-1,TimesTermSOC]) == 0; % not required but adding for faster convergence
sum(M,1) <= 1;
a1 >= 0;
a1(:,:,2:T) >= M(:,:,1:T-1) + M(:,:,2:T) - ones(NumUsers,NumDrivers,T-1); 
a1(:,:,2:T) <= M(:,:,1:T-1);
a1(:,:,2:T) <= M(:,:,2:T);
sum(a1,2) <= 1;
UserCharLim*ones(1,T).*u <= SelEnerCon;
u(:,2:T) <= reshape(sum(a1(:,:,2:T),2),NumUsers,T-1);
a2 >= 0;
a2 <= a1;
for j = 1:NumDrivers
    reshape(a2(:,j,:),NumUsers,T) >= reshape(a1(:,j,:),NumUsers,T) + u - ones(NumUsers,T);
    reshape(a2(:,j,:),NumUsers,T) <= u;
end
e_dis == reshape(sum(UserCharLim.*ones(NumUsers,NumDrivers,T).*a2,1),NumDrivers,T);
x = reshape(XDes'*reshape(M,NumUsers,NumDrivers*T),NumDrivers,T);
y = reshape(YDes'*reshape(M,NumUsers,NumDrivers*T),NumDrivers,T);
e_tr(:,1) == 0;
e_tr(:,2:T) >= EnerPerMile* (abs(x(:,2:T)-x(:,1:T-1)) + abs(y(:,2:T)-y(:,1:T-1)));
e_ch >= 0; 
m = ones(NumDrivers,T) - reshape(sum(M,1),NumDrivers,T);
e_ch(:,2:T) <= EVCharLim*min(m(:,1:T-1),m(:,2:T));
e_ch_total == sum(e_ch,1);
SOC(:,2:T) == SOC(:,1:T-1) + eta_ch*e_ch(:,2:T) - (1/eta_dis)*e_dis(:,2:T) - e_tr(:,2:T);
SOC <= SOCMax;
SOC >= SOCMin;
SOC(:,1) == InitSOC;
SOC(:,TimesTermSOC) == InitSOC*ones(1,length(TimesTermSOC)); 
cvx_end
diary off

Profit = - cvx_optval

% filename = ['MixedSetup','NumDrivers',num2str(NumDrivers),'Month',num2str(month),'MIP'];
% save(filename)
drawnow
end
end

%%
