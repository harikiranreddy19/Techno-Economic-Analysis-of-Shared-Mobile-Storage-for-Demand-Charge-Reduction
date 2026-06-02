%%
clc; clear; close all;
load('PGEDataMod.mat')
%%
NumUsers = size(EnerCon,1);
PeakEnerCon = zeros(NumUsers,NumMonths);
for month = 1:NumMonths
    PeakEnerCon(:,month) = max(EnerCon(:,MonthStartMet(month):MonthEndMet(month)),[],2);
end
%%
PeakPowCon = 4*PeakEnerCon;
MeanPeakPowCon = mean(PeakPowCon,2);
UserType = zeros(NumUsers,1);
%%
for user = 1:NumUsers
    if MeanPeakPowCon(user) <= 350 
        UserType(user) = 2;
    else
        UserType(user) = 3;
    end
end
%%
MonthType = [repmat('S',5,1);repmat('W',6,1);repmat('S',6,1);repmat('W',6,1);repmat('S',4,1)];
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

%%
ucParam.Delta = 0.25;

% c&i user parameters
ciParam.NumUsers = NumUsers;
ciParam.XDes = UserInfo.xRel(1:NumUsers);
ciParam.YDes = UserInfo.yRel(1:NumUsers);
ciParam.S = 5;
ciParam.UserCharLim = (15/4)*ones(NumUsers,1);  
ciParam.UserCharLim(UserType == 3,1) = (30/4);

% ev driver parameters
evParam.SOCMax = 30;
evParam.SOCMin = 10;
evParam.EnerPerMile = 40/150;
evParam.DepCostPerkWh = (2420/12330)*(1/evParam.EnerPerMile);
evParam.eta_rt = 0.87;
evParam.eta_ch = sqrt(evParam.eta_rt);
evParam.eta_dis = evParam.eta_ch;
evParam.EVCharLim = 30/4;

for month = 1:12  

for NumDrivers = 1:5                 % 1:5
evParam.NumDrivers = NumDrivers;
evParam.InitSOC = 20*ones(NumDrivers,1);

ucParam.T = 96*DaysInMonth(month);
evParam.TimesTermSOC = [1,95:96:ucParam.T-1,96:96:ucParam.T];
ucParam.PeakHours = find(repmat(SummerPeakIndentifier,DaysInMonth(month),1)==1);
if MonthType(month) == 'S'
    ucParam.PartialPeakHours = find(repmat(SummerPeakIndentifier,DaysInMonth(month),1)==2);
    PiMaxUser(UserType == 2,1) = 11.41;
    PiMaxUser(UserType == 3,1) = 14.42;
    % PiMaxUser(UserType == 3,2) = 0;
    % PiMaxUser(UserType == 3,3) = 0;
    PiT = 0.115*ones(96,1);
    ucParam.PiTEV = repmat(PiT,DaysInMonth(month),1);
    PiMaxEV = 11.41;
    ucParam.PiMaxEVPrime = PiMaxEV/ucParam.Delta;
else
    ucParam.PartialPeakHours = 1;
    PiMaxUser(UserType == 2,1) = 21.63;
    PiMaxUser(UserType == 3,1) = 17.10;
    % PiMaxUser(UserType == 3,2) = 0;
    % PiMaxUser(UserType == 3,3) = 0;
    PiT = 0.124*ones(96,1);
    ucParam.PiTEV = repmat(PiT,DaysInMonth(month),1);
    PiMaxEV = 21.63;
    ucParam.PiMaxEVPrime = PiMaxEV/ucParam.Delta;
end
ucParam.PiMaxUser = PiMaxUser(1:NumUsers,:);
ucParam.PiMaxUserPrime = PiMaxUser(1:NumUsers,:)/ucParam.Delta;
ciParam.EnerCon = EnerCon(1:NumUsers,MonthStartMet(month):MonthEndMet(month));
ciParam.PeakEnerCon = max(ciParam.EnerCon,[],2);
ciParam.OrgDC = ucParam.PiMaxUserPrime(:,1).*max(ciParam.EnerCon,[],2) + ucParam.PiMaxUserPrime(:,2).*max(ciParam.EnerCon(:,ucParam.PeakHours),[],2) + ucParam.PiMaxUserPrime(:,3).*max(ciParam.EnerCon(:,ucParam.PartialPeakHours),[],2);
M = zeros(ciParam.NumUsers, evParam.NumDrivers, ucParam.T);
u = zeros(ciParam.NumUsers,ucParam.T);
s_max = ciParam.S*ones(ciParam.NumUsers,1);
P = zeros(ciParam.NumUsers,ciParam.S);
P = PriorityValues(ucParam,ciParam,evParam,M,u,s_max,P,1:ciParam.NumUsers);
while sum(round(P(:))) ~= 0 
    [i_star,s_star,T_ser,T_tr] = HighestPriorityUserAndServices(ciParam,u,M,P,evParam,ucParam);
    Jf = FeasibleEV(ucParam,ciParam,evParam,M,i_star,T_ser,T_tr);
    if ~isempty(Jf)
        j_star = MostProfitableEV(ucParam,ciParam,evParam,M,u,i_star,T_ser,T_tr,Jf);
        if j_star ~= 0
            M(i_star,j_star,[T_ser,T_tr]) = 1;
            u(i_star,T_ser) = 1;
            [P,s_max] = PriorityValues(ucParam,ciParam,evParam,M,u,s_max,P,i_star);            
        else
            P(i_star,s_star) = 0; 
        end
    else
        P(i_star,s_star:ciParam.S) = 0;
        s_max(i_star) = (s_star-1);
        P = PriorityValues(ucParam,ciParam,evParam,M,u,s_max,P,i_star);        
    end
end
%%
OrgCost = sum(ciParam.OrgDC);
ModEnerCon = ciParam.EnerCon - ciParam.UserCharLim.*u;
PiMaxUserPrime = PiMaxUser/ucParam.Delta;

CICost = PiMaxUserPrime(1:NumUsers,1)'*max(ModEnerCon,[],2) + PiMaxUserPrime(1:NumUsers,2)'*max(ModEnerCon(:,ucParam.PeakHours),[],2) + PiMaxUserPrime(1:NumUsers,3)'*max(ModEnerCon(:,ucParam.PartialPeakHours),[],2);
[e_dis,e_tr,x,y] = DischargeTransitEnergy(ucParam,ciParam,evParam,M,u);
[~,CC,e_ch,SOC] = ChargingCost(ucParam,evParam,M,e_dis,e_tr);
EVCost = evParam.DepCostPerkWh*(sum(e_dis(:) + e_tr(:)))+ CC;
NewCost = CICost + EVCost;

Profit = OrgCost-NewCost;
filename = ['PSEMixedSetup','NumDrivers',num2str(NumDrivers),'Month',num2str(month),'AlgoNew'];
save(filename)
drawnow
end
end