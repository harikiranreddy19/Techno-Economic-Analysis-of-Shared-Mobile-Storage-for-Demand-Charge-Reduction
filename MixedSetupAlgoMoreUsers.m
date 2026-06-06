%%
clc; clear; close all;
load('PGEModMoreUsers');
%%
MaxNumUsers = size(EnerCon,1);
UserType = UserInfo.Type;
MonthType = ['W';repmat('S',4,1);repmat('W',8,1);repmat('S',4,1);repmat('W',8,1);repmat('S',2,1)];
%% User tariffs
PiMaxUser = zeros(MaxNumUsers,3);
PiMaxUser(UserType == 2,1) = 22.16;
PiMaxUser(UserType == 3,1) = 39.22;

SummerPeakIndentifier = zeros(96,1);
SummerPeakIndentifier(65:84) = 1;
SummerPeakIndentifier([57:64,85:92]) = 2;

PiT = zeros(96,1);
PiT(1:36) = 0.189;
PiT(37:56) = 0.162;
PiT(57:64) = 0.189;
PiT(65:84) = 0.38;
PiT(85:96) = 0.189;
PiMaxEV = 1.24;
%%
ucParam.Delta = 0.25;
ucParam.PiMaxEVPrime = PiMaxEV/ucParam.Delta;
ucParam.PiT = PiT;

% ev driver parameters
evParam.SOCMax = 30;
evParam.SOCMin = 10;
evParam.EnerPerMile = 40/150;
evParam.DepCostPerkWh = (2420/12330)*(1/evParam.EnerPerMile);
evParam.eta_rt = 0.87;
evParam.eta_ch = sqrt(evParam.eta_rt);
evParam.eta_dis = evParam.eta_ch;
evParam.EVCharLim = 30/4;

UserCount = [50,100,136,150,200,250,300,350,400,450,470];
RunTimes = zeros(length(UserCount),2);
for UC = 1:length(UserCount) 

NumUsers = UserCount(UC);

% c&i user parameters
ciParam.NumUsers = NumUsers;
ciParam.XDes = UserInfo.xRel(1:NumUsers);
ciParam.YDes = UserInfo.yRel(1:NumUsers);
ciParam.S = 5;
ciParam.UserCharLim = (15/4)*ones(NumUsers,1);  
ciParam.UserCharLim(UserType(1:NumUsers) == 3,1) = (30/4);

for month = 1:2  

NumDrivers = 3;
evParam.NumDrivers = NumDrivers;
evParam.InitSOC = 20*ones(NumDrivers,1);

ucParam.T = 96*DaysInMonth(month);
evParam.TimesTermSOC = [1,95:96:ucParam.T-1,96:96:ucParam.T];
ucParam.PiTEV = repmat(PiT,DaysInMonth(month),1);
ucParam.PeakHours = find(repmat(SummerPeakIndentifier,DaysInMonth(month),1)==1);
if MonthType(month) == 'S'
    ucParam.PartialPeakHours = find(repmat(SummerPeakIndentifier,DaysInMonth(month),1)==2);
    PiMaxUser(UserType == 3,2) = 54.17;
    PiMaxUser(UserType == 3,3) = 11.75;
else
    ucParam.PartialPeakHours = 1;
    PiMaxUser(UserType == 3,2) = 3.2;
    PiMaxUser(UserType == 3,3) = 0;
end
ucParam.PiMaxUser = PiMaxUser(1:NumUsers,:);
ucParam.PiMaxUserPrime = PiMaxUser(1:NumUsers,:)/ucParam.Delta;
ciParam.EnerCon = EnerCon(1:NumUsers,MonthStartMet(month):MonthEndMet(month));
ciParam.PeakEnerCon = max(ciParam.EnerCon,[],2);
ciParam.OrgDC = ucParam.PiMaxUserPrime(:,1).*max(ciParam.EnerCon,[],2) + ucParam.PiMaxUserPrime(:,2).*max(ciParam.EnerCon(:,ucParam.PeakHours),[],2) + ucParam.PiMaxUserPrime(:,3).*max(ciParam.EnerCon(:,ucParam.PartialPeakHours),[],2);

tic
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


RunTimes(UC,month) = toc;
save('MixedSetupMoreUsersResults','RunTimes');

end

end