function [FeasCheck,CC,e_ch,SOC] = ChargingCost(ucParam,evParam,M,e_dis,e_tr)

    T = ucParam.T;
    NumDrivers = evParam.NumDrivers;
    eta_ch = evParam.eta_ch;
    eta_dis = eta_ch;
    
    m = ones(NumDrivers,T) - reshape(sum(M,1),NumDrivers,T);
    cvx_begin quiet
    cvx_solver Gurobi_2
    variables e_ch(NumDrivers,T) e_ch_total(1,T) SOC(NumDrivers,T);
    CC = ucParam.PiMaxEVPrime*max(e_ch_total) + e_ch_total*ucParam.PiTEV;
    minimize CC
    subject to
    e_ch >= 0; 
    e_ch(:,2:T) <= evParam.EVCharLim*min(m(:,1:T-1),m(:,2:T));
    e_ch_total == sum(e_ch,1);
    SOC(:,2:T) == SOC(:,1:T-1) + eta_ch*e_ch(:,2:T) - (1/eta_dis)*e_dis(:,2:T) - e_tr(:,2:T);
    SOC <= evParam.SOCMax;
    SOC >= evParam.SOCMin;
    SOC(:,1) == evParam.InitSOC;
    SOC(:,evParam.TimesTermSOC) == evParam.InitSOC*ones(1,length(evParam.TimesTermSOC)); 
    cvx_end

    if string(cvx_status) == 'Solved'
        FeasCheck = 1;
    else
        FeasCheck = 0;
    end



end