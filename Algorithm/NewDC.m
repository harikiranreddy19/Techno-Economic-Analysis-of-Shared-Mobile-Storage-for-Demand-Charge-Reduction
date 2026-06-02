function [DC_plus,T_ser] = NewDC(user,s,ciParam,ucParam,evParam,u,M)

    T = ucParam.T;
    UserCharLim = ciParam.UserCharLim;
    EnerCon = ciParam.EnerCon(user,:);
    ModEnerCon = (EnerCon - UserCharLim(user)*u(user,:))';
    PiMaxUserPrime = ucParam.PiMaxUserPrime;
    t1 = ucParam.PeakHours;
    t2 = ucParam.PartialPeakHours;
    LowestDemandPossible = max([min([max(EnerCon(1,:)),max(EnerCon(1,t1)),max(EnerCon(1,t2))]) - UserCharLim(user),UserCharLim(user)]);
    t3 = [find(reshape(sum(M,[1,2]),1,[]) == evParam.NumDrivers),evParam.TimesTermSOC,find(u(user,:)==1),find(ModEnerCon(user,:)<=LowestDemandPossible)];

    cvx_begin quiet
    cvx_solver Gurobi_2
    variables y_tilde(T,1)
    variable u_plus(T,1) binary
    DC_plus = PiMaxUserPrime(user,1)*max(y_tilde) + PiMaxUserPrime(user,2)*max(y_tilde(t1)) + PiMaxUserPrime(user,3)*max(y_tilde(t2));
    minimize DC_plus
    subject to
    y_tilde == ModEnerCon-UserCharLim(user)*u_plus;
    sum(u_plus) <= s;
    u_plus(t3) == 0;
    cvx_end
    T_ser = find(round(u_plus)==1);
    

end