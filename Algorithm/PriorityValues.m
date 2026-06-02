function [P,s_max] = PriorityValues(ucParam,ciParam,evParam,M,u,s_max,P,SelUsers)

    UserCharLim = ciParam.UserCharLim;
    PiMaxUserPrime = ucParam.PiMaxUserPrime;
    ModEnerCon = ciParam.EnerCon - UserCharLim.*u;
    ModDC = PiMaxUserPrime(:,1).*max(ModEnerCon,[],2) + PiMaxUserPrime(:,2).*max(ModEnerCon(:,ucParam.PeakHours),[],2) + PiMaxUserPrime(:,3).*max(ModEnerCon(:,ucParam.PartialPeakHours),[],2);

    for user = SelUsers
        if s_max(user) >= 1
        for s = 1:s_max(user)
            [DC_plus, ~] = NewDC(user,s,ciParam,ucParam,evParam,u,M);
            P(user,s) = (ModDC(user)-DC_plus)/s;
        end
        if sum(P(user,:)) == 0
            s_max(user) = 0;
        end
        end
    end

end