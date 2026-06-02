function [j_star] = MostProfitableEV(ucParam,ciParam,evParam,M,u,i_star,T_ser,T_tr,Jf)

    EnerCon = ciParam.EnerCon;
    UserCharLim = ciParam.UserCharLim;
    ModEnerCon = EnerCon-UserCharLim.*u;
    PiMaxUserPrime = ucParam.PiMaxUserPrime;
    ModDC = PiMaxUserPrime(:,1).*max(ModEnerCon,[],2) + PiMaxUserPrime(:,2).*max(ModEnerCon(:,ucParam.PeakHours),[],2) + PiMaxUserPrime(:,3).*max(ModEnerCon(:,ucParam.PartialPeakHours),[],2);
    DepCostPerkWh = evParam.DepCostPerkWh;

    u_plus = u;
    u_plus(i_star,T_ser) = 1;
    y_tilde_plus = ModEnerCon(i_star,:) - UserCharLim(i_star)*u_plus(i_star,:);
    MUS = ModDC(i_star) - (ucParam.PiMaxUserPrime(i_star,1)*max(y_tilde_plus) + ucParam.PiMaxUserPrime(i_star,2)*max(y_tilde_plus(ucParam.PeakHours)) + ucParam.PiMaxUserPrime(i_star,3)*max(y_tilde_plus(ucParam.PartialPeakHours)));
    [e_dis,e_tr] = DischargeTransitEnergy(ucParam,ciParam,evParam,M,u);
    [~,CC] = ChargingCost(ucParam,evParam,M,e_dis,e_tr);
    OC = DepCostPerkWh*(sum(e_dis(:) + e_tr(:)))+ CC;
    
    % MOC = zeros(length(Jf),1);
    % for j = 1:length(Jf)
    %     M_plus = M;
    %     M_plus(i_star,Jf(j),[T_ser,T_tr]) = 1; 
    %     [e_dis_plus,e_tr_plus] = DischargeTransitEnergy(ucParam,ciParam,evParam,M_plus,u_plus);
    %     [FeasCheck,CC_plus] = ChargingCost(ucParam,evParam,M_plus,e_dis_plus,e_tr_plus);
    %     if FeasCheck == 1
    %         OC_plus = DepCostPerkWh*(sum(e_dis_plus(:) + e_tr_plus(:)))+ CC_plus;
    %         MOC(j,1) = OC_plus-OC;
    %     else
    %         MOC(j,1) = inf;
    %     end 
    % end
    % [MOC_star,j_star_ind] = min(MOC);
    % j_star = Jf(j_star_ind);

    [~,j_star_ind] = min(sum(e_dis(Jf',:) + e_tr(Jf',:),2));
    j_star = Jf(j_star_ind);
    M_plus = M;
    M_plus(i_star,j_star,[T_ser,T_tr]) = 1; 
    [e_dis_plus,e_tr_plus] = DischargeTransitEnergy(ucParam,ciParam,evParam,M_plus,u_plus);
    [FeasCheck,CC_plus] = ChargingCost(ucParam,evParam,M_plus,e_dis_plus,e_tr_plus);
    if FeasCheck == 1
        OC_plus = DepCostPerkWh*(sum(e_dis_plus(:) + e_tr_plus(:)))+ CC_plus;
        MOC_star = OC_plus-OC;
    else
        MOC_star = inf;
    end 

    if MOC_star > MUS
        j_star = 0;
    end

end