% add path for Sims Tools
%addpath('E:\Working folder\Sims Tools') %paths to sims tools. optimization
%toolbox???
global M_ options_ oo_


switch options_.policy_regime %chosen in mod regime is saved globally


    case 0 % No MAP


        %Compute welfare implied by the solution obtained by fmincon
        load('final_mp_only_param_coop','final_mp_only_param_coop');
        xopt_fmc = final_mp_only_param_coop;

        set_param_value('rho_ib',xopt_fmc(1));
        set_param_value('phi_pie',xopt_fmc(2));

        var_list_={'W_EMU'};
        %info=stoch_simul(var_list_);
        %fhat_fmc=oo_.mean;

        [info, oo_, options_] = stoch_simul(M_, options_, oo_, var_list_); %get decision rules and moments
        fhat_fmc=oo_.mean(strmatch('W_EMU',var_list_,'exact'));

        %Find csminwel solution for the welfare criterion specified in the respective case
        %in welfare_objective_bound_coop
        xopt=xopt_fmc;
        npar=length(xopt);
        [fhat,xhat,ghat,Hhat,itct,fcount,retcodehat] = csminwel('welfare_objective_bound_coop',xopt,eye(npar)*.5,[] ,1e-8,100);

        final_mp_only_param_0_coop = xhat;
        save('final_mp_only_param_0_coop','final_mp_only_param_0_coop');

        %Check whether csminwel solves for a higher point than fmincon
        %Good if positive (?)
        check_fmc = -fhat - fhat_fmc


        %Compute the SS consumption level needed to reach a steady state welfare equal to the mean welfare under the respective policy regime
        set_param_value('rho_ib',final_mp_only_param_0_coop(1));
        set_param_value('phi_pie',final_mp_only_param_0_coop(2));

        %var_list_=char('W_S', 'W_B');
        var_list_={'W_S', 'W_B'};


        %info=stoch_simul(var_list_);
        [info, oo_, options_] = stoch_simul(M_, options_, oo_, var_list_); %get decision rules and moments

        Chat_S = exp(((1-beta_p))*oo_.mean(1));

        Chat_B = exp(((1-beta_e))*oo_.mean(2));

        MP_only_Chat_0_coop = [Chat_S, Chat_B];

        save('MP_only_Chat_0_coop','MP_only_Chat_0_coop');



        %----------------------------------------------------------------------------------%

    case 1 % MAP exists: Time varying capital reqirements



        %Compute welfare implied by the solution obtained by fmincon
        load('active_map_param_coop','active_map_param_coop');
        xopt_fmc = active_map_param_coop;

        set_param_value('rho_vi',xopt_fmc(1));
        set_param_value('phi_vi',xopt_fmc(2));
        set_param_value('rho_ib',xopt_fmc(3));
        set_param_value('phi_pie',xopt_fmc(4));

        var_list_={'W_EMU'};
        %info=stoch_simul(var_list_);
        %fhat_fmc=oo_.mean;

        [info, oo_, options_] = stoch_simul(M_, options_, oo_, var_list_); %get decision rules and moments
        fhat_fmc=oo_.mean(strmatch('W_EMU',var_list_,'exact'));

        %Find csminwel solution for the welfare criterion specified in the respective case
        %in welfare_objective_bound_coop
        xopt=xopt_fmc;
        npar=length(xopt);
        [fhat,xhat,ghat,Hhat,itct,fcount,retcodehat] = csminwel('welfare_objective_bound_coop',xopt,eye(npar)*.5,[] ,1e-8,100);

        active_map_param_1_coop = xhat;
        save('active_map_param_1_coop','active_map_param_1_coop');

        %Check whether csminwel solves for a higher point than fmincon
        check_fmc = -fhat - fhat_fmc


        %Compute the SS consumption level needed to reach a steady state welfare equal to the mean welfare under the respective policy regime
        set_param_value('rho_vi',active_map_param_1_coop(1));
        set_param_value('phi_vi',active_map_param_1_coop(2));
        set_param_value('rho_ib',active_map_param_1_coop(3));
        set_param_value('phi_pie',active_map_param_1_coop(4));

        %var_list_=char('W_S', 'W_B');
        var_list_={'W_S', 'W_B'};


        %info=stoch_simul(var_list_);
        [info, oo_, options_] = stoch_simul(M_, options_, oo_, var_list_); %get decision rules and moments

        Chat_S = exp(((1-beta_p))*oo_.mean(1));

        Chat_B = exp(((1-beta_e))*oo_.mean(2));

        MAP_Chat_1_coop = [Chat_S, Chat_B];

        save('MAP_Chat_1_coop','MAP_Chat_1_coop');

    case 2 % LTV cap

        %Compute welfare implied by the solution obtained by fmincon
        load('ltv_cap_coop','ltv_cap_coop');
        xopt_fmc = ltv_cap_coop;

        set_param_value('rho_m_e',xopt_fmc(1));
        set_param_value('rho_ib',xopt_fmc(2));
        set_param_value('phi_pie',xopt_fmc(3));

        var_list_={'W_EMU'};
        %info=stoch_simul(var_list_);
        %fhat_fmc=oo_.mean;

        [info, oo_, options_] = stoch_simul(M_, options_, oo_, var_list_); %get decision rules and moments
        fhat_fmc=oo_.mean(strmatch('W_EMU',var_list_,'exact'));

        %Find csminwel solution for the welfare criterion specified in the respective case
        %in welfare_objective_bound_coop
        xopt=xopt_fmc;
        npar=length(xopt);
        [fhat,xhat,ghat,Hhat,itct,fcount,retcodehat] = csminwel('welfare_objective_bound_coop',xopt,eye(npar)*.5,[] ,1e-8,100);

        ltv_cap_2_coop = xhat;
        save('ltv_cap_2_coop','ltv_cap_2_coop');

        %Check whether csminwel solves for a higher point than fmincon
        check_fmc = -fhat - fhat_fmc


        %Compute the SS consumption level needed to reach a steady state welfare equal to the mean welfare under the respective policy regime
        set_param_value('rho_m_e',ltv_cap_2_coop(1));
        set_param_value('rho_ib',ltv_cap_2_coop(2));
        set_param_value('phi_pie',ltv_cap_2_coop(3));


        %var_list_=char('W_S', 'W_B');
        var_list_={'W_S', 'W_B'};


        %info=stoch_simul(var_list_);
        [info, oo_, options_] = stoch_simul(M_, options_, oo_, var_list_); %get decision rules and moments

        Chat_S = exp(((1-beta_p))*oo_.mean(1));

        Chat_B = exp(((1-beta_e))*oo_.mean(2));

        LTV_Chat_2_coop = [Chat_S, Chat_B];

        save('LTV_Chat_2_coop','LTV_Chat_2_coop');


    case 3 % Weights

        %Compute welfare implied by the solution obtained by fmincon
        load('weights_param_coop','weights_param_coop');
        xopt_fmc = weights_param_coop;

        set_param_value('rho_w_bu',xopt_fmc(1));
        set_param_value('chi_w_bu',xopt_fmc(2));
        set_param_value('rho_w_b',xopt_fmc(3));
        set_param_value('chi_w_b',xopt_fmc(4));
        set_param_value('rho_ib',xopt_fmc(5));
        set_param_value('phi_pie',xopt_fmc(6));

        var_list_={'W_EMU'};
        %info=stoch_simul(var_list_);
        %fhat_fmc=oo_.mean;

        [info, oo_, options_] = stoch_simul(M_, options_, oo_, var_list_); %get decision rules and moments
        fhat_fmc=oo_.mean(strmatch('W_EMU',var_list_,'exact'));

        %Find csminwel solution for the welfare criterion specified in the respective case
        %in welfare_objective_bound_coop
        xopt=xopt_fmc;
        npar=length(xopt);
        [fhat,xhat,ghat,Hhat,itct,fcount,retcodehat] = csminwel('welfare_objective_bound_coop',xopt,eye(npar)*.5,[] ,1e-8,100);

        weights_3_coop = xhat;
        save('weights_3_coop','weights_3_coop');

        %Check whether csminwel solves for a higher point than fmincon
        check_fmc = -fhat - fhat_fmc


        %Compute the SS consumption level needed to reach a steady state welfare equal to the mean welfare under the respective policy regime
        set_param_value('rho_w_bu',weights_3_coop(1));
        set_param_value('chi_w_bu',weights_3_coop(2));
        set_param_value('rho_w_b',weights_3_coop(3));
        set_param_value('chi_w_b',weights_3_coop(4));
        set_param_value('rho_ib',weights_3_coop(5));
        set_param_value('phi_pie',weights_3_coop(6));


        %var_list_=char('W_S', 'W_B');
        var_list_={'W_S', 'W_B'};


        %info=stoch_simul(var_list_);
        [info, oo_, options_] = stoch_simul(M_, options_, oo_, var_list_); %get decision rules and moments

        Chat_S = exp(((1-beta_p))*oo_.mean(1));

        Chat_B = exp(((1-beta_e))*oo_.mean(2));

        weights_Chat_3_coop = [Chat_S, Chat_B];

        save('weights_Chat_3_coop','weights_Chat_3_coop');

    case 4 % LAW: augmented with borrowings

        %Compute welfare implied by the solution obtained by fmincon
        load('law_borrowing_coop','law_borrowing_coop');
        xopt_fmc = law_borrowing_coop;

        set_param_value('rho_ib',xopt_fmc(1));
        set_param_value('phi_pie',xopt_fmc(2));
        set_param_value('phi_B',xopt_fmc(3));


        var_list_={'W_EMU'};
        %info=stoch_simul(var_list_);
        %fhat_fmc=oo_.mean;

        [info, oo_, options_] = stoch_simul(M_, options_, oo_, var_list_); %get decision rules and moments
        fhat_fmc=oo_.mean(strmatch('W_EMU',var_list_,'exact'));

        %Find csminwel solution for the welfare criterion specified in the respective case
        %in welfare_objective_bound_coop
        xopt=xopt_fmc;
        npar=length(xopt);
        [fhat,xhat,ghat,Hhat,itct,fcount,retcodehat] = csminwel('welfare_objective_bound_coop',xopt,eye(npar)*.5,[] ,1e-8,100);

        law_borr_4_coop = xhat;
        save('law_borr_4_coop','law_borr_4_coop');

        %Check whether csminwel solves for a higher point than fmincon
        check_fmc = -fhat - fhat_fmc


        %Compute the SS consumption level needed to reach a steady state welfare equal to the mean welfare under the respective policy regime
        set_param_value('rho_ib',law_borr_4_coop(1));
        set_param_value('phi_pie',law_borr_4_coop(2));
        set_param_value('phi_B',law_borr_4_coop(3));

        %var_list_=char('W_S', 'W_B');
        var_list_={'W_S', 'W_B'};


        %info=stoch_simul(var_list_);
        [info, oo_, options_] = stoch_simul(M_, options_, oo_, var_list_); %get decision rules and moments

        Chat_S = exp(((1-beta_p))*oo_.mean(1));

        Chat_B = exp(((1-beta_e))*oo_.mean(2));

        law_borr_Chat_4_coop = [Chat_S, Chat_B];

        save('law_borr_Chat_4_coop','law_borr_Chat_4_coop');


    case 5 % LAW: augmented capital prices

        %Compute welfare implied by the solution obtained by fmincon
        load('law_asset_prices_coop','law_asset_prices_coop');
        xopt_fmc = law_asset_prices_coop;

        set_param_value('rho_ib',xopt_fmc(1));
        set_param_value('phi_pie',xopt_fmc(2));
        set_param_value('phi_AP',xopt_fmc(3));


        var_list_={'W_EMU'};
        %info=stoch_simul(var_list_);
        %fhat_fmc=oo_.mean;

        [info, oo_, options_] = stoch_simul(M_, options_, oo_, var_list_); %get decision rules and moments
        fhat_fmc=oo_.mean(strmatch('W_EMU',var_list_,'exact'));

        %Find csminwel solution for the welfare criterion specified in the respective case
        %in welfare_objective_bound_coop
        xopt=xopt_fmc;
        npar=length(xopt);
        [fhat,xhat,ghat,Hhat,itct,fcount,retcodehat] = csminwel('welfare_objective_bound_coop',xopt,eye(npar)*.5,[] ,1e-8,100);

        law_asset_prices_coop = xhat;
        save('law_asset_prices_coop','law_asset_prices_coop');

        %Check whether csminwel solves for a higher point than fmincon
        check_fmc = -fhat - fhat_fmc


        %Compute the SS consumption level needed to reach a steady state welfare equal to the mean welfare under the respective policy regime
        set_param_value('rho_ib',law_asset_prices_coop(1));
        set_param_value('phi_pie',law_asset_prices_coop(2));
        set_param_value('phi_y',law_asset_prices_coop(3));



        %var_list_=char('W_S', 'W_B');
        var_list_={'W_S', 'W_B'};


        %info=stoch_simul(var_list_);
        [info, oo_, options_] = stoch_simul(M_, options_, oo_, var_list_); %get decision rules and moments

        Chat_S = exp(((1-beta_p))*oo_.mean(1));

        Chat_B = exp(((1-beta_e))*oo_.mean(2));

        law_asset_prices_Chat_5_coop = [Chat_S, Chat_B];

        save('law_asset_prices_Chat_5_coop','law_asset_prices_Chat_5_coop');

end