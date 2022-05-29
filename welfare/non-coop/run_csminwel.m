% add path for Sims Tools
%addpath('E:\Working folder\Sims Tools') %paths to sims tools. optimization
%toolbox???
global M_ options_ oo_


switch options_.policy_regime %chosen in mod regime is saved globally


    case 0 % No MAP


        %Compute welfare implied by the solution obtained by fmincon
        load('final_mp_only_param','final_mp_only_param');
        xopt_fmc = final_mp_only_param;

        set_param_value('rho_ib',xopt_fmc(1));
        set_param_value('phi_pie',xopt_fmc(2));
        %set_param_value('phi_y',xopt_fmc(3));

        var_list_={'W_EMU'};
        %info=stoch_simul(var_list_);
        %fhat_fmc=oo_.mean;

        [info, oo_, options_] = stoch_simul(M_, options_, oo_, var_list_); %get decision rules and moments
        fhat_fmc=oo_.mean(strmatch('W_EMU',var_list_,'exact'));

        %Find csminwel solution for the welfare criterion specified in the respective case
        %in welfare_objective_bound
        xopt=xopt_fmc;
        npar=length(xopt);
        [fhat,xhat,ghat,Hhat,itct,fcount,retcodehat] = csminwel('welfare_objective_bound',xopt,eye(npar)*.5,[] ,1e-8,100);

        final_mp_only_param_0 = xhat;
        save('final_mp_only_param_0','final_mp_only_param_0');

        %Check whether csminwel solves for a higher point than fmincon
        %Good if positive (?)
        check_fmc = -fhat - fhat_fmc


        %Compute the SS consumption level needed to reach a steady state welfare equal to the mean welfare under the respective policy regime
        set_param_value('rho_ib',final_mp_only_param_0(1));
        set_param_value('phi_pie',final_mp_only_param_0(2));
        %set_param_value('phi_y',final_mp_only_param_0(3));

        %var_list_=char('W_S', 'W_B');
        var_list_={'W_S', 'W_B'};


        %info=stoch_simul(var_list_);
        [info, oo_, options_] = stoch_simul(M_, options_, oo_, var_list_); %get decision rules and moments

        Chat_S = exp(((1-beta_p))*oo_.mean(1));

        Chat_B = exp(((1-beta_e))*oo_.mean(2));

        MP_only_Chat_0 = [Chat_S, Chat_B];

        save('MP_only_Chat_0','MP_only_Chat_0');



        %----------------------------------------------------------------------------------%

    case 1 % MAP exists: Time varying capital reqirements



        %Compute welfare implied by the solution obtained by fmincon
        load('active_map_param','active_map_param');
        xopt_fmc = active_map_param;

        set_param_value('rho_vi',xopt_fmc(1));
        set_param_value('phi_vi',xopt_fmc(2));

        var_list_={'W_EMU'};
        %info=stoch_simul(var_list_);
        %fhat_fmc=oo_.mean;

        [info, oo_, options_] = stoch_simul(M_, options_, oo_, var_list_); %get decision rules and moments
        fhat_fmc=oo_.mean(strmatch('W_EMU',var_list_,'exact'));

        %Find csminwel solution for the welfare criterion specified in the respective case
        %in welfare_objective_bound
        xopt=xopt_fmc;
        npar=length(xopt);
        [fhat,xhat,ghat,Hhat,itct,fcount,retcodehat] = csminwel('welfare_objective_bound',xopt,eye(npar)*.5,[] ,1e-8,100);

        active_map_param_1 = xhat;
        save('active_map_param_1','active_map_param_1');

        %Check whether csminwel solves for a higher point than fmincon
        check_fmc = -fhat - fhat_fmc


        %Compute the SS consumption level needed to reach a steady state welfare equal to the mean welfare under the respective policy regime
        set_param_value('rho_vi',active_map_param_1(1));
        set_param_value('phi_vi',active_map_param_1(2));

        %var_list_=char('W_S', 'W_B');
        var_list_={'W_S', 'W_B'};


        %info=stoch_simul(var_list_);
        [info, oo_, options_] = stoch_simul(M_, options_, oo_, var_list_); %get decision rules and moments

        Chat_S = exp(((1-beta_p))*oo_.mean(1));

        Chat_B = exp(((1-beta_e))*oo_.mean(2));

        MAP_Chat_1 = [Chat_S, Chat_B];

        save('MAP_Chat_1','MAP_Chat_1');

    case 2 % LTV cap

        %Compute welfare implied by the solution obtained by fmincon
        load('ltv_cap','ltv_cap');
        xopt_fmc = ltv_cap;

        set_param_value('rho_m_e',xopt_fmc(1));

        var_list_={'W_EMU'};
        %info=stoch_simul(var_list_);
        %fhat_fmc=oo_.mean;

        [info, oo_, options_] = stoch_simul(M_, options_, oo_, var_list_); %get decision rules and moments
        fhat_fmc=oo_.mean(strmatch('W_EMU',var_list_,'exact'));

        %Find csminwel solution for the welfare criterion specified in the respective case
        %in welfare_objective_bound
        xopt=xopt_fmc;
        npar=length(xopt);
        [fhat,xhat,ghat,Hhat,itct,fcount,retcodehat] = csminwel('welfare_objective_bound',xopt,eye(npar)*.5,[] ,1e-8,100);

        ltv_cap_2 = xhat;
        save('ltv_cap_2','ltv_cap_2');

        %Check whether csminwel solves for a higher point than fmincon
        check_fmc = -fhat - fhat_fmc


        %Compute the SS consumption level needed to reach a steady state welfare equal to the mean welfare under the respective policy regime
        set_param_value('rho_m_e',ltv_cap_2(1));


        %var_list_=char('W_S', 'W_B');
        var_list_={'W_S', 'W_B'};


        %info=stoch_simul(var_list_);
        [info, oo_, options_] = stoch_simul(M_, options_, oo_, var_list_); %get decision rules and moments

        Chat_S = exp(((1-beta_p))*oo_.mean(1));

        Chat_B = exp(((1-beta_e))*oo_.mean(2));

        LTV_Chat_2 = [Chat_S, Chat_B];

        save('LTV_Chat_2','LTV_Chat_2');


    case 3 % Weights

        %Compute welfare implied by the solution obtained by fmincon
        load('weights_param','weights_param');
        xopt_fmc = weights_param;

        set_param_value('rho_w_bu',xopt_fmc(1));
        set_param_value('chi_w_bu',xopt_fmc(2));
        set_param_value('rho_w_b',xopt_fmc(3));
        set_param_value('chi_w_b',xopt_fmc(4));

        var_list_={'W_EMU'};
        %info=stoch_simul(var_list_);
        %fhat_fmc=oo_.mean;

        [info, oo_, options_] = stoch_simul(M_, options_, oo_, var_list_); %get decision rules and moments
        fhat_fmc=oo_.mean(strmatch('W_EMU',var_list_,'exact'));

        %Find csminwel solution for the welfare criterion specified in the respective case
        %in welfare_objective_bound
        xopt=xopt_fmc;
        npar=length(xopt);
        [fhat,xhat,ghat,Hhat,itct,fcount,retcodehat] = csminwel('welfare_objective_bound',xopt,eye(npar)*.5,[] ,1e-8,100);

        weights_3 = xhat;
        save('weights_3','weights_3');

        %Check whether csminwel solves for a higher point than fmincon
        check_fmc = -fhat - fhat_fmc


        %Compute the SS consumption level needed to reach a steady state welfare equal to the mean welfare under the respective policy regime
        set_param_value('rho_w_bu',weights_3(1));
        set_param_value('chi_w_bu',weights_3(2));
        set_param_value('rho_w_b',weights_3(3));
        set_param_value('chi_w_b',weights_3(4));


        %var_list_=char('W_S', 'W_B');
        var_list_={'W_S', 'W_B'};


        %info=stoch_simul(var_list_);
        [info, oo_, options_] = stoch_simul(M_, options_, oo_, var_list_); %get decision rules and moments

        Chat_S = exp(((1-beta_p))*oo_.mean(1));

        Chat_B = exp(((1-beta_e))*oo_.mean(2));

        weights_Chat_3 = [Chat_S, Chat_B];

        save('weights_Chat_3','weights_Chat_3');

    case 4 % LAW: augmented with borrowings

        %Compute welfare implied by the solution obtained by fmincon
        load('law_borrowing','law_borrowing');
        xopt_fmc = law_borrowing;

        set_param_value('rho_ib',xopt_fmc(1));
        set_param_value('phi_pie',xopt_fmc(2));
        set_param_value('phi_B',xopt_fmc(3));
        %set_param_value('phi_y',xopt_fmc(4));

        var_list_={'W_EMU'};
        %info=stoch_simul(var_list_);
        %fhat_fmc=oo_.mean;

        [info, oo_, options_] = stoch_simul(M_, options_, oo_, var_list_); %get decision rules and moments
        fhat_fmc=oo_.mean(strmatch('W_EMU',var_list_,'exact'));

        %Find csminwel solution for the welfare criterion specified in the respective case
        %in welfare_objective_bound
        xopt=xopt_fmc;
        npar=length(xopt);
        [fhat,xhat,ghat,Hhat,itct,fcount,retcodehat] = csminwel('welfare_objective_bound',xopt,eye(npar)*.5,[] ,1e-8,100);

        law_borr_4 = xhat;
        save('law_borr_4','law_borr_4');

        %Check whether csminwel solves for a higher point than fmincon
        check_fmc = -fhat - fhat_fmc


        %Compute the SS consumption level needed to reach a steady state welfare equal to the mean welfare under the respective policy regime
        set_param_value('rho_ib',law_borr_4(1));
        set_param_value('phi_pie',law_borr_4(2));
        set_param_value('phi_B',law_borr_4(3));
        %set_param_value('phi_y',law_borr_4(4));

        %var_list_=char('W_S', 'W_B');
        var_list_={'W_S', 'W_B'};


        %info=stoch_simul(var_list_);
        [info, oo_, options_] = stoch_simul(M_, options_, oo_, var_list_); %get decision rules and moments

        Chat_S = exp(((1-beta_p))*oo_.mean(1));

        Chat_B = exp(((1-beta_e))*oo_.mean(2));

        law_borr_Chat_4 = [Chat_S, Chat_B];

        save('law_borr_Chat_4','law_borr_Chat_4');


    case 5 % LAW: augmented capital prices

        %Compute welfare implied by the solution obtained by fmincon
        load('law_asset_prices','law_asset_prices');
        xopt_fmc = law_asset_prices;

        set_param_value('rho_ib',xopt_fmc(1));
        set_param_value('phi_pie',xopt_fmc(2));
        set_param_value('phi_AP',xopt_fmc(3));

        var_list_={'W_EMU'};
        %info=stoch_simul(var_list_);
        %fhat_fmc=oo_.mean;

        [info, oo_, options_] = stoch_simul(M_, options_, oo_, var_list_); %get decision rules and moments
        fhat_fmc=oo_.mean(strmatch('W_EMU',var_list_,'exact'));

        %Find csminwel solution for the welfare criterion specified in the respective case
        %in welfare_objective_bound
        xopt=xopt_fmc;
        npar=length(xopt);
        [fhat,xhat,ghat,Hhat,itct,fcount,retcodehat] = csminwel('welfare_objective_bound',xopt,eye(npar)*.5,[] ,1e-8,100);

        law_asset_prices = xhat;
        save('law_asset_prices','law_asset_prices');

        %Check whether csminwel solves for a higher point than fmincon
        check_fmc = -fhat - fhat_fmc


        %Compute the SS consumption level needed to reach a steady state welfare equal to the mean welfare under the respective policy regime
        set_param_value('rho_ib',law_asset_prices(1));
        set_param_value('phi_pie',law_asset_prices(2));
        set_param_value('phi_y',law_asset_prices(3));


        %var_list_=char('W_S', 'W_B');
        var_list_={'W_S', 'W_B'};


        %info=stoch_simul(var_list_);
        [info, oo_, options_] = stoch_simul(M_, options_, oo_, var_list_); %get decision rules and moments

        Chat_S = exp(((1-beta_p))*oo_.mean(1));

        Chat_B = exp(((1-beta_e))*oo_.mean(2));

        law_asset_prices_Chat_5 = [Chat_S, Chat_B];

        save('law_asset_prices_Chat_5','law_asset_prices_Chat_5');

end