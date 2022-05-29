% Final model for the dissertation

%Estimation works without data on EONIA.

%Comment 1: add inflation channel later on
%Comment 2: add detailed banking

%Setup Dynare launch options
% Added shocks
%Substitute m_e_ss by m_e
% Express mark-ups via elasticities

%Bayesian estimation
@#define bayesian = 0


% =1 for disables OccBin
@#define occbin = 1

% =1 for active MAP
@#define policy_regime = 0

%---------------------Endogenous variables---------------------------------------------%

var 
%----------------------------For Bayesian estimation: data variables-----------------------------------------------%

%Variables
data_C   %1
data_I    %2
data_D       % 3
data_B_E  % 4
data_B_EU   % 5

%Pie
data_PIE      % 6

%Rates
data_R_B % 7 dataseries
data_R_D   % 8 dataseries
data_R_BU  % 9 dataseries
data_EONIA  % 10 dataseries


%----------------------------For Bayesian estimation-----------------------------------------------%
m_e
mk_d      ${\kappa_{kb}}$ (long_name='Elasticity of substitution of deposits')
mk_b_ee   ${\kappa_{kb}}$ (long_name='Elasticity of substitution of secured loans')
mk_b_eeu   ${\kappa_{kb}}$ (long_name='Elasticity of substitution of unsecured loans')
ee_qk
eps_l
eps_K_b
%PIW ${PIW}$ (long_name='Wage inflation')

%----------------------------Modelling defaults-----------------------------------------------%

omega ${\Omega}$ (long_name='Macro-variable of economy indebtness')

%----------------------------Welfare-----------------------------------------------%

%W_EMU ${W^{EMU}}$ (long_name='Welfare deviation to steady state in the economy')

%------------------Utility and Welfare-------------------------------------------%

%V_S ${V^S}$ (long_name='Recursive utility of savers in the core')
%V_B ${V^B}$ (long_name='Recursive utility of borrowers in the core')
%W_S ${W^S}$ (long_name='Welfare deviation to steady state of savers in the core')
%W_B ${W^B}$ (long_name='Welfare deviation to steady state of borrowers in the core')

%---------------------Households---------------------------------------------%

c_p   ${c_p}$ (long_name='Patient HHs consumption')   
d_p   ${d_p}$ (long_name='Patient HHs deposits')   
lam_p ${\lambda_p}$ (long_name='Patient HHs Lagrange multiplier')  
l_p   ${l_p}$ (long_name='Patient HHs supply of labor') 
 

 %---------------------Entrepreneurs---------------------------------------------%

c_e   ${c_e}$ (long_name='Entrepreneurs consumption')   
k_e   ${k_e}$ (long_name='Entrepreneurs capital')    
b_ee  ${b_{ee}}$ (long_name='Entrepreneurs secured loans')
b_eeu ${b_{eeu}}$ (long_name='Entrepreneurs unsecured loans')
lam_e ${\lambda_e}$ (long_name='Entrepreneurs Lagrange multiplier') 
s_e   ${s_e}$ (long_name='Entrepreneurs Lagrange multiplier of borrowing constraint')  
l_pd  ${l_{pd}}$ (long_name='Entrepreneurs demand for labor')  
y_e   ${y_e}$ (long_name='Entrepreneurs output')  
r_k   ${r_k}$ (long_name='Entrepreneurs return on capital')
deltae ${\delta^E}$ (long_name='Proportion of non-repaid unsecured debt')

%---------------------Production, CGPs, Retailers---------------------------------------------%

pie   ${\pi}$ (long_name='Inflation rate')   
mc_E  ${mc_E}$ (long_name='Relative competitive price of wholesale good') 
J_R   ${J_R}$ (long_name='Retailers profits')  
q_k   ${q_k}$ (long_name='Price of capital')  
x     ${x}$ (long_name='Inverse of mc_E')  
I     ${I}$ (long_name='Capital investment')

%---------------------Aggregation & Equilibrium---------------------------------------------%

C     ${C}$ (long_name='Aggregate consumption')  
Y     ${Y}$ (long_name='Aggregate output')    
w_p   ${w_p}$ (long_name='Wage rate')    
B     ${B}$ (long_name='Aggregate debt')    
D     ${D}$ (long_name='Aggregate deposits')    
K     ${K}$ (long_name='Aggregate capital')  
 
%---------------------Monetary Policy---------------------------------------------%

r_ib  ${r_{ib}}$ (long_name='Policy rate')   
 
%---------------------Banks---------------------------------------------%
K_b   ${K_b}$ (long_name='Bank capital')
J_B   ${J_B}$ (long_name='Bank profit')   
r_d   ${r_d}$ (long_name='Retail rate on deposits')
r_b   ${r_b}$ (long_name='Retail rate on secured loans')  
r_bu  ${r_{bu}}$ (long_name='Retail rate on unsecured loans')   
R_b        ${R_b}$     (long_name='Rate on deposits')  
R_bu       ${R_bu}$    (long_name='Wholesale rate on unsecured loans')
lev        ${lev}$     (long_name='Bank leverage')  
%rr         ${rr}$      (long_name='Real rate on loans')  
Y1         ${Y1}$      (long_name='Definition of output used in policy rules')  %//Auxiliary variable 1
B_rw       ${B_{rw}}$  (long_name='Risk-weighted loans') 
weight_bu  ${w_{bu}}$  (long_name='Risk weight on unsecured loans') 
B_eeu      ${B_{eeu}}$ (long_name='Unsecured loans in bank')  
weight_b   ${w_{b}}$   (long_name='Risk weight on secured loans') 
B_ee       ${B_{ee}}$  (long_name='Secured loans in bank') 
 


 %---------------------Exogenous processes---------------------------------------------%

mk_y  ${mk_y}$ (long_name='Retailers mark-up')
A_e   ${A_e}$ (long_name='Total factor productivity') 

%---------------------Macropudential Policy---------------------------------------------%

%Passive
@#if policy_regime == 0
vi ${\vi}$ (long_name='Macroprudential policy instrument')
@#endif

%Active
@#if policy_regime == 1
vi ${\vi}$ (long_name='Macroprudential policy instrument')
@#endif

;     

%---------------------Exogenous variables---------------------------------------------%
varexo  e_A_e e_mk_y 
%---------------------For Bayesian---------------------------------------------%
        e_m_e e_mk_d e_mk_b_ee e_mk_b_eeu e_qk e_l e_eps_K_b e_r_ib %last is Taylor rule shock
;

%---------------------Parameters---------------------------------------------%

parameters  
            %----------------------------Modelling defaults-----------------------------------------------%
            gamma     ${\gamma}$ (long_name='Default amplification in Omega')
            omega_r   ${\omega_r}$ (long_name='Credit to GDP amplification in Omega')
            kappa_u   ${\kappa_u}$ (long_name='Cost of renegotiation of unsecured debt')
            omega_ss  ${\omega_{ss}}$ (long_name='Steady state level of \Omega') %if use steady_state(omega), then there is a collinearity problem
            r_bu_ss   ${\r^{bu}_{ss}}$ (long_name='Steady state level on unsecured retail rate')
            %---------------------Households and entrepreneurs---------------------------------------------%
            beta_p ${\beta_p}$ (long_name='Discount factor of HHs')
            beta_e ${\beta_e}$ (long_name='Discount factor of entrepreneurs')
            phi    ${\phi}$ (long_name='Inverse elasticity of labor supply')
           % m_e_ss ${m^e_{ss}}$ (long_name='LTV cap')  
            gamma_p ${\gamma_p}$ (long_name='Relative HH weight')
            gamma_e ${\gamma_e}$ (long_name='Relative entrepreneurs weight')
            %---------------------Banks---------------------------------------------%
            vi_ss ${\nu_{ss}}$ (long_name='Capital requirement in SS')  
            %mcspread ${mcspread}$ (long_name='Loan-deposit spread')  
            delta_b  ${\delta_b}$ (long_name='Bank capital depreciation rate')
            %---------------------Production and retailers---------------------------------------------%
            eps_y    ${\epsilon_y}$ (long_name='For mk_y_ss calculations')
            mk_y_ss  ${\delta_k}$ (long_name='Retailers mark-up over wholesale price in steady state')
            ksi      ${\delta_k}$ (long_name='Fraction of capital used in production')
            kappa_p  ${\delta_k}$ (long_name='Cost of adjustment of retailers')
            kappa_i  ${\delta_k}$ (long_name='Cost of adjustment of capital goods producers')
            deltak   ${\delta_k}$ (long_name='Capital depreciation rate')
            piss     ${\pi_{ss}}$ (long_name='Steady state inflation')
            %---------------------Monetary policy---------------------------------------------%
            r_ib_ss ${r^{ib}_{ss}}$ (long_name='Steady state policy rate')
            rho_ib  ${\rho^{ib}}$ (long_name='MP policy rate inertia')
            phi_pie ${\phi^{\pi}}$ (long_name='MP inflation inertia')
            phi_y   ${\phi^{y}}$ (long_name='MP output gap inertia')
            phi_AP  ${\phi^{AP}}$ (long_name='MP asset prices gap')
            phi_B   ${\phi^{B}}$ (long_name='MP borrowing gap inertia')
            %---------------------Shocks---------------------------------------------%        
            rho_A_e  ${\rho^A_e}$ (long_name='AR(1) coefficient of TFP process')
            rho_mk_y ${\rho^{mk}_y}$ (long_name='AR(1) coefficient of mark-up process')            
            rho_w_bu ${\rho^{bu}_w}$ (long_name='Unsecured risk weight inertia')            
            chi_w_bu ${\chu^{bu}_w}$ (long_name='Unsecured risk weight output gap inertia')            
            rho_w_b  ${\rho^{b}_w}$ (long_name='Secured risk weight inertia')
            chi_w_b  ${\chu^{b}_w}$ (long_name='Secured risk weight output gap inertia') 


            %---------------------For bayesian estimation---------------------------------------------%     
            m_e_ss
            rho_m_e
            mk_d_ss    
            mk_b_ee_ss   
            mk_b_eeu_ss
            rho_mk_d
            rho_mk_b_eeu
            rho_mk_b_ee
            rho_ee_qk
            rho_eps_l
            rho_eps_K_b
            
            r_b_ss
            r_d_ss

            %---------------------Macroprudential policy---------------------------------------------%
            @#if policy_regime == 0 
            rho_vi ${\rho_vi}$ (long_name='MAP inertia')
            phi_vi ${\phi_vi}$ (long_name='MAP output gap inertia')
            @#endif

            @#if policy_regime == 1 
            rho_vi ${\rho_vi}$ (long_name='MAP inertia')
            phi_vi ${\phi_vi}$ (long_name='MAP output gap inertia')
            @#endif

            ind_d     ${ind_d}$ (long_name='Indexation, deposit rates')    
            ind_be    ${ind_be}$ (long_name='Indexation, rates on loans to firms')
            ind_bh    ${ind_bh}$ (long_name='Indexation, rates on unsecured loans to firms')

            kappa_d     ${\kappa_d}$       (long_name='Adjustment cost parameter for deposit rates')
            kappa_b_ee  ${\kappa_{b_ee}}$  (long_name='Adjustment cost parameter for rates on secured loans')
            kappa_b_eeu ${\kappa_{b_eeu}}$ (long_name='Ddjustment cost parameter for rates on unsecured loans')

            kappa_kb    ${\kappa_{kb}}$ (long_name='Adjustment cost parameter for Banking Capital (Basel II)')

             eps_d       ${\kappa_{kb}}$ (long_name='Elasticity of substitution of deposits')
             eps_b_ee    ${\kappa_{kb}}$ (long_name='Elasticity of substitution of secured loans')
             eps_b_eeu   ${\kappa_{kb}}$ (long_name='Elasticity of substitution of unsecured loans')
           ;

%------------------------------------Assignment of parameter values------------------------------------------------------------------------------------------%           
%---------------------For bayesian estimation---------------------------------------------%     
eps_d          = -1.46025;                          
eps_b_ee       = 3.154542134;
eps_b_eeu      = 2.117626240126008;

rho_m_e = 0.976;
mk_d_ss = eps_d   / (eps_d  - 1) ;
mk_b_ee_ss = eps_b_ee  / (eps_b_ee - 1) ;   
mk_b_eeu_ss = eps_b_eeu  / (eps_b_eeu - 1) ;
rho_mk_d= 0.6036;
rho_mk_b_eeu= 0.7623;
rho_mk_b_ee= 0.8042;

rho_ee_qk = 0.5259;
rho_eps_l = 0.1199;
rho_eps_K_b = 0.2070;

%Macro-variable rule and default
omega_r  = 0.4958;
gamma    = 1.5041;
omega_ss = 311.0046484546306;
kappa_u  = 38.8674;

%Parameters related to the main model block           
beta_p   = 0.996;
beta_e   = 0.975;
%beta_e   = 0.93571986;
phi      = 1;           % Robustness: High: 1.5, Low: 0.5
m_e_ss   = 0.35;        % Hi-LTV case: m_e_ss = 0.70
gamma_p  = 1;
gamma_e  = 1 ;
eps_y    = 6; 
mk_y_ss  = eps_y   / (eps_y  - 1);
ksi    = 0.20;
kappa_p  =  288.4960;      % Robustness: High: 100, Low: 15
kappa_i  = 4.6531;           % For Tech. Shock. For cost-push shock, subst with kappa_i = 0.05
                        % Robustness: Tech.Sh.: kappa_i = 0.05, 0.5, 10 
                        % Cost push shock: kappa_i = 5, 0.5 0.005
deltak    = 0.050;
piss     = 1.00;
vi_ss           = 0.09; % Robustness: High: 0.15, Low: 0.045
%mcspread = 0.0050;      % = 2% annual
delta_b      =  0.0721; % Calibrated such that the multiplier on b_ee is positive, hence b_ee_u is positive/non-negative
rho_A_e      = 0.9806;    % Robustness: rho_A_e  = 0, rho_A_e  = 0.50
rho_mk_y     = 0.50;    % Robustness: rho_mk_y = 0, rho_mk_y = 0.75
r_ib_ss  = (piss/beta_p - 1) * (1-eps_d)/(-eps_d) ;
r_bu_ss = (1/beta_e) - 1;
r_d_ss = (1/beta_p - 1);
r_b_ss = r_ib_ss * eps_b_ee/(eps_b_ee-1); 
ind_d        = 0.0;              
ind_be       = 0.0;              
ind_bh       = 0.0;              
kappa_d      = 14.5261;              
kappa_b_ee     = 28.5744;             
kappa_b_eeu     = 112.9446;            
kappa_kb     = 1.5947;               


%Parameters related to MAP policy
@#if policy_regime == 0
rho_vi      =   0    ; 
phi_vi      =   0    ; 
@#endif

@#if policy_regime == 1
rho_vi      =   0.8     ; 
phi_vi      =   2.35    ; 
@#endif

%Parameters related to MP policy
rho_ib   = 0.7555;        % Robustness: rho_ib = 0
phi_pie  = 2.8725;         % MP RULE parameter, to be changed for simulations
phi_y    = 0.2121;           % MP RULE parameter, to be changed for simulations
phi_AP   = 0;           % MP RULE parameter, to be changed for simulations
phi_B    = 0;           % MP RULE parameter, to be changed for simulations

%Parameters related to risk weights
rho_w_bu = 0.93;
chi_w_bu = -12.5;
rho_w_b  = 0.93;
chi_w_b  = -12.5;

%--------------------------------Model Block-----------------------------------------------------------------------------------------------------------------%
model;

%------------------------------------Adjustment costs and weights----------------------------------------------------%

%Adjustments for wholesale unit
%#adj_bee = - kappa_kb * ( exp(K_b)/exp(B_rw) - exp(vi) ) * ((exp(K_b)/exp(B_rw))^2) * exp(weight_b);
%#adj_beeu = - kappa_kb * ( exp(K_b)/exp(B_rw) - exp(vi) ) * ((exp(K_b)/exp(B_rw))^2) * exp(weight_bu); 
%Disutility from default for entrepreneurs
%#quadratic_cost = 0.5*exp(omega)*((1+exp(r_bu))*exp(b_eeu)*exp(deltae))^2;
%Adjustments for retail units
%#adj_retail_d_p = - kappa_d  * ( exp(r_d)/exp(r_d(-1)) - 1  )  * exp(r_d)/exp(r_d(-1))
%+ beta_p * ( exp(lam_p(+1))/exp(lam_p) ) * kappa_d  * ( exp(r_d(+1))/exp(r_d) - ( exp(r_d)/exp(r_d(-1)))^ind_d )   * ( (exp(r_d(+1))/exp(r_d))^2 )   * (exp(D(+1))/exp(D));
%#adj_retail_b_ee = - kappa_b_ee * (exp(r_b)/exp(r_b(-1)) - 1 ) * exp(r_b)/exp(r_b(-1))
 %  + beta_p * ( exp(lam_p(+1))/exp(lam_p) ) * kappa_b_ee * ( exp(r_b(+1))/exp(r_b) - ( exp(r_b)/exp(r_b(-1)))^ind_be ) * ( (exp(r_b(+1))/exp(r_b))^2 ) * (exp(B_ee(+1))/exp(B_ee));
%#adj_retail_b_eeu = - kappa_b_eeu * (exp(r_bu)/exp(r_bu(-1)) - 1 ) * exp(r_bu)/exp(r_bu(-1)) 
%+ beta_p * ( exp(lam_p(+1))/exp(lam_p) ) * kappa_b_eeu * ( exp(r_bu(+1))/exp(r_bu) - ( exp(r_bu)/exp(r_bu(-1)))^ind_bh ) * ( (exp(r_bu(+1))/exp(r_bu))^2 ) * (exp(B_eeu(+1))/exp(B_eeu));

[name = 'Risk weight for secured loans']
%exp(weight_b) = (1 - rho_w_b) * 1 + (1 - rho_w_b) * chi_w_b * (exp(Y) - exp(Y(-4))) + rho_w_b *  exp(weight_b(-1));
exp(weight_b) = (1);
[name = 'Risk weight for unsecured loans']
%exp(weight_bu) = (1 - rho_w_bu) * 1 + (1 - rho_w_bu) * chi_w_bu * (exp(Y) - exp(Y(-4))) + rho_w_bu *  exp(weight_bu(-1));
exp(weight_bu) = (1);
%------------------------------------Modelling defaults----------------------------------------------------%

[name = 'Macro-variable']
%A-la Peiris specification
%omega = steady_state(omega) * ( ( steady_state(b_eeu) * Y * (1 + steady_state(r_bu)) ) / ( steady_state(Y) * b_eeu * (1+r_bu) ) )^omega_r 
% * ( steady_state(deltae) / (deltae) )^gamma;

 %Simplified specification
%omega = steady_state(omega) * (( steady_state(b_eeu)  / b_eeu )^(omega_r)) * (( steady_state(deltae) / deltae )^gamma);
%omega = steady_state(omega) * (( (b_eeu)  / b_eeu(-1) )^(omega_r)) * (( (deltae) / deltae(-1) )^gamma);
exp(omega) = omega_ss * (( exp(b_eeu)  / exp(b_eeu(-1)) )^(omega_r)) * (( exp(deltae) / exp(deltae) )^gamma);

%For macro-variable only case
%deltae = deltae_ss;

% Defaults related block for ENTREPRENEURS
[name = 'Ent: FOC with respect to delta_e']
-exp(omega) * exp(b_eeu(-1)) * exp(deltae) * (1 + exp(r_bu(-1))) + exp(lam_e)/exp(pie) = 0; 

[name = 'Ent: FOC with respect to unsecured debt']
exp(lam_e) * (1 - kappa_u * ( exp(b_eeu) - exp(steady_state(b_eeu)) ) ) 
+ beta_e * ( -exp(omega(+1)) * (exp(deltae(+1)))^2 * exp(b_eeu) * (1 + exp(r_bu))^2 
- exp(lam_e) * (1 - exp(deltae(+1))) * (1 + exp(r_bu)) )/exp(pie)  = 0; 

[name = 'Ent: Budget constraint']
exp(c_e) + (1+exp(r_b(-1))) * exp(b_ee(-1))/exp(pie)  + (1-exp(deltae)) * (1 + exp(r_bu(-1))) * exp(b_eeu(-1)) + exp(w_p)*exp(l_pd)  + exp(q_k) * exp(k_e) 
   = exp(y_e) / exp(x) + exp(b_ee) + exp(b_eeu) + exp(q_k) * (1-deltak) * exp(k_e(-1)) + 0.5*kappa_u*(exp(b_eeu) - exp(steady_state(b_eeu)))^2 ; 

%------------------------------------Households----------------------------------------------------%c, l, lam, d

[name = 'HH: Marginal utility of consumption']
(exp(c_p))^(-1) = exp(lam_p);

[name = 'HH: Euler equation']
exp(lam_p) = beta_p*exp(lam_p(+1))*(1+exp(r_d));

[name = 'HH: Labor supply decision']
exp(l_p)^phi = ( exp(eps_l) / (1 - exp(eps_l)) ) * exp(lam_p)*exp(w_p);

[name = 'HH: Budget constraint'] %Added profits from banks
exp(c_p) + exp(d_p) = exp(w_p)*exp(l_p)+(1+exp(r_d(-1)))*exp(d_p(-1))/exp(pie) +exp(J_R)/gamma_p;% + (J_B)/gamma_p;
 

%------------------------------------Entrepreneurs----------------------------------------------------% lam_e, c, s_e, b_u, b_ee, delta, k

[name = 'Ent: Marginal utility of consumption']
(exp(c_e))^(-1) = exp(lam_e); 

[name = 'Ent: Labor demand']
exp(w_p) = (1-ksi)*exp(y_e)/(exp(l_pd)*exp(x)); 

[name = 'Ent: Euler equation']
exp(s_e) * exp(m_e) * exp(q_k(+1))*(1-deltak)*exp(pie(+1))/(1+exp(r_b)) + beta_e*exp(lam_e(+1))*(exp(q_k(+1))*(1-deltak) + exp(r_k(+1))) 
    = exp(lam_e) * exp(q_k) ;

[name = 'Ent: Consumption Euler equation']    
exp(lam_e)-exp(s_e) = beta_e*exp(lam_e(+1))*(1+exp(r_b))/exp(pie(+1));

%[name = 'Ent: Budget constraint']
%(c_e) + (1+r_b(-1)) * (b_ee(-1))  +  (w_p)*(l_pd)  + (q_k) * (k_e) 
%   = (y_e) / (x) + (b_ee) + (q_k) * (1-deltak) * (k_e(-1))  ;    

[name = 'Ent: Production function']   
exp(y_e) = exp(A_e) * ( exp(k_e(-1)) )^ksi * ( exp(l_pd) ) ^ (1-ksi);

[name = 'Ent: Borrowing constraint']
(1+exp(r_b)) * exp(b_ee) = exp(m_e) * exp(q_k(+1))  * exp(k_e) * (1-deltak)*exp(pie(+1));

[name = 'Ent: Return to capital']
exp(r_k) = ksi * exp(A_e) * (exp(k_e(-1)))^(ksi-1) * ( exp(l_pd) ) ^ (1-ksi) /exp(x); 

%------------------------------------Capital Producers----------------------------------------------------% 

[name = 'Capital accumulation']
exp(K) = (1-deltak) * exp(K(-1)) + ( 1 - kappa_i/2 * (exp(I)*exp(ee_qk)/exp(I(-1)) - 1)^2 ) * exp(I) ;

[name = 'Capital producers FOC']
1 = exp(q_k) * ( 1 -  kappa_i/2 * (exp(I)*exp(ee_qk)/exp(I(-1)) - 1)^2  - kappa_i * (exp(I)*exp(ee_qk)/exp(I(-1)) - 1) * exp(I)*exp(ee_qk)/exp(I(-1)) ) 
  + beta_e * exp(lam_e(+1)) / exp(lam_e) * exp(q_k(+1)) *   kappa_i * (exp(I(+1)*exp(ee_qk(+1)))/exp(I) - 1)  * (exp(I(+1)*exp(ee_qk))/exp(I))^2 ;

    
%------------------------------------Retailers----------------------------------------------------% 

[name = 'NK Phillips curve']
1 - exp(mk_y)/(exp(mk_y)-1)+ exp(mk_y)/(exp(mk_y)-1)*exp(mc_E) - kappa_p*(exp(pie) - 1)*exp(pie) 
   + beta_p*exp(lam_p(+1))/exp(lam_p)*kappa_p*(exp(pie(+1))-1)*exp(pie(+1))*(exp(Y(+1))/exp(Y)) = 0;

[name = 'mc_E definition']
exp(mc_E) = 1 / exp(x);

[name = 'Retailer profits']
exp(J_R) = exp(y_e) * ((1 - exp(mc_E)) - 0.5*kappa_p*(exp(pie) - 1)^2);

%------------------------------------Monetary Policy----------------------------------------------------% 

[name = 'Interest rate rule'] %Substitute r_ib_ss for steady_state(r_ib_ss)
(1+exp(r_ib)) = (1 + e_r_ib) * (1+r_ib_ss)^(1-rho_ib)*(1+exp(r_ib(-1)))^rho_ib*((exp(pie)/piss)^phi_pie*(exp(Y1)/exp(steady_state(Y1)))^phi_y  
           *(exp(q_k)/exp(steady_state(q_k)))^phi_AP*(exp(B)/exp(steady_state(B)))^phi_B)^(1-rho_ib);



%------------------------------------Banks----------------------------------------------------% 

[name = 'Aggregate bank profits'] %add adjustment costs
exp(J_B)  = (exp(r_b)*exp(b_ee) + exp(r_bu)*exp(b_eeu)*(1-exp(deltae)) - exp(r_d)*exp(d_p)
- kappa_kb/2*(exp(K_b)/exp(B_rw) - exp(vi))^2 *exp(K_b))
           - kappa_d/2  * ( (exp(r_d)/exp(r_d(-1))-1)^2)   * exp(r_d) *exp(d_p) 
           - kappa_b_ee/2 * ( (exp(r_b)/exp(r_b(-1))-1)^2) * exp(r_b)*exp(b_ee) 
           - kappa_b_eeu/2 * ( (exp(r_bu)/exp(r_bu(-1))-1)^2) * exp(r_bu)*exp(b_eeu);


[name = 'Risk-weighted assets']
exp(B_rw) = exp(weight_bu) * exp(B_eeu) + exp(weight_b) * exp(B_ee);

[name = 'Wholesale rate on secured loans']
exp(R_b) = exp(r_ib) - kappa_kb * ( exp(K_b)/exp(B_rw) - exp(vi) ) * ((exp(K_b)/exp(B_rw))^2) * exp(weight_b) ; 

[name = 'Wholesale rate on unsecured loans']
exp(R_bu) = exp(r_ib) - kappa_kb * ( exp(K_b)/exp(B_rw) - exp(vi) ) * ((exp(K_b)/exp(B_rw))^2) * exp(weight_bu) ; 

[name = 'Retail rate on unsecured loans']
(1 - exp(deltae)) * (1 - exp(mk_b_eeu)/(exp(mk_b_eeu)-1)) + exp(mk_b_eeu)/(exp(mk_b_eeu)-1)  * exp(R_bu)/exp(r_bu) 
- kappa_b_eeu * (exp(r_bu)/exp(r_bu(-1)) - 1 ) * exp(r_bu)/exp(r_bu(-1)) 
+ beta_p * ( exp(lam_p(+1))/exp(lam_p) ) * kappa_b_eeu * ( exp(r_bu(+1))/exp(r_bu) - ( exp(r_bu)/exp(r_bu(-1)))^ind_bh ) * ( (exp(r_bu(+1))/exp(r_bu))^2 ) * (exp(B_eeu(+1))/exp(B_eeu)) = 0;

[name = 'Retail rate on secured loans']
+ 1 - exp(mk_b_ee)/(exp(mk_b_ee)-1)  +  exp(mk_b_ee)/(exp(mk_b_ee)-1) * exp(R_b)/exp(r_b) 
- kappa_b_ee * (exp(r_b)/exp(r_b(-1)) - 1 ) * exp(r_b)/exp(r_b(-1))
   + beta_p * ( exp(lam_p(+1))/exp(lam_p) ) * kappa_b_ee * ( exp(r_b(+1))/exp(r_b) - ( exp(r_b)/exp(r_b(-1)))^ind_be ) * ( (exp(r_b(+1))/exp(r_b))^2 ) * (exp(B_ee(+1))/exp(B_ee)) = 0;

[name = 'Retail rate on deposits']
- 1 + exp(mk_d)/(exp(mk_d)-1)  - exp(mk_d)/(exp(mk_d)-1) * exp(r_ib)/exp(r_d) 
- kappa_d  * ( exp(r_d)/exp(r_d(-1)) - 1  )  * exp(r_d)/exp(r_d(-1))
+ beta_p * ( exp(lam_p(+1))/exp(lam_p) ) * kappa_d  * ( exp(r_d(+1))/exp(r_d) - ( exp(r_d)/exp(r_d(-1)))^ind_d )   * ( (exp(r_d(+1))/exp(r_d))^2 )   * (exp(D(+1))/exp(D)) = 0; 

[name = 'Aggregate bank capital']
exp(K_b)*exp(pie)  = ((1-delta_b)*exp(K_b(-1))/exp(eps_K_b) + exp(J_B(-1)));

%[name = 'Real rate']
%rr = r_b - (pie);

[name = 'Output definition for policies']
exp(Y1)   = exp(C)+exp(I);

[name = 'Leverage']
exp(lev)  =  exp(B)/exp(K_b);

%------------------------------------Aggregation and Equilibrium----------------------------------------------------% 

[name = 'Secured loan market clearing'] 
exp(B_ee) =  exp(b_ee);
[name = 'Unsecured loan market clearing'] 
exp(B_eeu) = exp(b_eeu);
[name = 'Aggregate loan market clearing'] 
exp(B)              = gamma_e*exp(b_ee) + gamma_e*exp(b_eeu);
[name = 'Aggregate consumption']
exp(C)              = gamma_p*exp(c_p) + gamma_e*exp(c_e); 
[name = 'Labor market clearing']
gamma_e * exp(l_pd) = gamma_p*exp(l_p); 
[name = 'Deposit market clearing']
exp(D)              = gamma_p*exp(d_p) ; 
[name = 'Capital market clearing']
exp(K)              = gamma_e*exp(k_e); 
[name = 'Bank balance sheet']
exp(B)              = exp(D) + exp(K_b);  
[name = 'Output clearing']
exp(Y)              = gamma_e*exp(y_e); 
%[name = 'Aggregate outpute definition']
%(Y)              = (C) + (q_k)*((K) - (1-deltak)*(K(-1))) + delta_b*(K_b(-1)) - deltae*b_eeu; 


%---------------------------------------Equations for MAP--------------------------------------------------------------------------------%

%Passive
@#if policy_regime == 0
[name = 'Capital requirements rule']
exp(vi) = vi_ss;
@#endif

%Active
@#if policy_regime == 1

[name = 'Capital requirements rule']
vi = vi_ss^(1 - rho_vi) * vi(-1)^rho_vi * 
             ((((B)/(Y1))/((steady_state(B))/(steady_state(Y1))))^phi_vi)^(1-rho_vi);

@#endif


%---------------------------------------Equations for welfare--------------------------------------------------------------------------------%

% [name = 'Recursive utility of households (savers)']
% V_S = (log((c_p)) - ((l_p)^(1+phi)/(1+phi))) + beta_p*V_S(+1);
% 
% [name = 'Recursive utility of entrepreneurs (borrowers)']
% V_B = (log(c_e) - quadratic_cost) + beta_e*V_B(+1);
% 
% [name = 'Welfare criterion of households (savers)']
% W_S = V_S -steady_state(V_S);
% 
% [name = 'Welfare criterion of entrepreneurs (borrowers)']
% W_B = V_B -steady_state(V_B);
% 
% [name = 'Aggregate welfare criterion in the economy']
% W_EMU = 0.5*W_S + 0.5*W_B;

%------------------------------------Exogenous processes (shocks, should be 9)----------------------------------------------------% 

[name = 'Technology shock']
exp(A_e)  = 1 - rho_A_e + rho_A_e*exp(A_e(-1)) + e_A_e;

[name = 'Goods mark-up shock']
exp(mk_y) = (1-rho_mk_y)*mk_y_ss + rho_mk_y*exp(mk_y(-1)) + e_mk_y;


%------------------------------------For Bayesian estimation----------------------------------------------------% 
[name = 'Shock to LTV of entrepreneurs shock']
exp(m_e) = (1-rho_m_e)*m_e_ss + rho_m_e*exp(m_e(-1)) + e_m_e;

[name = 'Shock to deposit markdown']
exp(mk_d)     = (1-rho_mk_d)   * mk_d_ss       + rho_mk_d   * exp(mk_d(-1))    + e_mk_d;%

[name = 'Shock to secured loans markup']
exp(mk_b_ee)    = (1-rho_mk_b_ee)  * mk_b_ee_ss      + rho_mk_b_ee  * exp(mk_b_ee(-1))   + e_mk_b_ee;%

[name = 'Shock to unsecured loans markup']
exp(mk_b_eeu)    = (1-rho_mk_b_eeu)  * mk_b_eeu_ss      + rho_mk_b_eeu  * exp(mk_b_eeu(-1))   + e_mk_b_eeu;%

[name = 'Investment specific technology shock']
exp(ee_qk)    =  (1-rho_ee_qk)   *    1          + rho_ee_qk  * exp(ee_qk(-1))   + e_qk;%

[name = 'Labor supply shock']
exp(eps_l)    = (1-rho_eps_l)  * 0.5      + rho_eps_l  * exp(eps_l(-1))   + e_l;%

[name = 'Bank capital shock']
exp(eps_K_b)  = (1-rho_eps_K_b)*    1          + rho_eps_K_b* exp(eps_K_b(-1)) + e_eps_K_b;%

%----------------------------MEASUREMENT EQUATIONS (VARIABLES TAKEN TO THE DATA)------------------------------%

%Variables
data_C       =  C         - steady_state(C); %1
data_I       =  I         - steady_state(I); %2
data_D       =  D         - steady_state(D);   % 3
data_B_E     =  B_ee      - steady_state(B_ee) ;% 4
data_B_EU    =  B_eeu     - steady_state(B_eeu) ;% 5

%Pie
data_PIE     =  exp(pie)  - piss;  % 6

%Rates
data_R_B     =  exp(r_b)  - r_b_ss;  % 7 dataseries
data_R_D     =  exp(r_d)  - r_d_ss;  % 8 dataseries
data_R_BU    =  exp(r_bu) - r_bu_ss;  % 9 dataseries
data_EONIA   =  exp(r_ib) - r_ib_ss; % 10 dataseries
%data_R_B     =  (r_b)  - log(r_b_ss);  % 7 dataseries
%data_R_D     =  (r_d)  - log(r_d_ss);  % 8 dataseries
%data_R_BU    =  (r_bu) - log(r_bu_ss);  % 9 dataseries
%data_EONIA   =  (r_ib) - log(r_ib_ss); % 10 dataseries


end;

initval;
%------------------Bayesian------------------%
%Variables
data_C       =  C; %1
data_I       =  I; %2
data_D       =  0;   % 3
data_B_E     =  0;% 4
data_B_EU    = 0;% 5

%Pie
data_PIE     =  0;  % 6

%Rates
data_R_B     =  0;  % 7 dataseries
data_R_D     =  0;  % 8 dataseries
data_R_BU    = 0;  % 9 dataseries
data_EONIA   =  0; % 10 dataseries

%--------------------------------------------%

m_e             =  		 log(m_e_ss);
mk_d          =  		 log( mk_d_ss);
mk_b_ee       =  		 log(mk_b_ee_ss);
mk_b_eeu       =  		 log( mk_b_eeu_ss);
ee_qk        =  		 0;
eps_l          =  		 log(0.5);
eps_K_b       =  		 0;
%PIW             =  		 0;
omega           =  		 log(omega_ss);
% W_EMU         =  		 log(    		 0
% V_S            =  		 log(   		 -118.83
% V_B           =  		 log(    		 -109.46
% W_S               		 0
% W_B               		 0
c_p               =  		 log(    		 0.926072);
d_p               =  		 log(    		 0.854046);
lam_p             =  		 log(    		 1.07983);
l_p               =  		 log(    		 0.892769);
c_e               =  		 log(    		 0.0822021);
k_e               =  		 log(    		 2.61888);
b_ee             =  		 log(     		 0.862237);
b_eeu             =  		 log(    		 0.0762755);
lam_e             =  		 log(    		 12.1651);
s_e               =  		 log(    		 0.186624);
l_pd                =  		 log(  		 0.892769);
y_e              =  		 log(     		 1.10717);
r_k              =  		 log(     		 0.0704607);
deltae           =  		 log(     		 0.5);
pie               =  		 0;
mc_E             =  		 log( 1 + (1 - mk_y_ss)/mk_y_ss);
J_R              =  		 log(     		 0.184528);
q_k             =  		 log(      		 1);
x                =  		 log(     		 1.2);
I               =  		 log(      		 0.130944);
C               =  		 log(      		 1.00827);
Y               =  		 log(      		 1.10717);
w_p             =  		 log(      		 0.826769);
B               =  		 log(      		 0.938513);
D               =  		 log(      		 0.854046);
K               =  		 log(      		 2.61888);
r_ib                 =  		 log(r_ib_ss);
K_b               =  		 log(    		 0.0844661);
J_B               =  		 log(    		 0.00609001);
r_d               =  		 log(r_d_ss);
r_b               =  		 log(r_b_ss);
r_bu              =  		 log(r_bu_ss);
R_b                =  		 log(r_ib_ss);
R_bu               =  		 log(r_ib_ss);
lev                =  		 log( 11.1111);
Y1                 =  		 log( 1.13922);
B_rw              	   =  		 log( 	 0.938513);
weight_bu         	   =  		 0;
B_eeu             	   =  		 log( 	 0.0762755);
weight_b             =  		0;
B_ee                 =  		 log(0.862237);
mk_y              	   =  		 log(mk_y_ss);
A_e               	   =  		 log(1);
vi                	   =  		 log(vi_ss);
end;

steady(solve_algo=0);
resid;
model_diagnostics;


shocks;
%var e_A_e            =  1;
var e_A_e            =  (0.0072)^2;
var e_mk_y               =  (0.2572)^2;
var e_m_e                =  (0.0173)^2;
var e_mk_d               =  (0.5905)^2;
var e_mk_b_ee            =  (0.0117)^2;
var  e_mk_b_eeu           =  (1.4289)^2;
var e_qk                 =  (0.0053)^2;
var e_l                  =  (0.0552)^2;
var e_eps_K_b            =  (0.7376)^2;
var e_r_ib               =  (0.0034)^2;
end;
%   
% 
options_.noprint=0;
options_.nofunctions=1;
options_.nograph=1;
stoch_simul(order=1, irf=0, periods=10000) C I pie r_ib r_b r_bu d_p r_d b_ee b_eeu;



%-----------------------------------------------------BAYESIAN ESTIMATION---------------------------------------------------%  

@#if bayesian == 1
% Defining which dataseries are used in estimation

%varobs data_C data_I data_rBE data_rD data_rIB data_D data_BE data_PIE;% data_R_BU;
varobs      data_C
            data_I      
            data_D     
            data_B_E   
            data_B_EU  

            %Pie
            data_PIE   

            %Rates
            data_R_B   
            data_R_D   
            data_R_BU  
            %data_EONIA 
            ;

load(['data',filesep,'final_data.mat']); % load dataseries from file
%load(['data',filesep,'end_period_MRO.txt']); % load dataseries from file

%mro = (1+end_period_MRO./100).^0.25 - mean( (1+end_period_MRO./100).^0.25 );

% delete initial/final dates if there are NaN in the data used ...

BEG_OBS          = 2;                          % 1998:Q1
END_OBS          = length(DATES);              % 2009:Q1    53 quarterly observations
DATES            = DATES(BEG_OBS:END_OBS,:);
datadescr.DATES  = DATES;
datadescr.MODEL  = M_.fname;

% give the dataseries the names listed in the 'varobs' command ...

dataseries.data_C    = C_HP (BEG_OBS:END_OBS);
datadescr.data_C     = 'Real Consumption';

dataseries.data_I    = I_HP  (BEG_OBS:END_OBS);
datadescr.data_I     = 'Real Investment';

dataseries.data_D    = D_HP  (BEG_OBS:END_OBS);
datadescr.data_D     = 'Real Deposits';

dataseries.data_B_E    = B_E_HP  (BEG_OBS:END_OBS);
datadescr.data_B_E     = 'Real secured loans';

dataseries.data_B_EU    = B_EU_HP  (BEG_OBS:END_OBS);
datadescr.data_B_EU     = 'Real unsecured loans';

dataseries.data_PIE  = PIE (BEG_OBS:END_OBS) - mean(PIE (BEG_OBS:END_OBS));
datadescr.data_PIE   = 'Inflation (q/q)';

dataseries.data_R_B  = R_B (BEG_OBS:END_OBS) - mean(R_B (BEG_OBS:END_OBS));
datadescr.data_R_B   = 'Interest rate on secured loans';

dataseries.data_R_BU  = R_BU (BEG_OBS:END_OBS) - mean(R_BU (BEG_OBS:END_OBS));
datadescr.data_R_BU   = 'Interest rate on unsecured loans';

dataseries.data_R_D  = R_D (BEG_OBS:END_OBS) - mean(R_D (BEG_OBS:END_OBS));
datadescr.data_R_D   = 'Interest rate on deposits';

dataseries.data_EONIA  = EONIA (BEG_OBS:END_OBS) - mean(EONIA (BEG_OBS:END_OBS));
datadescr.data_EONIA   = 'Short-term interest rate';

%--------------------OLD--------------------%
% 
% dataseries.data_rBE  = RBE (BEG_OBS:END_OBS) - mean(RBE (BEG_OBS:END_OBS));
% datadescr.data_rBE   = 'Interest rate on loans to firms';
% 
% dataseries.data_rD   = RD  (BEG_OBS:END_OBS) - mean(RD  (BEG_OBS:END_OBS));
% datadescr.data_rD    = 'Interest rate on deposits';
% 
% dataseries.data_rIB  = EONIA (BEG_OBS:END_OBS) - mean(EONIA (BEG_OBS:END_OBS));
% datadescr.data_rIB   = 'Short-term interest rate';
% 
% dataseries.data_D    = DR_HP  (BEG_OBS:END_OBS);
% datadescr.data_D     = 'Real deposits';
% 
% dataseries.data_BE   = BER_HP (BEG_OBS:END_OBS);
% datadescr.data_BE    = 'Real loans to firms';
% 
% %dataseries.data_R_BU = R_BU (BEG_OBS:END_OBS)  - mean(R_BU (BEG_OBS:END_OBS));
% %datadescr.data_R_BU   = 'Interest rate on unsecured loans';
% 
% dataseries.data_PIE  = PIE (BEG_OBS:END_OBS) - mean(PIE (BEG_OBS:END_OBS));
% datadescr.data_PIE   = 'Inflation (q/q)';



% saving the renamed series for Dynare
save ('data/DYNARE_estimdata','-STRUCT','dataseries')


estimated_params ;

%                 Initial conditions    PRIOR shape     MEAN     STD
%stderr e_z       ,       0.010,       inv_gamma_pdf,   0.0100,   0.05  ;      
stderr e_A_e      ,       0.007,       inv_gamma_pdf,   0.0100,   0.05  ;  
%stderr e_j       ,       0.020,       inv_gamma_pdf,   0.0100,   0.05  ;  
stderr e_m_e      ,      0.0030,       inv_gamma_pdf,   0.0100,   0.05  ;  
%stderr e_mi      ,      0.0041,       inv_gamma_pdf,   0.0100,   0.05  ;  
stderr e_mk_d     ,        0.09,       inv_gamma_pdf,   0.0100,   0.05  ;  
stderr e_mk_b_eeu   ,      0.0049,       inv_gamma_pdf,   0.0100,   0.05  ;  
stderr e_mk_b_ee    ,      0.0835,       inv_gamma_pdf,   0.0100,   0.05  ;	//10
stderr e_qk      ,      0.0136,       inv_gamma_pdf,   0.0100,   0.05  ;      
stderr e_r_ib    ,      0.0020,       inv_gamma_pdf,   0.0100,   0.05  ;  
stderr e_mk_y    ,        0.90,       inv_gamma_pdf,   0.0100,   0.05  ;	
stderr e_l       ,       0.348,       inv_gamma_pdf,   0.0100,   0.05  ;	
stderr e_eps_K_b ,       0.008,       inv_gamma_pdf,   0.0100,   0.05  ;	
%rho_ee_z         ,        0.55,            beta_pdf,     0.80,   0.10  ;                   
rho_A_e          ,        0.95,            beta_pdf,     0.80,   0.10  ;                   
%rho_ee_j         ,        0.98,            beta_pdf,     0.80,   0.10  ;
rho_m_e           ,        0.95,            beta_pdf,     0.80,   0.10  ;
%rho_mi           ,        0.95,            beta_pdf,     0.80,   0.10  ;                   
rho_mk_d          ,        0.95,            beta_pdf,     0.80,   0.10  ;                   
rho_mk_b_eeu      ,        0.85,            beta_pdf,     0.80,   0.10  ;                   
rho_mk_b_ee       ,        0.85,            beta_pdf,     0.80,   0.10  ;                   
rho_ee_qk         ,        0.55,            beta_pdf,     0.80,   0.10  ;                   
rho_mk_y          ,        0.80,            beta_pdf,     0.80,   0.10  ;                   
rho_eps_l        ,        0.60,            beta_pdf,     0.80,   0.10  ;                   
rho_eps_K_b      ,        0.80,            beta_pdf,     0.80,   0.10  ;     
kappa_p 		 ,          28.65,           gamma_pdf,       50,    20   ;   	
kappa_u          ,          10,           gamma_pdf,     10,    5    ;
%kappa_w 		 ,          50,           gamma_pdf,       50,    20   ;			
kappa_i 		 ,         5,           gamma_pdf,      2.5,    1.0  ;			
kappa_d 		 ,          110,           gamma_pdf,       100,    50  ;
kappa_b_ee		 ,          50,           gamma_pdf,        50,    25  ;
kappa_b_eeu		 ,           60,           gamma_pdf,        60,    30  ;
kappa_kb         ,           10,           gamma_pdf,     10.0,    5.0  ;
phi_pie          ,         2.0,           gamma_pdf,      2.0,    0.5  ; 
rho_ib           ,        0.77,            beta_pdf,     0.75,    0.10 ;   
phi_y            ,         0.0,          normal_pdf,     0.10,    0.15 ;
gamma            ,        1.54,         inv_gamma_pdf,   1.5,     0.25;
omega_r          ,        0.682,        inv_gamma_pdf,    0.5,    0.25;
%ind_p            ,        0.20,            beta_pdf,     0.50,    0.15 ;
%ind_w            ,        0.25,            beta_pdf,     0.50,    0.15 ;
%a_i              ,         0.6,            beta_pdf,     0.50,    0.10 ; 

end;

estimation(datafile='data/DYNARE_estimdata'
                ,mode_compute = 9 %9
                ,mode_file    = 'data/final_exp_mh_mode'
                ,mh_jscale    = 0.30            
                ,presample    = 1
                ,prefilter    = 0
                ,prior_trunc  = 1e-14
                ,mh_replic    = 100000 %put 500 for test and a big number (like 100000) when serious
                ,mh_nblocks   = 5   %put 1   for test and 5 or 10 when serious
                ,filtered_vars
                ,lik_init     = 1
                ,order        = 1
                ,mode_check
                )                
                C I Y1 r_ib r_b r_d pie B D J_B K_b;
             
save([M_.fname '_results.mat'],'oo_','M_','estim_params_','options_','dataseries','datadescr');

@#endif

%------------------------------------Welfare analysis----------------------------------------------------% 

%Make sure Dynare does not print out stuff during runs
% options_.nocorr=1;
% options_.noprint=1;
% options_.nofunctions=1;
% options_.nograph=1;
% options_.verbosity=0;
%stoch_simul(order=2,irf=0, nofunctions);

@#if occbin == 0

%rho_vi phi_vi
@#if policy_regime == 0

    %x_opt_name={'rho_vi',0,Inf
     %   'phi_vi',0,Inf
      %  };
    options_.policy_regime=0;
    x_start=[0.8,2.35]; %use optimal values as starting point
    A = []; %arguments for fmincon
    b = []; %arguments for fmincon
    Aeq = []; %arguments for fmincon
    beq = []; %arguments for fmincon
    lb = [0,0]; %sets lower and upper bounds for policy rule coefficients
    ub = [1,5];
    [xhat,fhat] = fmincon(@welfare_objective,x_start,A,b,Aeq,beq,lb,ub); %calls separate file to solve for optimal policy with declared arguments
    %structure above returns argument xhat and value fhat of the function f(x)
    nk_gs14_xopt_fmc_0=xhat;
    save ('nk_gs14_xopt_fmc_0','nk_gs14_xopt_fmc_0'); %saves optimal policy coefficients
    @#endif


%rho_vi phi_vi
@#if policy_regime == 1
    options_.policy_regime=1;
    x_start=[0.8,2.35]; %use optimal values as starting point
    A = []; %arguments for fmincon
    b = []; %arguments for fmincon
    Aeq = []; %arguments for fmincon
    beq = []; %arguments for fmincon
    lb = [0,0]; %sets lower and upper bounds for policy rule coefficients
    ub = [1,5];
    [xhat,fhat] = fmincon(@welfare_objective,x_start,A,b,Aeq,beq,lb,ub); %calls separate file to solve for optimal policy with declared arguments
    %structure above returns argument xhat and value fhat of the function f(x)
    nk_gs14_xopt_fmc_1=xhat;
    save ('nk_gs14_xopt_fmc_1','nk_gs14_xopt_fmc_1'); %saves optimal policy coefficients
    @#endif

    @#endif