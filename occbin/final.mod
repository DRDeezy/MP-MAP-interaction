% Final model for the dissertation

%Setup Dynare launch options

%Need to tune initial values for phi_B and phi_AP, depending on whether households have profits from banks or not

%OccBin setup
%------------
%Turn on OccBin; enable graph plotting; set sign of shocks; =1 enables
@#define occbin = 1
@#define graph_plotting = 0
@#define positive_shock = 0

%Estimated parameters
%------------
%Use parameters obtained after estimation; =1 enables
@#define estimated_param = 0

%Welfare setup
%------------
% = 1 Enables welfare analysis
@#define welfare_analysis = 0
% = 1 Enables cooperative policies
@#define cooperative = 0

%Readme:
%%%%%
%Set different types of shocks
% 1. TFP shock (A_e)
% 2. Cost push shock (mk_y)
% 3. Shock to LTV ratio
% 4. Three shocks bank mark-ups (a la Gerali et al.)
% 5. Shock to banking capital
%%%%%
@#define policy_regime = 2
%Setting the active policy
% = 0 Monetary Policy (Classic Taylor rule) only; not used in cooperative
% = 1 MP + Time varying capital requirement MAP
% = 2 MP + Time varying LTV cap
% = 3 MP + Time varying risk weights
% = 4 Lean against the wind MP: total borrowing
% = 5 Lean against the wind MP: asset prices
%%%%%


%---------------------Endogenous variables---------------------------------------------%

var 

W_EMU ${W^{EMU}}$ (long_name='Welfare deviation to steady state in the economy')    %1
Y1         ${Y1}$      (long_name='Definition of output used in policy rules')      %2
pie   ${\pi}$ (long_name='Inflation rate')                                          %3
B_rw       ${B_{rw}}$  (long_name='Risk-weighted loans') %4 
r_ib  ${r_{ib}}$ (long_name='Policy rate')   %5
vi     ${\nu}$ (long_name='Macroprudential policy instrument') %6

%-------------------------Log-deviations-------------------------%
log_r_ib

@#if occbin == 1
%----------------------------OccBin block-----------------------------------------------%
log_b_ee
log_q_k
log_k_e
one_plus_r_b
log_one_plus_r_b
q_k_exp ${q_k_exp}$ (long_name='Auxilary varaible')
log_q_k_exp
log_c_e %8
@#endif


%-----------------------Macropudential Policies---------------------------------------------%

%Declaring variables
@#if policy_regime == 0 || policy_regime == 1 || policy_regime == 2 || policy_regime == 3 || policy_regime == 4 || policy_regime == 5 
m_e    ${LTV^E}$ (long_name='Time-varying LTV cap')
%vi     ${\nu}$ (long_name='Macroprudential policy instrument')
weight_b   ${w_{b}}$   (long_name='Risk weight on secured loans') 
weight_bu  ${w_{bu}}$  (long_name='Risk weight on unsecured loans') 
@#endif


%----------------------------Modelling defaults-----------------------------------------------%

omega ${\Omega}$ (long_name='Macro-variable of economy indebtness')

%----------------------------Welfare-----------------------------------------------%

%W_EMU ${W^{EMU}}$ (long_name='Welfare deviation to steady state in the economy')

%------------------Utility and Welfare-------------------------------------------%

V_S ${V^S}$ (long_name='Recursive utility of savers in the core')
V_B ${V^B}$ (long_name='Recursive utility of borrowers in the core')
W_S ${W^S}$ (long_name='Welfare deviation to steady state of savers in the core')
W_B ${W^B}$ (long_name='Welfare deviation to steady state of borrowers in the core')

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

%pie   ${\pi}$ (long_name='Inflation rate')   
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

%r_ib  ${r_{ib}}$ (long_name='Policy rate')   
 
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
%Y1         ${Y1}$      (long_name='Definition of output used in policy rules')  %//Auxiliary variable 1
%B_rw       ${B_{rw}}$  (long_name='Risk-weighted loans')  
B_eeu      ${B_{eeu}}$ (long_name='Unsecured loans in bank')   
B_ee       ${B_{ee}}$  (long_name='Secured loans in bank') 
 


 %---------------------Exogenous processes---------------------------------------------%

mk_y  ${mk_y}$ (long_name='Retailers mark-up')
A_e   ${A_e}$ (long_name='Total factor productivity') 
eps_K_b

;     

%---------------------Exogenous variables---------------------------------------------%
varexo  e_A_e e_mk_y e_eps_K_b;% e_m_e;

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
            m_e_ss ${m^e_{ss}}$ (long_name='LTV cap')  
            gamma_p ${\gamma_p}$ (long_name='Relative HH weight')
            gamma_e ${\gamma_e}$ (long_name='Relative entrepreneurs weight')
            %---------------------Banks---------------------------------------------%
            vi_ss ${\nu_{ss}}$ (long_name='Capital requirement in SS')  
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
            %---------------------Shocks---------------------------------------------%        
            rho_A_e  ${\rho^A_e}$ (long_name='AR(1) coefficient of TFP process')
            rho_mk_y ${\rho^{mk}_y}$ (long_name='AR(1) coefficient of mark-up process') 
            rho_eps_K_b

            %---------------------Macroprudential policy---------------------------------------------%
            @#if policy_regime == 0 || policy_regime == 1 || policy_regime == 2 || policy_regime == 3 || policy_regime == 4 || policy_regime == 5  
            rho_m_e    ${\rho_m}$ (long_name='LTV cap inertia')
            rho_vi ${\rho_vi}$ (long_name='MAP inertia')
            phi_vi ${\phi_vi}$ (long_name='MAP output gap inertia')
            rho_w_bu ${\rho^{bu}_w}$ (long_name='Unsecured risk weight inertia')            %Estimate/calibrate
            chi_w_bu ${\chu^{bu}_w}$ (long_name='Unsecured risk weight output gap inertia') %Estimate/calibrate           
            rho_w_b  ${\rho^{b}_w}$ (long_name='Secured risk weight inertia')               %Estimate/calibrate
            chi_w_b  ${\chu^{b}_w}$ (long_name='Secured risk weight output gap inertia')    %Estimate/calibrate
            phi_AP  ${\phi^{AP}}$ (long_name='Asset prices augmented MP: asset prices gap')  
            phi_B   ${\phi^{B}}$ (long_name='Borrowings augmented MP: borrowing gap inertia')
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

%-----------Starting values of parameters------------%
@#if estimated_param == 0
%omega_ss = 311.0046484546306;
omega_ss = 3.130498674553582e+02;
%Parameters related to the main model block           
beta_p   = 0.996;
beta_e   = 0.975;
phi      = 1;           % Robustness: High: 1.5, Low: 0.5
m_e_ss   = 0.35;        % Hi-LTV case: m_e_ss = 0.70
gamma_p  = 1;
gamma_e  = 1 ;
eps_y    = 6; 
mk_y_ss  = eps_y   / (eps_y  - 1);
ksi      = 0.20;
kappa_p  =  28.65;      % Robustness: High: 100, Low: 15
kappa_i  = 5;           % For Tech. Shock. For cost-push shock, subst with kappa_i = 0.05
                        % Robustness: Tech.Sh.: kappa_i = 0.05, 0.5, 10 
                        % Cost push shock: kappa_i = 5, 0.5 0.005
deltak         = 0.050;
piss           = 1.00;
vi_ss          = 0.09; % Robustness: High: 0.15, Low: 0.045
delta_b        =  0.0721; % Calibrated such that the multiplier on b_ee is positive, hence b_ee_u is positive/non-negative
rho_A_e        = 0.95;    % Robustness: rho_A_e  = 0, rho_A_e  = 0.50
rho_mk_y       = 0.50;    % Robustness: rho_mk_y = 0, rho_mk_y = 0.75

rho_eps_K_b = 0.2070;

r_bu_ss        = (1/beta_e) - 1; 
ind_d          = 0.0;              
ind_be         = 0.0;              
ind_bh         = 0.0;                           
eps_d          = -1.46025;  % Might be estimated/calibrated                       
eps_b_ee       = 3.154542134; % Might be estimated/calibrated                       
eps_b_eeu      = 2.117626240126008; % Might be estimated/calibrated                       
r_ib_ss        = (piss/beta_p - 1) * (1-eps_d)/(-eps_d) ;

%---------------------------Costs of adjustments-----------------------%
% To be estimated/calibrated
@#if policy_regime == 0 || policy_regime == 1 || policy_regime == 3 || policy_regime == 4 || policy_regime == 5
kappa_u   = 0.5; % For MAP cap req rule, MP only rule
@#endif
@#if policy_regime == 2 
kappa_u  = 0.2; % For LTV rule
@#endif


% Old parameters, from Gerali et al.
% kappa_d      = 110;              
% kappa_b_ee   = 50;             
% kappa_b_eeu  = 60;            
% kappa_kb     = 10;  
kappa_d        = 110;              
kappa_b_ee     = 50;             
kappa_b_eeu    = 10; % To be estimated/calibrated            
kappa_kb       = 10; 

%-------------------------------------------------------Policies-----------------------------------------------------------------------%

% Macro-variable rule
omega_r  = 0.682;
gamma    = 1.54;

% Monetary Policy rule
%rho_ib   = 0.5;           % Robustness: rho_ib = 0
%phi_pie  = 25.5;          % MP RULE parameter, to be changed for simulations
%Parameters from MP optimization
% rho_ib   = 0.5457;         % Robustness: rho_ib = 0
% phi_pie  = 4.4470;         % MP RULE parameter, to be changed for simulations
% phi_y    = 0.221;              % MP RULE parameter, to be changed for simulationsv
@#endif


%Parameters related to MAP policy
@#if policy_regime == 0
rho_ib   = 0.8166;         % Robustness: rho_ib = 0
phi_pie  = 1.9750;         % MP RULE parameter, to be changed for simulations
phi_y    = 0.0;              % MP RULE parameter, to be changed for simulationsv
rho_m_e     =   0       ;
%Parameters related to MAP
rho_vi      =   0    ; 
phi_vi      =   0    ;
%Parameters related to Risk weights
rho_w_bu = 0.0;
chi_w_bu = 0.0;
rho_w_b  = 0.0;
chi_w_b  = 0.0;
%Parameters related to LAW MP
phi_AP   = 0.0;           % Asset prices
phi_B    = 0.0;           % Indebtness
@#endif

%Parameters related to MAP policy
@#if policy_regime == 1
rho_ib   = 0.9376;         % Robustness: rho_ib = 0
phi_pie  = 1.7799;         % MP RULE parameter, to be changed for simulations
phi_y    = 0.0;              % MP RULE parameter, to be changed for simulationsv
rho_m_e     =   0       ;
%Parameters related to MAP
rho_vi      =   0.9248     ; 
phi_vi      =   2.2138    ;
%Parameters related to Risk weights
rho_w_bu = 0.0;
chi_w_bu = 0.0;
rho_w_b  = 0.0;
chi_w_b  = 0.0;
%Parameters related to LAW MP
phi_AP   = 0.0;           % Asset prices
phi_B    = 0.0;           % Indebtness
@#endif


%Parameters related to LTV policy
@#if policy_regime == 2
rho_ib   = 0.6533;         % Robustness: rho_ib = 0
phi_pie  = 1.6552;         % MP RULE parameter, to be changed for simulations
phi_y    = 0.0;              % MP RULE parameter, to be changed for simulationsv
rho_m_e     = -2.7191;
%rho_m_e     = -0.7191;
rho_vi      =   0    ; 
phi_vi      =   0    ;
%Parameters related to risk weights
rho_w_bu = 0.0;
chi_w_bu = 0.0;
rho_w_b  = 0.0;
chi_w_b  = 0.0;
%Parameters related to LAW MP
phi_AP   = 0.0;           % Asset prices
phi_B    = 0.0;           % Indebtness
@#endif

%Parameters related to risk weights policy
@#if policy_regime == 3
rho_ib   = 0.7313;         % Robustness: rho_ib = 0
phi_pie  = 2.2695;         % MP RULE parameter, to be changed for simulations
phi_y    = 0.0;              % MP RULE parameter, to be changed for simulationsv
rho_m_e     =   0   ;
rho_vi      =   0    ; 
phi_vi      =   0    ;
%Parameters related to risk weights
rho_w_bu = 0.9789;
chi_w_bu = -10.0997;
rho_w_b  = 0.0001;
chi_w_b  = -19.9999;
%Parameters related to LAW MP
phi_AP   = 0.0;           % Asset prices
phi_B    = 0.0;           % Indebtness
@#endif

%Parameters related to borrowing augmented MP
@#if policy_regime == 4
rho_ib   = 0.9327;         % Robustness: rho_ib = 0
phi_pie  = 1.4076;         % MP RULE parameter, to be changed for simulations
phi_y    = 0.0;              % MP RULE parameter, to be changed for simulationsv    
rho_m_e     =   0    ;
rho_vi      =   0    ; 
phi_vi      =   0    ;
%Parameters related to risk weights
rho_w_bu = 0.0;
chi_w_bu = 0.0;
rho_w_b  = 0.0;
chi_w_b  = 0.0;
%Parameters related to LAW MP
phi_B    = -0.3165;           % Indebtness, less than 1
phi_AP   = 0.0;            % Asset prices
@#endif

%Parameters related to borrowing augmented MP
@#if policy_regime == 5
rho_ib   = 0.9777;         % Robustness: rho_ib = 0
phi_pie  = 1.6931;         % MP RULE parameter, to be changed for simulations
phi_y    = 0.0;              % MP RULE parameter, to be changed for simulationsv
rho_m_e     =   0   ;
rho_vi      =   0    ; 
phi_vi      =   0    ;
%Parameters related to risk weights
rho_w_bu = 0.0;
chi_w_bu = 0.0;
rho_w_b  = 0.0;
chi_w_b  = 0.0;
%Parameters related to LAW MP
phi_B    = 0.0;             % Indebtness
phi_AP   = 4.9844;            % Asset prices
@#endif


%--------------------------------Model Block-----------------------------------------------------------------------------------------------------------------%

model;

%----------------------------OccBin block-----------------------------------------------%
@#if occbin == 1 
log_b_ee = log(b_ee);
log_q_k = log(q_k);
log_k_e = log(k_e);
one_plus_r_b = 1 + r_b;
log_one_plus_r_b = log(one_plus_r_b);

log_c_e = log(c_e/steady_state(c_e));

[name = 'Auxilary variable']
q_k_exp = q_k(+1);

log_q_k_exp = log(q_k_exp);

@#endif

%--------------------------------Policies-----------------------------------------------------------------------------------------------------------------%
%Active MAP
@#if policy_regime == 1
[name = 'Capital requirements rule']
vi = vi_ss^(1 - rho_vi) * vi(-1)^rho_vi * 
             ((((B)/(Y1))/((steady_state(B))/(steady_state(Y1))))^phi_vi)^(1-rho_vi);
@#else
[name = 'Capital requirements rule']
vi = vi_ss;
@#endif

%LTV cap
@#if policy_regime == 2
[name = 'LTV cap rule']    
m_e = m_e_ss * (b_eeu/steady_state(b_eeu))^rho_m_e;
@#else
[name = 'LTV cap rule']
m_e = m_e_ss;
%m_e = m_e(-1) * (1 + rho_m_e) + (-rho_m_e)*m_e_ss + e_m_e;
@#endif

%Risk weights
@#if policy_regime == 3
[name = 'Risk weight for secured loans']
weight_b = (1 - rho_w_b) * 1 + (1 - rho_w_b) * chi_w_b * (Y - Y(-4)) + rho_w_b *  weight_b(-1);
[name = 'Risk weight for unsecured loans']
weight_bu = (1 - rho_w_bu) * 1 + (1 - rho_w_bu) * chi_w_bu * (Y - Y(-4)) + rho_w_bu *  weight_bu(-1);
@#else
[name = 'Risk weight for secured loans']
weight_b = 1;
[name = 'Risk weight for unsecured loans']
weight_bu = 1;
@#endif


%------------------------------------Adjustment costs and weights----------------------------------------------------%

%Adjustments for wholesale unit
#adj_bee = - kappa_kb * ( (K_b)/B_rw - vi ) * ((K_b/B_rw)^2) * weight_b;
#adj_beeu = - kappa_kb * ( (K_b)/B_rw - vi ) * ((K_b/B_rw)^2) * weight_bu; 
%Disutility from default for entrepreneurs
%#quadratic_cost = 0.5*omega*((1+r_bu)*b_eeu*deltae)^2; %no inflation in penalty
#quadratic_cost = 0.5*omega*(((1+r_bu)*b_eeu*deltae)^2)/pie(+1);
%Adjustments for retail units
#adj_retail_d_p = - kappa_d  * ( (r_d)/(r_d(-1)) - 1  )  * (r_d)/(r_d(-1))
+ beta_p * ( (lam_p(+1))/(lam_p) ) * kappa_d  * ( (r_d(+1))/(r_d) - ( (r_d)/(r_d(-1)))^ind_d )   * ( ((r_d(+1))/(r_d))^2 )   * ((D(+1))/(D));
#adj_retail_b_ee = - kappa_b_ee * ((r_b)/(r_b(-1)) - 1 ) * (r_b)/(r_b(-1))
   + beta_p * ( (lam_p(+1))/(lam_p) ) * kappa_b_ee * ( (r_b(+1))/(r_b) - ( (r_b)/(r_b(-1)))^ind_be ) * ( ((r_b(+1))/(r_b))^2 ) * ((B_ee(+1))/(B_ee));
#adj_retail_b_eeu = - kappa_b_eeu * ((r_bu)/(r_bu(-1)) - 1 ) * (r_bu)/(r_bu(-1)) 
+ beta_p * ( (lam_p(+1))/(lam_p) ) * kappa_b_eeu * ( (r_bu(+1))/(r_bu) - ( (r_bu)/(r_bu(-1)))^ind_bh ) * ( ((r_bu(+1))/(r_bu))^2 ) * ((B_eeu(+1))/(B_eeu));

%------------------------------------Modelling defaults----------------------------------------------------%

[name = 'Macro-variable']
%A-la Peiris specification
%omega = omega_ss * ( ( steady_state(b_eeu) * Y * (1 + steady_state(r_bu)) ) / ( steady_state(Y) * b_eeu * (1+r_bu) ) )^omega_r 
% * ( steady_state(deltae) / (deltae) )^gamma;

%Simplified specification
omega = omega_ss * (( steady_state(b_eeu)  / b_eeu(-1) )^(omega_r)) * (( steady_state(deltae) / deltae )^gamma);

%------------------------------------Households----------------------------------------------------%

[name = 'HH: Marginal utility of consumption']
(c_p)^(-1) = (lam_p);

[name = 'HH: Euler equation']
(lam_p) = beta_p*(lam_p(+1))*(1+r_d)/(pie(+1));

[name = 'HH: Labor supply decision']
(l_p)^phi = (lam_p)*(w_p);

[name = 'HH: Budget constraint'] %Added profits from banks
(c_p) + (d_p) = (w_p)*(l_p)+(1+r_d(-1))*(d_p(-1))/(pie)+(J_R)/gamma_p + J_B;
 

%------------------------------------Entrepreneurs----------------------------------------------------% 

[name = 'Ent: Marginal utility of consumption']
(c_e) ^(-1) = (lam_e); 

[name = 'Ent: Labor demand']
(w_p) = (1-ksi)*(y_e)/((l_pd)*(x)); 

[name = 'Ent: Euler equation']
(s_e) * m_e * (q_k(+1))*(1-deltak)*(pie(+1))/(1+r_b) + beta_e*(lam_e(+1))*((q_k(+1))*(1-deltak) + (r_k(+1))) 
    = (lam_e) * (q_k) ;

[name = 'Ent: Consumption Euler equation']    
(lam_e)-(s_e) = beta_e*(lam_e(+1))*(1+r_b)/(pie(+1));   

[name = 'Ent: Production function']   
(y_e) = (A_e) * ( (k_e(-1)) )^ksi * ( (l_pd) ) ^ (1-ksi);

%Borrowing constraint determined for OccBin
%OccBin block
@#if occbin == 1
[name = 'Ent: Borrowing constraint', relax='BC']
%log_b_ee - (log(m_e_ss) + log_q_k_exp + log(1-deltak) + log_k_e - log_one_plus_r_b ) = 0;
log_b_ee - (log(m_e) + log_q_k_exp + log(1-deltak) + log_k_e - log_one_plus_r_b ) = 0;
[name = 'Ent: Borrowing constraint', bind='BC']
s_e = 0;
@#else
[name = 'Ent: Borrowing constraint']
(1+r_b) * (b_ee) = m_e * (q_k(+1))  * (k_e) * (1-deltak) * pie(+1);
@#endif

[name = 'Ent: Return to capital']
(r_k) = ksi * (A_e) * (k_e(-1))^(ksi-1) * ( (l_pd) ) ^ (1-ksi) /(x); 

% Defaults related block for ENTREPRENEURS
[name = 'Ent: FOC with respect to delta_e']
%-omega * (b_eeu(-1)) * (deltae) * (1 + r_bu(-1)) + (lam_e)/pie = 0; 
-omega * (b_eeu(-1)) * (deltae) * (1 + r_bu(-1))/pie + (lam_e)/pie = 0; %if pie is included into penalty

[name = 'Ent: FOC with respect to unsecured debt']
%(lam_e) * (1 - kappa_u * ( (b_eeu) - steady_state(b_eeu) ) ) 
%+ beta_e * ( -omega(+1) * (deltae(+1))^2 * (b_eeu) * ((1 + r_bu)^2) 
%- (lam_e) * (1 - deltae(+1)) * (1 + r_bu)/(pie(+1)) )  = 0; 
(lam_e) * (1 - kappa_u * ( (b_eeu) - steady_state(b_eeu) ) ) %if pie is included into penalty
+ beta_e * ( -omega(+1) * (deltae(+1))^2 * (b_eeu) * ((1 + r_bu)^2)/pie(+1) 
- (lam_e) * (1 - deltae(+1)) * (1 + r_bu)/(pie(+1)) )  = 0; 


[name = 'Ent: Budget constraint']
(c_e) + (1+r_b(-1)) * (b_ee(-1))/(pie)  + (1-deltae) * (1 + r_bu(-1)) * (b_eeu(-1))/(pie) + (w_p)*(l_pd)  + (q_k) * (k_e) 
   = (y_e) / (x) + (b_ee) + (b_eeu) + (q_k) * (1-deltak) * (k_e(-1)) + 0.5*kappa_u*(b_eeu - steady_state(b_eeu))^2 ; 

%------------------------------------Capital Producers----------------------------------------------------% 

[name = 'Capital accumulation']
(K) = (1-deltak) * (K(-1)) + ( 1 - kappa_i/2 * ((I)/(I(-1)) - 1)^2 ) * (I) ;

[name = 'Capital producers FOC']
1 = (q_k) * ( 1 -  kappa_i/2 * ((I)/(I(-1)) - 1)^2  - kappa_i * ((I)/(I(-1)) - 1) * (I)/(I(-1)) ) 
  + beta_e * (lam_e(+1)) / (lam_e) * (q_k(+1)) *   kappa_i * ((I(+1))/(I) - 1)  * ((I(+1))/(I))^2 ;

    
%------------------------------------Retailers----------------------------------------------------% 

[name = 'NK Phillips curve']
1 - (mk_y)/((mk_y)-1)+ (mk_y)/((mk_y)-1)*(mc_E) - kappa_p*((pie) - 1)*(pie) 
   + beta_p*(lam_p(+1))/(lam_p)*kappa_p*((pie(+1))-1)*(pie(+1))*((Y(+1))/(Y)) = 0;

[name = 'mc_E definition']
(mc_E) = 1 / (x);

[name = 'Retailer profits']
J_R = y_e * ((1 - mc_E) - 0.5*kappa_p*(pie - 1)^2);

%------------------------------------Monetary Policy----------------------------------------------------% 

[name = 'Interest rate rule'] %Substitute r_ib_ss for steady_state(r_ib_ss)
 (1+r_ib) = (1+r_ib_ss)^(1-rho_ib)*(1+r_ib(-1))^rho_ib*(((pie)/piss)^phi_pie*((Y1)/(steady_state(Y1)))^phi_y  
            *((q_k)/(steady_state(q_k)))^phi_AP*((B)/(steady_state(B)))^phi_B)^(1-rho_ib);
%(1+r_ib) = (1+r_ib_ss)^(1-rho_ib)*(1+r_ib(-1))^rho_ib*(((pie)/piss)^phi_pie*((Y1)/((Y1(-1))))^phi_y  
%           *((q_k)/(steady_state(q_k)))^phi_AP*((B)/(steady_state(B)))^phi_B)^(1-rho_ib);


%------------------------------------Banks----------------------------------------------------% 

[name = 'Aggregate bank profits'] %add adjustment costs
%(J_B)  = r_b*(b_ee) + r_bu*(b_eeu)*(1-deltae) - r_d*(d_p) + adj_bee + adj_beeu;% + adj_retail_d_p + adj_retail_b_ee + adj_retail_b_eeu;

(J_B)  = r_b*(b_ee) + r_bu*(b_eeu)*(1-deltae) - r_d*(d_p)
           - kappa_d/2  * ( ((r_d)/(r_d(-1))-1)^2)   * (r_d) *(d_p) 
           - kappa_b_ee/2 * ( ((r_b)/(r_b(-1))-1)^2) * (r_b)*(b_ee) 
           - kappa_b_eeu/2 * ( ((r_bu)/(r_bu(-1))-1)^2) * (r_bu)*(b_eeu)
           - kappa_kb/2 * ( ((K_b) / (B_rw)  - vi ) ^2) * (K_b);


[name = 'Risk-weighted assets']
B_rw = weight_bu * B_eeu + weight_b * B_ee;

[name = 'Wholesale rate on secured loans']
R_b = (r_ib) + adj_bee ; 

[name = 'Wholesale rate on unsecured loans']
R_bu = (r_ib) + adj_beeu ; 

[name = 'Retail rate on unsecured loans']
(1 - deltae) * (1 - eps_b_eeu) + eps_b_eeu  * (R_bu)/(r_bu) + adj_retail_b_eeu = 0;

[name = 'Retail rate on secured loans']
+ 1 - eps_b_ee + eps_b_ee * (R_b)/(r_b) + adj_retail_b_ee = 0;

[name = 'Retail rate on deposits']
-1 + eps_d - eps_d * (r_ib)/(r_d) + adj_retail_d_p = 0; 

[name = 'Aggregate bank capital']
(K_b)*pie  = (1-delta_b)*(K_b(-1))/eps_K_b + (J_B(-1));

%[name = 'Real rate']
%rr = r_b - (pie);

[name = 'Output definition for policies']
(Y1)   = (C)+(I);

[name = 'Leverage']
(lev)  =  (B)/(K_b);

%------------------------------------Aggregation, equilibrium and market clearing--------------------------------------------% 

[name = 'Secured loan market clearing'] 
B_ee =  b_ee;
[name = 'Unsecured loan market clearing'] 
B_eeu = b_eeu;
[name = 'Aggregate loan market clearing'] 
(B)              = gamma_e*(b_ee) + gamma_e*(b_eeu);
[name = 'Aggregate consumption']
(C)              = gamma_p*(c_p) + gamma_e*(c_e); 
[name = 'Labor market clearing']
gamma_e * (l_pd) = gamma_p*(l_p); 
[name = 'Deposit market clearing']
(D)              = gamma_p*(d_p) ; 
[name = 'Capital market clearing']
(K)              = gamma_e*(k_e); 
[name = 'Bank balance sheet']
(B)              = (D) + (K_b);  
[name = 'Output clearing']
(Y)              = gamma_e*(y_e); 
%[name = 'Aggregate outpute definition']
%(Y)              = (C) + (q_k)*((K) - (1-deltak)*(K(-1))) + delta_b*(K_b(-1)) - deltae*b_eeu; %Make sure not redundant in the model
%(Y)              = (C) + (q_k)*((K) - (1-deltak)*(K(-1)));

%---------------------------------------Equations for welfare--------------------------------------------------------------------------------%

[name = 'Recursive utility of households (savers)']
V_S = (log((c_p)) - ((l_p)^(1+phi)/(1+phi))) + beta_p*V_S(+1);

[name = 'Recursive utility of entrepreneurs (borrowers)']
V_B = (log(c_e) - quadratic_cost) + beta_e*V_B(+1);

[name = 'Welfare criterion of households (savers)']
W_S = V_S -steady_state(V_S);

[name = 'Welfare criterion of entrepreneurs (borrowers)']
W_B = V_B -steady_state(V_B);

[name = 'Aggregate welfare criterion in the economy']
W_EMU = 0.5*W_S + 0.5*W_B;

%------------------------------------Log-deviations----------------------------------------------------% 

log_r_ib = log(r_ib/steady_state(r_ib));

%------------------------------------Exogenous processes----------------------------------------------------% 

[name = 'Technology shock']
(A_e)  = 1 - rho_A_e + rho_A_e*(A_e(-1)) - e_A_e;

[name = 'Bank capital shock']
(eps_K_b)  = (1-rho_eps_K_b)*    1  + rho_eps_K_b* (eps_K_b(-1)) + e_eps_K_b;%

@#if occbin == 0
[name = 'Goods mark-up shock'] %Retailers, for NON OccBin case
(mk_y) = (1-rho_mk_y)*mk_y_ss + rho_mk_y*(mk_y(-1)) + e_mk_y;
@#endif


@#if occbin == 1
    @#if positive_shock == 1
        [name = 'Goods mark-up shock']
        (mk_y) = (1-rho_mk_y)*mk_y_ss + rho_mk_y*(mk_y(-1)) + e_mk_y;
        @#endif
        @#if positive_shock == 0
            [name = 'Goods mark-up shock']
            (mk_y) = (1-rho_mk_y)*mk_y_ss + rho_mk_y*(mk_y(-1)) - e_mk_y;
            @#endif
@#endif


end;

 
%OccBin block for borrowing constraint
@#if occbin == 1
occbin_constraints;
%name 'BC'; bind s_e<0-1e8; relax log_b_ee + 0+1e8 > log(m_e_ss) + log_q_k_exp + log(1-deltak) + log_k_e - log_one_plus_r_b; 
name 'BC'; bind s_e<0-1e8; relax log_b_ee + 0+1e8 > log(m_e) + log_q_k_exp + log(1-deltak) + log_k_e - log_one_plus_r_b; 
end;
@#endif

steady(solve_algo=0);
resid;
model_diagnostics;

shocks;
var e_A_e            =  (1)^2; %Small shocks, o/w negative values for endogenous parameters OR need to reestimate parameters
%var e_A_e            =  (0.0072)^2;
var e_mk_y           =  (0.2572)^2;
%var e_eps_K_b        =  (0.7376)^2;
end;

%stoch_simul(order=1, irf=1000, periods=10000) c_e;
%----------------------------------------------------------------------------------------------------------------%
%-----------------------------------------------OccBin Block-----------------------------------------------------% 


@#if occbin == 1

    shocks(surprise); %shock to TFP; value of the shock should be larger than from estimation, as estimation assumes binding constraint
    var e_mk_y;
    periods 1:9, 10:16, 16:20, 21:25;%, 10:15, 16:21, 22:27;% 10, 50, 90, 130, 131:169;
    values 0.5, -0.5, 0.5, -0.5;%, 0.7, 0.3, 0;% -0.01,-0.02, 0.01, 0.02, 0;
    var e_A_e;
    periods 1:9, 10:16, 16:20, 21:25;%, 10:15, 16:21, 22:27;% 10, 50, 90, 130, 131:169;
    values 0.5, -0.5, -0.5, 0.5;%, 0.7, 0.3, 0;% -0.01,-0.02, 0.01, 0.02, 0;
    end;

occbin_setup;
occbin_solver(simul_periods=200,simul_check_ahead_periods=200);

%occbin_graph(noconstant) c_e log_c_e;
@#if graph_plotting == 1
    @#if positive_shock == 0
        neg_shock = oo_.occbin.simul.piecewise(1:end,8); %After 1:end set the index of a desired variable, as in var block
        @#endif
        @#if positive_shock == 1
            pos_shock = oo_.occbin.simul.piecewise(1:end,8); %After 1:end set the index of a desired variable, as in var block
            @#endif
            @#endif

@#endif


