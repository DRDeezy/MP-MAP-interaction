%Loads alternative welfare improving rules

load 'MAP_Chat_1.mat'
load 'LTV_Chat_2.mat'
load 'weights_Chat_3.mat'
load 'law_borr_Chat_4.mat'
load 'law_asset_prices_Chat_5.mat'

%Loads initial values for the baseline model (basic rule)

load 'MP_only_Chat_0.mat'


%Calculates changes in consumption

change_1 = (MAP_Chat_1-MP_only_Chat_0)./MP_only_Chat_0;
change_2 = (LTV_Chat_2-MP_only_Chat_0)./MP_only_Chat_0;
change_3 = (weights_Chat_3-MP_only_Chat_0)./MP_only_Chat_0;
change_4 = (law_borr_Chat_4-MP_only_Chat_0)./MP_only_Chat_0;
change_5 = (law_asset_prices_Chat_5-MP_only_Chat_0)./MP_only_Chat_0;

n = 1;
lambd = 0.5;

%Calculates changes in welfare

change_EMU_1 = n*(lambd*change_1(1) + (1-lambd)*change_1(2));
change_EMU_2 = n*(lambd*change_2(1) + (1-lambd)*change_2(2));
change_EMU_3 = n*(lambd*change_3(1) + (1-lambd)*change_3(2));
change_EMU_4 = n*(lambd*change_4(1) + (1-lambd)*change_4(2));
change_EMU_5 = n*(lambd*change_5(1) + (1-lambd)*change_5(2));

%Stores pairs of consumption change / changes in welfare 
%From zero to 1
welfare_result_1 = [change_1, change_EMU_1];
welfare_result_2 = [change_2, change_EMU_2];
welfare_result_3 = [change_3, change_EMU_3];
welfare_result_4 = [change_4, change_EMU_4];
welfare_result_5 = [change_5, change_EMU_5];


%Creates and saves table
%1000 is a scaling factor
welfare_results = 1000 * [welfare_result_1;welfare_result_2;welfare_result_3;welfare_result_4;welfare_result_5];
rounded_results = round(welfare_results,5);
welfare_results_table = array2table(rounded_results, 'VariableNames',{'Savers','Borrowers','Weighted'}, ...
    'RowNames',{'MAP','LTV', 'Weights', 'LAW_B', 'LAW_P'} );
%save('welfare_result','welfare_result');
save('welfare_results_table', 'welfare_results_table')
welfare_results_table