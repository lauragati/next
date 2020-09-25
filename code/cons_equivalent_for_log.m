function ce = cons_equivalent_for_log(param,Wa,Wb)
% Consumption equivalent = how much do I have to increase consumption in model A to be indifferent
% between model A and B? (B is the alternative policy)
% For log utility only
% 24 Sept 2020
% Wa and Wb must be expected utility under models A and B
bet=param.bet;

ce = exp((1-bet)*(Wb-Wa))-1;
