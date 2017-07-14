function [ measurement_prediction,num_of_data ] = measurement_predction_func2nd( measurement_history )

measurement_history(isnan(measurement_history)) = [];
num_of_data = 0;
if isempty(measurement_history)
    measurement_prediction = NaN;
else
    
t = 0:1:length(measurement_history)-1;

num_of_data = length(t);

%% polyfit
p = polyfit(t,measurement_history,2);

%% One sample measurement prediction
measurement_prediction = p(3)+p(2)*(t(end)+1)+p(1)*(t(end)+1)^2;

end

