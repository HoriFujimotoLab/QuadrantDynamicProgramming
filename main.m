Xmax = 200;         % Maximum inter-vehicular distance [m]
Vmax = 120/3.6;     % Maximum velocity of ego vehicle [m/s]

Xnum = 100 + 1;     % Count of discritized distance
Vnum = 120 + 1;     % Count of discritized velocity
Pnum = 60  + 1;     % Count of preceding velocity

[QDPData,~] = qdp(@qdp_cost,Xmax,Vmax,Xnum,Vnum,Pnum);

save("QDPData", "QDPData");


function cost = qdp_cost(time,velocity,acceleration,distance)
% QDP Cost - Cost function for QDP
% 
% cost = qdp_cost(time,velocity,acceleration,distance)
%
% Inputs and the output are vectors of same size
% to avoid for loop and enhance computation time. 

% Dummy cost function
cost = time .* velocity + acceleration.^2 + distance;
% Use .* ./ .^ to calculate the cost element by element.

end