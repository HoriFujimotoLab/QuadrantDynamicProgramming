function [QDPData, AllCosts] = qdp(costfun,Xmax,Vmax,Xnum,Vnum,Pnum)
% QDP Calculate optimal table for adaptive cruise control.
%
%   Using a method called quadrant dynamic programming (qdp),
%   the table for optimal velocity trajectory is calculated.
%   
%   [QDPData,AllCosts] = qdp(costfun,Xmax,Vmax,Xnum,Vnum,Pnum)

%% Create Four Main Vectors
x_vec = linspace(0,Xmax,Xnum);  % Vector of ego distance
v_vec = linspace(0,Vmax,Vnum);  % Vector of ego velocity
p_vec = linspace(0,Vmax,Pnum);  % Vector of preceding velocity
d_vec = p_vec * 2;              % Vector of final distance

%% Initialize Result Matrix
AllCosts = NaN(Xnum,Vnum,Pnum);
Xindices = NaN(Xnum,Vnum,Pnum); % Results of next x index
Vindices = NaN(Xnum,Vnum,Pnum); % Results of next v index

%% Start Quadrant Dynamic Programming
%parfor p_idx = 1:Pnum          % Use Parallel Computing Toolbox
for p_idx = 1:Pnum
    %% Create Final Condition
    x_fin = d_vec(p_idx);       % Final ego distance
    v_fin = p_vec(p_idx);       % Final ego velocity
    [~,x_fin_idx] = min(abs(x_vec-x_fin));  % Final index of x_vec
    [~,v_fin_idx] = min(abs(v_vec-v_fin));  % Final index of v_vec
    
    %% Create Matrix of Velocity, Time, and Accel (v_from,v_next)
    % example: v_next = Vnext(?,v_next_idx) v_from = Vfrom(v_from_idx,?)
    [Vnext, Vfrom] = meshgrid(v_vec);                               % Matrices for velocity change
    Velocity = (Vnext + Vfrom)/2;                                   % Matrix of velocity
    Time = (x_vec(2) - x_vec(1)) ./ abs(Velocity - v_fin);          % Matrix of time
    Time(Vnext > v_fin & Vfrom < v_fin & Velocity >= v_fin) = NaN;  % Distance getting shorter
    Time(Vnext < v_fin & Vfrom > v_fin & Velocity <= v_fin) = NaN;  % Distance getting longer
    Accel = (Vnext - Vfrom) ./ Time;                                % Matrix of accel
    
    %% Create Result Table for p_idx
    Cost = NaN(Xnum,Vnum);
    Xidx = NaN(Xnum,Vnum);
    Vidx = NaN(Xnum,Vnum);
    Cost(x_fin_idx,v_fin_idx) = 0;
    Xidx(x_fin_idx,v_fin_idx) = 0;
    Vidx(x_fin_idx,v_fin_idx) = 0;
    
    %% Calculate Four Quadrants
    % Far & Fast
    for i = x_fin_idx+1:1:Xnum
        for j = Vnum:-1:v_fin_idx
            [Cost(i,j),Xidx(i,j),Vidx(i,j)] = calculate_cell(costfun,i,j,i-1,v_fin_idx:Vnum,x_vec,Time,Velocity,Accel,Cost);
        end
    end
    
    % Far & Slow
    for i = Xnum-1:-1:x_fin_idx
        for j = v_fin_idx-1:-1:1
            [Cost(i,j),Xidx(i,j),Vidx(i,j)] = calculate_cell(costfun,i,j,i+1,1:Vnum,x_vec,Time,Velocity,Accel,Cost);
        end
    end
    
    % Close & Fast
    for i = x_fin_idx-1:-1:1
        for j = v_fin_idx:-1:1
            [Cost(i,j),Xidx(i,j),Vidx(i,j)] = calculate_cell(costfun,i,j,i+1,1:v_fin_idx,x_vec,Time,Velocity,Accel,Cost);
        end
    end
    
    % Close & Slow
    for i = 2:1:x_fin_idx
        for j = Vnum:-1:v_fin_idx+1
            [Cost(i,j),Xidx(i,j),Vidx(i,j)] = calculate_cell(costfun,i,j,i-1,1:Vnum,x_vec,Time,Velocity,Accel,Cost);
        end
    end
    
    % Far & Fast (Again for NaNs)
    for i = x_fin_idx+1:1:Xnum
        for j = Vnum:-1:v_fin_idx
            [cost,x_idx_next,v_idx_next] = calculate_cell(costfun,i,j,i-1,v_fin_idx:Vnum,x_vec,Time,Velocity,Accel,Cost);
            if isnan(Cost(i,j)) || cost < Cost(i,j)
                Cost(i,j) = cost;
                Xidx(i,j) = x_idx_next;
                Vidx(i,j) = v_idx_next;
            end
        end
    end
    
    %% Apply Result
    Xidx(isnan(Xidx)) = -1;
    Vidx(isnan(Vidx)) = -1;
    AllCosts(:,:,p_idx) = Cost;
    Xindices(:,:,p_idx) = Xidx;
    Vindices(:,:,p_idx) = Vidx;
end

%% Create Return Struct
QDPData.Xindices = Xindices;
QDPData.Vindices = Vindices;
QDPData.x_vec = x_vec;
QDPData.v_vec = v_vec;
QDPData.p_vec = p_vec;
QDPData.d_vec = d_vec;
end

%% Define Funtion to Calculate Cell
function [cost,x_idx_next,v_idx_next] = calculate_cell(costfun,x_idx,v_idx,x_idx_next,v_idx_next_vec,x_vec,Time,Velocity,Accel,Cost)
    distance = (x_vec(x_idx) + x_vec(x_idx_next)) / 2;

    CostDiff = costfun(                             ...
        Time(v_idx,v_idx_next_vec),                 ...
        Velocity(v_idx,v_idx_next_vec),             ...
        Accel(v_idx,v_idx_next_vec),                ...
        distance * ones(1,length(v_idx_next_vec))   ...
    );

    [cost,idx] = min(Cost(x_idx_next,v_idx_next_vec)+CostDiff);
    if isnan(cost)
        x_idx_next = NaN;
        v_idx_next = NaN;
    else
        v_idx_next = v_idx_next_vec(idx);
    end
end
