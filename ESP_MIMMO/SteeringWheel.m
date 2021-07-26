function delta = SteeringWheel(t)

ti = 5; % [s] inital time
tf = 7; % [s] final time
A = 5/180*pi; %[rad] steering wheel maximum amplitude

type = 0;
% 0 = sin
% 1 = step

switch type
    case 0
        % SIN COMPUTATION
        DT = tf-ti; % [s] steering action duration
        w = 2*pi/DT;
        
        if t >= ti && t <= tf
            delta_f = A*sin(w*(t-ti));
        else
            delta_f = 0;
        end
    case 1
        if t >= ti 
            delta_f = A;
        else
            delta_f = 0;
        end
end

delta = delta_f*[1;1;0;0];