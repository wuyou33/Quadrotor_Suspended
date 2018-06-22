 function load_traj = get_load_traj3(t)
%% Circular Displacement  
f = 0.05;
a_x = 3;
a_y = 3;
a_z = 0;

%% Desired Load Position Generation
load_traj.xL = [a_x*(1 - cos(2*pi*f*t));
    a_y * sin(2*pi*f*t);
    a_z];
load_traj.dxL = 2*pi*[f*a_x*sin(2*pi*f*t);
    f*a_y * cos(2*pi*f*t);
    0];
load_traj.d2xL = (2*pi)^2*[f^2*a_x*cos(2*pi*f*t);
    f^2*a_y * -sin(2*pi*f*t);
    0];
load_traj.d3xL = (2*pi)^3*[f^3*a_x*-sin(2*pi*f*t);
    f^3*a_y * -cos(2*pi*f*t);
    0];
load_traj.d4xL = (2*pi)^4*[f^4*a_x*-cos(2*pi*f*t);
    f^4*a_y * sin(2*pi*f*t);
    0];
load_traj.d5xL = (2*pi)^5*[f^5*a_x*sin(2*pi*f*t);
    f^5*a_y * cos(2*pi*f*t);
    f^5*a_z * cos(2*pi*f*t)];
load_traj.d6xL = (2*pi)^6*[f^6*a_x*cos(2*pi*f*t);
    f^6*a_y * -sin(2*pi*f*t);
    0];    
end