% BME 445 Project - Example provided by Dr. Kawaji %%
 
currentDirectory = pwd;
addpath(fullfile(currentDirectory,'Support-Functions/'));

 
%%  Set up Membrane Ring                                                              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This part of code creates a 100-node ring 
Num = 100;
radius = 1.5/(2*pi);
d_theta = 2*pi/Num;
theta = d_theta*(0:Num-1);
x=radius*cos(theta);
y=radius*sin(theta);
ground = zeros(1,Num); %creates the ring at rest potential
dx = 0.015;
 
 
C_m = 1; 
R_i = 0.5*ones(1,Num);
a = 0.002*ones(1,Num);
sigma_i = 150;
sigma_e = 300;
 
 
 
 
%%  Set up Mode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This part of code divides up into four modes:
% 1) Healthy mode with two stimuli (demonstrates no reentry) 
%      Reentry is difficult to obtain!
% 2) Unilateral block
%      Change cable parameters for specific region:  How?  For loop with
%      coefficient.
% 3) Complex VF considering oxygen level
%      During VF, the oxygen level in the blood goes down.  What does this
%      do to the membrane potential?
% 4) Bonus Example: stray currents with VF
%      Stray stimulus current affects specific region 
 
 
Qop1 = input('Please enter a choice:\nHealthy          ... (press 1) \nUnidirectional   ... (press 2) \nVF model         ... (press 3) \nBonus feature    ... (press 4)\n  Please choose:  ');
defib = 0;
shock = 0;
if Qop1 == 3
    defib = input('Active Defibrillator? (1 for on, 0 for off) \nPlease choose:  ');
end
 
 
tdelay3 = -8888;
 
 
%%  Set up Initial Variables for Stimuli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I_stim = 250;
tdur = 1;
tdelay = 3;
Nstim1 = 1;
tdur3 = 0.5;
 
 
if Qop1 == 1
    T_tot = 40;
    
    S1 = 40;
 
 
% Create second Stimuli
    I_stim2 = 350;
    tdur2 = 1;
    tdelay2 = 25;
    Nstim2 = 1;
 
 
    ndt = 5000;
    
elseif Qop1 == 2 % Unidirectional block values
    Variant1 = 0;
    uni = 30;
    Nstim2 = 26;
    S1 = 30;
    I_stim = 0;
    T_tot = 30;
    I_stim2 = 350;
    tdur2 = 1;
    tdelay2 = 10.5;
    tdelay = 0;
    tdur = 1;
    ndt = 4000;
    
elseif Qop1 == 3 % VF model
    Variant1 = 0;
    uni = 30;
    Nstim2 = 15;
    S1 = 30;
    ndt = 7500;
    if defib == 0
        T_tot = 60;
    else
        T_tot = 90;
        I_defib = input('Defib current:  \nPlease choose:  ');
        tdur3 = 0.5;
        tdelay3 = 1000;
        Nstim3=1:100;
        
        shock = 0;
        problem = 0;
        count2 = 0;
        count3 = 0;
        dv = 0;
        dvpvs = 0;
        
        I_stim4 = 350;
        tdur4 = 0.5;
        tdelay4 = 75;
        Nstim4=1;
 
 
        % defib algorithm
        dtrack = 0;
        problem = 0;
    end
    
    I_stim2 = 350;
    tdur2 = 1;
    tdelay2 = 10.5;
    tdelay = 0;
    tdur = 1;
    
end
dt = 0.01;
numsteps = T_tot/dt;
 
 
 
 
 
 
%%  Set up HH model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Declaration of initial potential values.
vm = ground;
 
 
% Declaration of max conductance values.
g_L_sup = 0.3; 
g_Na_sup = 120;  
g_K_sup = 36;  
e_K = -12.26;  
e_Na = 117.56;
 
 
% Calculation of initial m, n, and h values.
m = alpha_m(vm)./(alpha_m(vm)+beta_m(vm));
n = alpha_n(vm)./(alpha_n(vm)+beta_n(vm));
h = alpha_h(vm)./(alpha_h(vm)+beta_h(vm));
 
 
% Calculation of initial conductivity and current values.
g_K = g_K_sup*n.^4;
g_Na = g_Na_sup*m.^3.*h;
g_L = g_L_sup*ones(1,Num);
 
 
% Finding e_L
e_L = (ik(vm,n,g_K,e_K) + ina(vm,m,h,g_Na,e_Na))./g_L;
 
 
V = vm(S1); 
Vmax = 0;
M = m(S1);  N = n(S1); H = h(S1); 
G_k = g_K(S1); G_na = g_Na(S1); 
 
 
I_na = ina(vm,m,h,g_Na,e_Na);
I_k = ik(vm,n,g_K,e_K);
I_l = il(vm,e_L);
 
 
I_Na = I_na;
G_Na = I_na./(vm-e_Na);
 
 
Iion = I_k + I_na + I_l;
I_S = 0;
I_s = ground;
 
 
INa1 = I_na(S1);
IK1 = I_k(S1);
IL1 = I_l(S1);
 
 
coeff = a./(2*R_i.*dx*dx);
 
 
EP = 0;
 
 
%%  Set up Matrix A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
b1 = 1;
b2 = 2;
A = [-b2, b1, zeros(1,Num-3), b1; b1, -b2, b1, zeros(1,Num-3)];
 
 
 
 
for k = 3:Num-1
    B = 0*linspace(1,Num,Num);
    B(k-1) = b1;
    B(k) = -b2;
    B(k+1) = b1;
    A = [A; B];
end
A = [A; [b1, [zeros(1,Num-3)], b1, -b2]];
 
 
if (Qop1 == 2 || Qop1 == 3) && Variant1 == 1
   for k = 1:25
       R_i(k+uni-12) = 1/(abs((k -12.5))/25+0.3)*R_i(k+uni-12);
       A(k+uni-1, k+uni-2) = abs((k -12.5))/12;
       A(k+uni-1, k+uni-1) = -1-abs((k -12.5))/12;
   end
   coeff = a./(2*R_i.*dx*dx);
elseif (Qop1 == 2 || Qop1 == 3) && Variant1 == 2
   for k = 1:25
       R_i(k+uni-12) = 1/(abs((k -12.5))/25+0.3)*R_i(k+uni-12);
   end 
   coeff = a./(2*R_i.*dx*dx);
elseif (Qop1 == 2 || Qop1 == 3) && Variant1 == 0;
   for k = 1:25
       R_i(k+uni-12) = (abs((k -12.5))/25+0.3)*R_i(k+uni-12);
       A(k+uni-12, k+uni-13) = abs((k -12.5))/12;
       A(k+uni-12, k+uni-12) = -1-abs((k -12.5))/12;
   end
   coeff = a./(2*R_i.*dx*dx);
end
%% Ventricular Fibrillation
 
 
% Oxygen level (Nondimensionalized time)
oxy_lvl = 1;  % value between 0 and 1
O = 1;        % archives oxy_lvl over time
scale_factor = 1000; 
 
 
Ext_pot = 0;
count = 0;
 
 
%%  Forward Euler's Method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
 
 
for k = 1:numsteps-1    
% Euler's method for membrane potential
    v0 = dt/C_m.*(coeff.*(A*vm')' - Iion + I_s);
    vm = vm+v0;
    V = [V; vm(S1)];
    Vmax = [Vmax; max(vm)];
% Runge-Kutta Algorithm for m
    m1 = f_m(m,vm);
    m2 = f_m(m+dt*m1/2,vm);
    m3 = f_m(m+dt*m2/2,vm);
    m4 = f_m(m+dt*m3,vm);
    m0 = (m1 + 2*m2 + 2*m3 + m4)/6;
    
% Runge-Kutta Algorithm for n
    n1 = f_n(n,vm);
    n2 = f_n(n+dt*n1/2,vm);
    n3 = f_n(n+dt*n2/2,vm);
    n4 = f_n(n+dt*n3,vm);
    n0 = (n1 + 2*n2 + 2*n3 + n4)/6;
    
% Runge-Kutta Algorithm for h
    h1 = f_h(h,vm);
    h2 = f_h(h+dt*h1/2,vm);
    h3 = f_h(h+dt*h2/2,vm);
    h4 = f_h(h+dt*h3,vm);
    h0 = (h1 + 2*h2 + 2*h3 + h4)/6;
 
 
% Solving for next timestep for m, n, h.
    m = m + dt*m0;
    n = n + dt*n0;
    h = h + dt*h0;
 
 
    M = [M; m(S1)];
    N = [N; n(S1)];
    H = [H; h(S1)];
 
 
% VF Model:    
if Qop1 == 3
    if ((max(vm)-min(vm)/100)) > 30
        count = count + 1;
    else
        count = count - 1;
        if count < 0
            count =0;
        end
    end
 
 
 
 
    if count >1000
        oxy_lvl = oxy_lvl - dt*((max(vm)-min(vm)/100)-60)/scale_factor;
    end
 
 
    if oxy_lvl > 1
        oxy_lvl = 1;
    elseif oxy_lvl < 0.8
        oxy_lvl = 0.8;
    end
    O = [O, oxy_lvl]; 
 
 
    g_Na_sup = 120*oxy_lvl;  
    g_K_sup = 36*oxy_lvl;  
end
 
 
 
 
% solving for ion currents
    g_Na = g_Na_sup*m.^3.*h;
    I_na = g_Na.*(vm-e_Na);
    G_na = [G_na; g_Na(S1)];
    INa1 = [INa1; I_na(S1)];
    
    g_K = g_K_sup*n.^4;
    I_k = g_K.*(vm-e_K);
    G_k = [G_k; g_K(S1)];
    INa1 = [INa1; I_na(S1)];
    I_l = g_L.*(vm-e_L);
    IL1 = [IL1; I_l(S1)];
 
 
    I_S = [I_S; I_s(S1)];
    I_s = ground;
    Iion = I_k + I_na + I_l;
 
 
% Draw figure
    if mod(k,20) == 0
        plot3(x,y,vm, x,y,ground,':b', 'LineWidth', 2);
        axis([-0.3 0.3 -0.3 0.3 -120 120]);
        title(sprintf('t  = %2.2f ms', k*dt));
        drawnow
    end
    
    
    waitbar(k/numsteps);
    
    % First stimulus current
    if mod(k,ndt) > tdelay/dt && mod(k,ndt) <= (tdelay+tdur)/dt
        I_s(Nstim1) = I_stim;
    elseif k > tdelay2/dt && k <= (tdelay2+tdur2)/dt            
        I_s(Nstim2) = I_stim2;
    elseif k > tdelay3/dt && k <= (tdelay3+tdur3)/dt;
        I_s(Nstim3) = I_defib;
    end
    
    % Bonus Feature Example: Addition of stray stimulus current less 
    % than 1/3 of regular stimulus current.
    
    if Qop1 == 4
        if mod(k,200) < 101 && I_s(Nstim1) == 0 
            I_s(Nstim1+ floor(25*rand(1))) = rand(1)*I_stim/3;
            I_s(Nstim1+ floor(25*rand(1))) = rand(1)*I_stim/3;
            I_s(Nstim1+ floor(25*rand(1))) = rand(1)*I_stim/3;
            I_s(Nstim1+ floor(25*rand(1))) = rand(1)*I_stim/3;
        end
    end
 
 
    % Code for Extracellular potential
    x_e =10; y_e = 0; z_e = 1; 
 
 
    M_m = ground;
    I_lump = ground;
    dist = ground;
 
 
    for k1 = 2:Num-1
        M_m(k1) = pi*a(k1)^2*sigma_i*(vm(k1-1)-2*vm(k1)+vm(k1+1))/dx^2;
    end
    M_m(Num) = pi*a(Num)^2*sigma_i*(vm(Num-1)-2*vm(Num)+vm(1))/dx^2;
    M_m(1) = pi*a(1)^2*sigma_i*(vm(Num)-2*vm(1)+vm(2))/dx^2;
 
 
    I_lump(2:99) = M_m(1:98)-M_m(3:100);
    I_lump(1)    = M_m(100)-5.*M_m(2);
    I_lump(100)  = M_m(99)-5.*M_m(1);
 
 
    
    for k1 = 1:Num
        dist(k1) = sqrt((x_e-x(k1))^2 + (y_e-y(k1))^2 + z_e^2);
    end
    
    M_div = I_lump./dist;
%     M_div = M_m./dist;
    if k1 == 1
        Ext_pvs = 0;
    else 
        Ext_pvs = Ext_pot;
    end
    Ext_pot = 1/(4*sigma_e)*sum(M_div);
    EP = [EP, Ext_pot];
 
 
    % defibrillator code
    if defib == 1 && shock ~= 1
        if EP > 0.2
            problem = 0;
            count2 = 0;
            count3 = 0;
        else 
            problem = problem+1;
            count2 = count2+1;
        end
 
 
        dvpvs = dv;
        dv = EP;
        
        if count2 > 0
            dv = Ext_pot-Ext_pvs;
            if dvpvs*dv < 0
                count3 = 0;
            else
                count3 = count3+1;
            end
        end
        if count3 < 10 && problem > 4000
            % activate defib
            shock = 1;
            tdelay3 = (k+500)*dt;
        end
    end
end
 
 
%% Draw Figure
time = dt:dt:T_tot;
 
 
figure(2);
subplot(2,1,1)
plot(time, EP);
title('Extracellular potential measured 1 cm above ring axis');
xlabel('time (ms)');
ylabel('Extracellular potential (mV)');
subplot(2,1,2)
plot(time, Vmax);
title('Maximum potential of all nodes');
xlabel('time (ms)');
ylabel('Maximum membrane potential (mV)');
 
 
pause
 
 
clf;
plot(time,M,time,N,time,H);
title(sprintf('m, n, h values at Node %2.0f over time', S1));
xlabel('time (ms)');
ylabel('m, n, and h values');
 
 
pause
 
 
clf;
plot(time,O, 'b', time, 0, ':b',  'LineWidth', 2);
title(sprintf('Oxygen fraction Node %2.0f over time', S1));
xlabel('time (ms)');
ylabel('fraction of Oxygen');
 
 
clf;
plot(time,G_k/max(G_k), 'r', time, G_na/max(G_na), 'g', time, O, 'b', 'LineWidth', 2);
title(sprintf('Normalized Conductance and O values at Node %2.0f over time', S1));
xlabel('time (ms)');
ylabel('Conductance (mS/uA^2)');
