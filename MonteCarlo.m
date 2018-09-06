% A basic Monte Carlo simulation of light transport in
% semi-infinite tissue
% A       : Fraction of light absorbed
% R       : Fraction of light reflected
% paths   : Pathlength traveled by each detected photon (cm)
% mua     : Tissue absorption coefficient (/cm)
% mus     : Tissue scattering coefficient (/cm)
% g       : Anisotrophy factor
% rho     : Source-Detector separation (cm)
% ntissue : Refractive index of tissue
% nout    : Refractive index of medium
% Nphotons: Number of photons for simulation
%% Initialization

%-- Initialize your Tissue optical properties here --%

close all;
clear all;
clc;

A =   0   ;   % Absorbtion matrix
R =   0   ;   % Resflectance matrix

mua  = 0.1 ;  % Absobtion Coefficient
lamda_a = 1/mua;
mus =  100 ;  % Scattering Coefficient
lamda_s = 1/mus;
mut = mua + mus;
lamda_t = 1/mut;
g   =   0.91 ;   %Anisotropy Factor
rho  =  5 ;      %Source Detector Separation
ntissue = 1.4;   %Refrative Index of Tissue
nout   = 1;      %Refractive Index of Air
Nphotons = 1e5;  %No. of Photons
Reflectance=0;
R_rho=0;
radius =0;
phi_cw=0;
phi_CW_n=0;
TOF_1=0;
TOF_2=0;
DistTravel_1=[];
DistTravel_2=[];
v=3e10;         % Speed of Light in cm^2/sec

%Max thickness in cm -- Set to a large value for semi-infinite
thickness = 50;

% To account for specular reflection
Rsp = ((nout - ntissue).^2) / ((nout + ntissue).^2);
Rs = Rsp.*Nphotons; % Fraction of specularly reflected photons;

% Simulation parameters to ensure photon survival
epsilon = 0.0001;   %Thershold weight
m = 10;             % for Roulette

% for Q2
rho_two= 1:0.5:5;    % Source Detector Separation
R_CW   = zeros(size(rho_two));  %Reflectance matrix for ring rho_two

h = waitbar(0,'Progress...');

%% Simulation Loop

for n=1:Nphotons % THIS 'FOR' LOOP COVERS SIMULATION OF ALL PHOTONS
    
    % INITIALIZE PHOTON PARAMETERS
    % Initialize photon position
    x=0;
    y=0;
    z=0;
    % Initialize photon direction
    mu_x=0;
    mu_y=0;
    mu_z=1;
    
    % Initialize photon weight (HINT: Photon weight will be reduced by Rsp due to specular reflection)
    w = 1 - Rsp;
    % Initialize photon path length (if desired)
    paths = 0;
    
    % Initialize scoring parameters (if desired)
    
    %Thus begins the long march of the photon...
    
    while(w >= epsilon && z>=0) % Enter check condition for photon still in tissue and has sufficient weight
        
        % Pick Photon Step Size
        %random=rand();
        step = (- lamda_t) * log(rand);
        paths=paths+step;   %total distance traveled
        
        % Move photon
        % Update path length
        
        x = x + step*mu_x;
        y = y + step*mu_y;
        z = z + step*mu_z;
        
        % Photon Absorption
        % Calculate weight absorbed
        
        W_absorbed = w * (mua / mut);
        
        % Update absorption matrix
        
        A = A + W_absorbed;
        
        % Update photon weight
        
        w = w - W_absorbed ;
        
        % Photon Scattering
        % Pick new direction
        
        % Use Heyney Greenstein formulation to find scattering angle
        % HINT: You need cos(theta) and sin(theta)
        
        cos_theta = (1/(2*g)) *(1 + (g.^2) - (  (1 - (g.^2)) / (1 - g + (2*g*rand)) )^2  );
        sin_theta=sqrt(1-(cos_theta^2));
        
        % Find new azimuthal angle
        phi = 2 * pi * rand;
        % Find the new direction cosines using the scattering and
        % azimuthal angles
        
        
        % HINT: Make sure that you pay attention to the uz=1 case
        if (abs(mu_z)>0.999) %uses a differnt formula
            mu_x = sin_theta * cos(phi);
            mu_y = sin_theta * sin(phi);
            mu_z = sign(mu_z)*cos_theta;
        else
            d=sqrt(1-(mu_z^2));
            mu_x= (sin_theta * (((mu_x*mu_z*cos(phi))-(mu_y*sin(phi)))/d) ) + (mu_x*cos_theta);
            mu_y= (sin_theta * (((mu_y*mu_z*cos(phi))-(mu_x*sin(phi)))/d) ) + (mu_y*cos_theta);
            mu_z= ((-sin_theta)*cos(phi)*d)+(mu_z*cos_theta);
            
        end
        
        
        %  % Check photon weight - for survival
        if( w < epsilon)
            % Photon's weight is below threshold
            
            % ENTER ROULETTE -- GOOD LUCK PHOTON
            if(rand()<=(1/m))
                % PHOTON SURVIVES with updated weight
                w = m*w;
            else
                % The photon's watch has ended
                w = 0;
                
            end
        end
    end
    
    % Check if Photon is reflected...
    
    if(z< 0) % Reflected
        % Update Reflection Matrix
        R = R + w ;  %Total Reflected photons
        radius = sqrt(x^2 + y^2);
        %Q1 Checking if photons are inside the rho = 5cm
        if( radius <= rho)
            % Photon is detected, i.e., it is within detection annular ring
            % Update R(rho)
            R_rho = R_rho+w;
        end
        
        
        %% Question 2
        
        for r = 1:numel(rho_two)
            
            drho=0.1*rho_two(r);  % +/- 10% of the ring
            rho_exit_lower  = rho_two(r) - drho;
            rho_exit_higher = rho_two(r) + drho;
            
            if ((rho_exit_lower<= radius) && (radius <= rho_exit_higher))
                
                R_CW(r) = R_CW(r) + w;   %update Reflection matrix for the ring of detection
                
            end
        end
        
        % Question 3
        % for 1cm separation
        if ((1.1 >= radius) && (radius >= 0.9))
            DistTravel_1=[DistTravel_1 paths];  %Storing total path
            TOF_1=DistTravel_1./v;              %time taken to travel the distance
        end
        % for 1cm separation
        if ((2.75 >= radius) && (radius >= 2.25))
            DistTravel_2=[DistTravel_2 paths];
            TOF_2=DistTravel_2./v;
        end
    end


waitbar(n/Nphotons,h,[num2str(100*n/Nphotons, '%.2f') '%']);
end

close(h);

%% Question 2

dr = 0.1*rho_two;
phi_CW = R_CW./(4*pi*rho_two.*dr*Nphotons);  %Fluence
phi_CW_n=phi_CW/(phi_CW(1));   %Normalized Fluence
% %% Fit
   start_point = [6];
   options = optimset('MaxFunEvals',1e10);
   fun=@(params)DCWmodel(params,rho_two,phi_CW_n);
   params = fminsearch(fun,start_point,options);
   [sse, FittedCurve] = DCWmodel(params,rho_two,phi_CW_n);
   
 figure;
    subplot(2,1,1);
    plot(rho_two,FittedCurve,'r',rho_two,phi_CW_n,'b.');
    title('Spatially Resolved Reflectance ');
    xlabel('Source detector separation in cm');
    ylabel('Fluence in J/cm^2');
    legend('Fit','Fluencce');

%% Question 3
Bin_edges= linspace (1e-11,1e-9,513);
%Time resolved calculations for 1cm
figure;
h1=histogram(TOF_1,100,'Normalization','pdf');

%Time resolved calculations for 1cm
figure;
h2=histogram(TOF_2,100,'Normalization','pdf');

% y_hist = h1.Values;
% BinEdges=0:0.05e-9:5e-9;
% xbin = BinEdges;
% Xbin=((h1.BinEdges(2:end)+h1.BinEdges(1:end-1))/2);
% t=h1.BinEdges(1,2:51);
% figure;
% plot(Xbin,y_hist/max(y_hist),'b.');
% h2=histogram(TOF_2,50,'Normalization','pdf');
% y_hist = h2.Values;
% %fitting
%  start_point_TD = [0.2,70 ];
%  options_TD = optimset('MaxFunEvals',1e10);
%  fun_TD=@(params_TD)TDmodelmonte(params_TD,rho,h.BinEdges(1,2:51),ntissue,nout,h.Values);
%  params_TD = fminsearch(fun_TD,start_point_TD,options_TD);
%  [sum, FittedCurve_TD] = TDmodelmonte(params_TD,rho,h.BinEdges(1,2:51),ntissue,nout,h.Values);
% figure;
%  plot(t,y_hist,'b.');
%  t,FittedCurve_TD,'r');
% Plot Time resolved result; fit to find mu_a and mu_s'


