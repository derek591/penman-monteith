function [RnModel, HgModel, H, HL] = PMmodel(avgPer,Kc,VWC,albedo_sat,albedo_dry,P,Rs,U,T2,RH2,Tg0,Tg5,Tg25)

%% Albedo as a function of VWC - Model from Idso et al.
albedo_model = VWC*(albedo_sat-albedo_dry)./0.2+albedo_dry;

%% Heat Capacity as a function of VWC - Model from Moene and van Dam
eta = VWC;
eta_s = 0.41;  % from Clapp and Horneberger Table 2, Loamy Sand
Cp = 1.8; % from Stull P 643 with Sand(1.24) and Rocks(2.13)
Cw = 4.186;
C_model = (1-eta_s)*Cp+eta*Cw;

%% Thermal Conductivity as a function of VWC - Model from McCumber and
% Pielkle Eq. 7 and 11 - for reference, our linear fit for the East Slope is K = 2.2*VWC+0.27
psi_s = 9.0;  % Loamy Sand Clapp and Horneberger Tble 2
b = 4.38;
psi = psi_s*(eta_s./eta).^b;  % cm
Pf = log10(psi);
Kmodel = exp(-(Pf+2.7));  % cal/s/cm/C
Kmodel = Kmodel*4.186*100;

%% Rn
% find modeled radiation values
SWinModel = Rs;
SWoutModel = -Rs*albedo_model;

% LW model
e0 = 0.6108*exp(17.27*T2./(T2+237.3)); %[kPa] saturation vapor pressure at T Allen Eq. 11
ea = e0.*RH2./100; % [kPa] actual vapor pressure Allen Eq. 54
relSWRadiation = 1;
sigma = (4.90310*10^-9); % Stefan-Boltzman [MJ/m^2/day];
LWnModelmean = sigma*((prctile(Tg0,5)+273.15)^4+(prctile(Tg0,95)+273.15)^4)/2*(0.34-0.14*nanmean(ea))*(1.35*relSWRadiation-0.35); % [MJ/m^2/day]
LWnModelmean = -LWnModelmean*1e6/86400;

sigma = 5.67e-8;  % Stefan-Boltzmann constant W/m^2/K^4
e = 0.95;
LWoutModel = -e*sigma*(Tg0+273.15).^4;

LWinModel = ones(size(LWoutModel))*(LWnModelmean-nanmean(LWoutModel));
LWnModel = LWoutModel+LWinModel;
RnModel = SWinModel+SWoutModel+LWinModel+LWoutModel;

%% Hg

% soil depths
z1 = 0.05; % 5 cm
z2 = 0.25; % 25 cm

Tstar = 0.5*(Tg0+Tg5); % Bryan recommends the average, 5cm may be more charectic of the actual layer though
% Tstar = Tg5;

% find storage
deltaHg = C_model*10^6*z1*gradient(Tstar,15*60);

% find Hg at 5cm
Hg5 = Kmodel*((z2-z1)^2*Tg0 + (z1^2-(z2-z1)^2)*Tg5 - z1^2*Tg25)./((z2-z1)*z2*z1);

% find Hg at surface
HgModel = Hg5 + deltaHg;

%% H and HL
% initialize output
HL = nan(size(U));
H = nan(size(U));

% iterate through time steps
for ii = 1:numel(H)
    
    % find Rn and Hg inputs in MJ/m^2/avgPer
    Rn = RnModel(ii)*60*60/1e6;  % J/s/m^2 * 1MJ/1e6J * (60*60)s/hr = MJ/m^2/hr
    Hg = HgModel(ii)*60*60/1e6;  % MJ/m^2/hr
    
    % find median T, RH and P from LEMs
    T = T2(ii);  % deg C
    RH = RH2(ii);
    
    U_2 = U(ii);
    
    % find constants
    Delta = 4098*0.6108*exp(17.27*T/(T+237.3))/(T+237.3)^2; % [ kPa/C] saturation slope vapor pressure curve at T - Allen Eq. 13
    cp = 1.013e-3; % [MJ/kg/C] specific heat at constant pressure
    epsilon = 0.622; % ratio of molecular weight of vapour/dry water
    lambda = 2.45; % [MJ/kg] latent heat of vaporization
    gamma = cp*P/epsilon/lambda; % [kPa/C] psychrometric constant Allen Eq. 8
    e0 = 0.6108*exp(17.27*T/(T+237.3)); %[kPa] saturation vapor pressure at T Allen Eq. 11
    ea = e0*RH/100; % [kPa] actual vapor pressure Allen Eq. 54
    
    % find reference evaporation
    E0 = (0.408*Delta*(Rn-Hg)+gamma*37/(T+273)*U_2*(e0-ea))/(Delta+gamma*(1+0.34*U_2)); % [mm/hr = kg/m^2/hr] Allen Eq. 53
    
    % find actual evaporation
    E = E0*Kc; % [mm/hr = kg/m^2/hr] Allen Eq. 56
    HL(ii) = E*lambda*10^6/(60*60);
    H(ii) =  RnModel(ii)-HgModel(ii)-HL(ii);
    
end