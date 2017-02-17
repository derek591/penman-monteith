% Sample PM model
clc; clearvars; close all;

% sample averaging period
avgPer = 5; % min

% specify crop coefficient
Kc = 0.25;

% specify VWC (g/g) at 5 cm. Use NaN for observations
VWC = 0.1;

% specify saturated and dry albedo, I used the mediatn +/- 2 STDs
albedo_sat = 0.16;
albedo_dry = 0.25;

% Pressure in kPa
P = 86; 

% veloctiy
U = ones(100,1);

% shortwave radiation
Rs = sin([0:pi/99:pi]') * 1000;

% ground temperature at 0, 5 and 25 cm
Tg0 = ones(100,1)*25 + sin([0:pi/99:pi]')*5;
Tg5 = ones(100,1)*25 + sin([0:pi/99:pi]')*3;
Tg25 = ones(100,1)*25 + sin([0:pi/99:pi]')*1;

% air temperature and RH at 2 m
T2 = ones(100,1)*23 + sin([0:pi/99:pi]')*7;
RH2 = ones(100,1)*40 - sin([0:pi/99:pi]')*7;

[Rn, Hg, H, HL] = PMmodel(avgPer,Kc,VWC,albedo_sat,albedo_dry,P,Rs,U,T2,RH2,Tg0,Tg5,Tg25);
figure
plot([Rn,Hg,H,HL])
legend('Rn','Hg','H','HL')
