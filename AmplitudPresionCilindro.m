close all
clear variables
clc

tau=10*(10^-9);         %Tau
deltat=500*(10^-12);    %periodo de muestreo
c=1500;                 %velocidad del sonido
ROl=0.0015875;
NpS=400;
NpF=150;
distCC=0;
RF=0.003/2;


%-------------------------- sección alpha 1300-------------------------
LongCil=3*(10^-3);
RebZ=ceil((LongCil*11)/(60*(10^-6)));
Alpha=766;

ZDisc=0.045;
tic;
[ t,P ] = CalcPresSensDiscXYv1DAFuenteCilindroLB(ROl,NpS,NpF,distCC,RF,Alpha,ZDisc,LongCil,RebZ,deltat,tau,c);
tC1=toc;
filename = 'P766.mat';
save(filename,'t','P')

clear t P

Alpha=233;
ZDisc=0.045;
tic;
[ t,P ] = CalcPresSensDiscXYv1DAFuenteCilindroLB(ROl,NpS,NpF,distCC,RF,Alpha,ZDisc,LongCil,RebZ,deltat,tau,c);
tC2=toc;
filename = 'P233.mat';
save(filename,'t','P')
clear t P

% 
% ZDisc=0.035;
% tic;
% [ t3,P3 ] = CalcPresSensDiscXYv1DAFuenteCilindroLBMat(ROl,NpS,NpF,distCC,RF,Alpha,ZDisc,LongCil,RebZ,deltat,tau,c);
% tC3=toc;
% filename = 'Pt3.mat';
% save(filename,'t3','P3')
% clear t3 P3
% 
% ZDisc=0.03;
% tic;
% [ t4,P4 ] = CalcPresSensDiscXYv1DAFuenteCilindroLBMat(ROl,NpS,NpF,distCC,RF,Alpha,ZDisc,LongCil,RebZ,deltat,tau,c);
% tC4=toc;
% filename = 'Pt4.mat';
% save(filename,'t4','P4')
% clear t4 P4
% 
% %-------------------------- sección alpha 25500-------------------------
% % LongCil=3*(10^-3);
% LongCil=0.00015;
% RebZ=ceil((LongCil*11)/(60*(10^-6)));
% Alpha=25500;
% 
% ZDisc=0.045;
% tic;
% [ t5,P5 ] = CalcPresSensDiscXYv1DAFuenteCilindroLBMat(ROl,NpS,NpF,distCC,RF,Alpha,ZDisc,LongCil,RebZ,deltat,tau,c);
% tC5=toc;
% filename = 'Pt5.mat';
% save(filename,'t5','P5')
% clear t5 P5
% 
% 
% ZDisc=0.04;
% tic;
% [ t6,P6 ] = CalcPresSensDiscXYv1DAFuenteCilindroLBMat(ROl,NpS,NpF,distCC,RF,Alpha,ZDisc,LongCil,RebZ,deltat,tau,c);
% tC6=toc;
% filename = 'Pt6.mat';
% save(filename,'t6','P6')
% clear t6 P6
% 
% 
% ZDisc=0.035;
% tic;
% [ t7,P7 ] = CalcPresSensDiscXYv1DAFuenteCilindroLBMat(ROl,NpS,NpF,distCC,RF,Alpha,ZDisc,LongCil,RebZ,deltat,tau,c);
% tC7=toc;
% filename = 'Pt7.mat';
% save(filename,'t7','P7')
% clear t7 P7
% 
% ZDisc=0.03;
% tic;
% [ t8,P8 ] = CalcPresSensDiscXYv1DAFuenteCilindroLBMat(ROl,NpS,NpF,distCC,RF,Alpha,ZDisc,LongCil,RebZ,deltat,tau,c);
% tC8=toc;
% filename = 'Pt8.mat';
% save(filename,'t8','P8')
% clear t8 P8

filename = 'WorkS.mat';
save(filename)

% %--------------------------------------------------------------------
% 
% c=1500;                 %velocidad del sonido
% d=(3*(.1^2))^(1/2);     %Distancia máxima
% % d=0.4;
% % d=0.25;
% deltat=500*(10^-12);    %periodo de muestreo
% n=d/c;  
% t=0:deltat:n-deltat;    %vector tiempo
% 
% % función anónima, -------------------------------------------------------
% %Respuesta al impulso del sensor usado
% Olympus=@(x)(...
%     (0.1 *exp(-0.5*(((x-16.5e-9)./0.9e-9).^2)))+(-0.5 *exp(-0.5*(((x-19e-9)./0.9e-9).^2)))+ ...
%     (0.6 *exp(-0.5*(((x-22e-9)./0.9e-9).^2)))+(-0.3 *exp(-0.5*(((x-25e-9)./0.9e-9).^2)))+...
%     (0.17*exp(-0.5*(((x-28e-9)./0.9e-9).^2)))+(-0.13*exp(-0.5*(((x-31e-9)./0.9e-9).^2)))+...
%     (0.05*exp(-0.5*(((x-34e-9)./0.9e-9).^2))));
% %calculando ek kernel con la respuesta al impulso
% % tk=0:1e-9:50e-9;
% tk=t(1:100);
% RIO=Olympus(tk);
% 
% VPuntualRIO=conv(P2,RIO,'same');
% figure();
% plot(ti2,VPuntualRIO)