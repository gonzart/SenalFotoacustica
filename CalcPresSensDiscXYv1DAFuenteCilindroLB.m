function [ t,P ] = CalcPresSensDiscXYv1DAFuenteCilindroLB(R,nPixS,nPixF,distCC,FR,Alpha,ZDisc,LongCil,RebZ,deltat,tau,c)
% ---------Argumentos de entrada---------
% R : Radio del sensor
% nPixS : número de arcos que discretizarán al sensor
% nPixF : número de arcos que discretizarán a la muestra

% -----------Generando el disco de fuente discretizado con arcos-----------
if FR<distCC
    RadioDeteccion=linspace(distCC-FR,distCC+FR,nPixF)';
else
    RadioDeteccion=linspace(FR/nPixF,distCC+FR,nPixF)';
%     ThetaInDegreesF = atan2d(norm(cross([0,0,0],[Fcen(1),Fcen(2),0])),dot([0,0,0],[Fcen(1),Fcen(2),0]));
end
% Pf=[cosd(ThetaInDegreesF).*discreteRF,sind(ThetaInDegreesF).*discreteRF,Fcen(3).*ones(size(discreteRF))];
Pf=zeros(2*nPixF,3);
ArcoF=zeros(nPixF,1);

% ------------------Sección cilindro-------------------------------------
cilindro=linspace(ZDisc-LongCil,ZDisc,RebZ)';
ILB=exp(Alpha.*(cilindro-ZDisc));
% ----------------------------------------------------------------

% Se genera el disco más externo-------------------
for i=1:nPixF
    [xoutF,youtF] = circcirc(0,0,RadioDeteccion(i),distCC,0,FR);
    if isnan(xoutF(1))
        Pf(i,:)=[distCC+RadioDeteccion(i),0,ZDisc];
        ArcoF(i)=2*pi*RadioDeteccion(i);
    else
        Pf(i,:)=[xoutF(1),youtF(1),ZDisc];
        %--arco fuente
%         x=(FR*FR-RadioDeteccion(i)*RadioDeteccion(i)+distCC*distCC)/(2*RadioDeteccion(i)*distCC);
        x=(RadioDeteccion(i)*RadioDeteccion(i)+distCC*distCC-FR*FR)/(2*RadioDeteccion(i)*distCC);
        alf=real(acos(x));%es un numero complejo por la precision numerica de las cotas [-1 1]
        ArcoF(i)=2*RadioDeteccion(i)*alf;
    end
end
Pf((nPixF+1:2*nPixF),:)=[Pf(1:nPixF,1),Pf(1:nPixF,1),(ZDisc-LongCil).*ones(nPixF,1)];

%-----Calculo de la ventana de tiempo, Tiempo de vuelo mínimo y máximo-----
%-----del disco más externo-----
% 1° se desplazan las coordenadas de los elemento diferenciales de la muestra al
% eje X
Pf(:,1)=sqrt(dot(Pf(:,1:2),Pf(:,1:2),2));
Pf(:,2)=0;

% -----------------------sección ventanas ajustadas------------------------

% Calculamos la distancia mínima de las particulas que estan contenidas en
% el plano del sensor y obtenemos la mínima y máxima distancia
Pfdist=Pf(:,1)<=R;
PfdistI=Pf(Pfdist,:);
MinI=min(PfdistI(:,3));
% Cálculo valor máximo
PfdistINMax=[PfdistI(:,1)+R,PfdistI(:,2),PfdistI(:,3)];
PfdistINMax=sqrt(dot(PfdistINMax,PfdistINMax,2));
MaxI=max(PfdistINMax);

% Calculamos la distancia mínima de las particulas que no estan contenidas
% el plano del sensor y obtenemos la mínima y máxima distancia
Pfdist=Pf(:,1)>R;
PfdistE=Pf(Pfdist,:);
% Cálculo valor mínimo
PfdistENMin=[PfdistE(:,1)-R,PfdistE(:,2),PfdistE(:,3)];
PfdistENMin=sqrt(dot(PfdistENMin,PfdistENMin,2));
MinE=min(PfdistENMin);
% Calculo valor máximo
PfdistENMax=[PfdistE(:,1)+R,PfdistE(:,2),PfdistE(:,3)];
PfdistENMax=sqrt(dot(PfdistENMax,PfdistENMax,2));
MaxE=max(PfdistENMax);
% calculamos los valores de tiempo mínimo y máximo
Ww=3;
if MinI>MinE
    if isempty(MinE)
        tmin=(MinI/c)-(Ww*tau);
    else
        tmin=(MinE/c)-(Ww*tau);
    end
else
    if isempty(MinI)
        tmin=(MinE/c)-(Ww*tau);
    else
        tmin=(MinI/c)-(Ww*tau);
    end
end
if MaxI>MaxE
    if isempty(MaxI)
        tmax=(MaxE/c)+(Ww*tau);
    else
        tmax=(MaxI/c)+(Ww*tau);
    end
else
    if isempty(MaxE)
        tmax=(MaxI/c)+(Ww*tau);
    else
        tmax=(MaxE/c)+(Ww*tau);
    end
end

clear Pfdist PfdistE PfdistI PfdistENMin PfdistENMax PfdistINMax
% ----------------------------------------------------------------


t=tmin:deltat:tmax;
P=zeros(size(t));

width=floor(2*tau/deltat);      %ancho en pixeles que genera la ventana
posmin=floor(tmin/deltat);     %indice que ocuparía el tiempo de vuelo 
%                                minimo en un vector completo

a=-4/(c^2*tau^2);

for d=1:RebZ
    for e=1:nPixF;

        if Pf(e,1)>R
    %         rVec=Pf(e,1)-R:R/nPix:Pf(e,1)+R;
            rVec=linspace(Pf(e,1)-R,Pf(e,1)+R,nPixS);
        else
    %         rVec=R/nPix:R/nPix:Pf(e,1)+R;% no puede comenzar en cero
            rVec=linspace(R/nPixS,(Pf(e,1)+R),nPixS);
        end
        for i=1:nPixS

            [xout,yout] = circcirc(0,0,R,Pf(e,1),0,rVec(i));    

            if isnan(xout(1))
                S=[Pf(e,1),rVec(i),0];%%%%%cambiar
                Vd=([Pf(e,1),Pf(e,2),cilindro(d)]-S);%vector distancia resultante de combinar elementos sensor-fuente
                Area=pi*rVec(i)*rVec(i);
            else
                S=[xout(1),yout(1),0];
                Vd=([Pf(e,1),Pf(e,2),cilindro(d)]-S);%vector distancia resultante de combinar elementos sensor-fuente
                %--------------Calculo diferencial de área----------------
    %             x=(rVec(i)*rVec(i)-R*R+Pf(e,1)*Pf(e,1))/(2*rVec(i)*Pf(e,1));
                x=((rVec(i)*rVec(i))-(R*R)+(Pf(e,1)*Pf(e,1)))/(2*(rVec(i))*Pf(e,1));
                alf=real(acos(x));%es un numero complejo por la precision numerica de las cotas [-1 1]
                PP=(dot([Pf(e,1),0],[xout(1),yout(1)],2)/(Pf(e,1)*R));
                Theta=real(acos(PP));
                TS=(abs(xout(1)))*(abs(yout(1)));
                if Theta>pi/2
                    PO=alf*rVec(i)*rVec(i); 
                    TO=rVec(i)*(cos(alf))*(abs(yout(1)));
                    LO=PO-TO;
                    PS=(pi-Theta)*R*R;
    %                 TS=(abs(xout(1)))*(abs(yout(1)));%sacar del if
                    LS=PS-TS;
                    Area=(pi*R*R)+LO-LS;
                else % Theta<=pi/2
                    PS=(Theta)*R*R;
    %                 TS=(abs(xout(1)))*(abs(yout(1)));
                    if alf<pi/2
                        PO=alf*rVec(i)*rVec(i);
                        TO=rVec(i)*(cos(alf))*(abs(yout(1)));
                        LO=PO-TO;
                        LS=PS-TS;
                        Area=LO+LS;
                    else % alf>=pi/2
                        PO=(pi-alf)*rVec(i)*rVec(i);
                        TO=rVec(i)*(cos(pi-alf))*(abs(yout(1)));
                        LO=PO-TO;
                        LS=PS-TS;
                        Area=(rVec(i)*rVec(i)*pi)+LS-LO;
                    end
                end

                %-------------------------------------------------------------
            end

            if i==1
                DA=Area;
            else
                DA=Area-AreaPast;
            end

            Dist=sqrt(dot(Vd,Vd,2));%vector distancia
            pii=floor((Dist./c)./deltat);
            pf=pii+width-posmin;
            pii=pii-width-posmin;
            tiempos=t(pii:pf);
            Coef = c.*tiempos-Dist;
            Pres=-(Coef./Dist).*exp(a.*Coef.^2).*DA.*ArcoF(e).*ILB(d);
            P(pii:pf)=P(pii:pf)+Pres;
            AreaPast=Area;
        end
    end
%     P=P.*ILB(d);
end


end