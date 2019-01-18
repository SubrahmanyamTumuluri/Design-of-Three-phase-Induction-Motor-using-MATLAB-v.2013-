Head=input('Enter the total head(in metres) needed to pump the water');
%Q discharge;
Q=input('Enter the discharge in liters per second\n');
%g=9.81m/s2 andSpecific gravity of water =1000Kg/m3
Waterpower=1000*9.81*Q*(Head*0.001)/1000;
Efficiencyofpump=input('Efficiency of pump');
Efficiencyofmotor=input('Enter the efficiency of motor');
%Work done in KW=WKW
WKW=Waterpower/(Efficiencyofpump*Efficiencyofmotor);
fprintf('Work done in KW is %f\n',WKW);
Hp=(WKW)/0.746;
fprintf('The  nearest rating of the Motor in HP is %f\n',ceil(Hp)); 
HP=ceil(Hp);
Kw=0.746*HP;
phases=3;
Es=415;
f=50;
poles=8;
% winding in Stator
        Esa=Es/1.732; 
 %Synchronous speed ---ns
ns=(2*f)/poles;
fprintf('Speed in RPS is = %d\n',ns);
kw=0.955;
eff=0.81;
 
pf=0.71;
disp('For 50Hz machines  design the value of Bav lies between 0.3 to 0.6 Wb/m2');
Bav=0.4;%Bav in Wb/m2
ac=13250;%ac is ampere conductors per metre.
%Output Coefficient ---Co
Co=11*kw*Bav*ac*10^(-3);
fprintf('Co=%f\n',Co);
%Q--KVA input;
eff*pf;
fprintf(' product eff*pf = %.2f\n',eff*pf);
Q=Kw/(eff*pf);
fprintf('KVA Input Q=%f\n',Q);
p1=Q/(Co*ns);
disp('For minimum cost L/tau varies from 1.5 to 2');
disp('Good power factor L/tau varies from 1 to 1.25');
disp('For Good Efficiency  L/tau ratio is 1.5');
disp('For good overall design L/tau ratio is 1');
ratio=1.5;
D=(p1/(pi*ratio/poles))^0.33;
fprintf('The value of D in metres  is %f\n',D);
L=(pi*ratio/poles)*D;
fprintf('The length in metres  =%f\n',L);
%Net Iron Length ---Li
Li=0.9*(L);
fprintf('Net Iron Length in metres = %f\n',Li);
tau=(pi*D)/poles;
fprintf('tau----Pole Pitch =%f\n',tau);
Fluxperpole =Bav*tau*L;
fprintf('Flux per pole in Wb =%f\n',Fluxperpole);
%Stator Turns per phase ---Ts
Ts=Esa/(4.44*f*kw*Fluxperpole);
fprintf('Ts=%d\n',ceil(Ts));
qs=2;
%(qs>=2)  that is two or more slots per pole per phase is good for design.
% Stator Slots ---Ss
Ss=phases*poles*qs;
coilspan=Ss/poles;
fprintf('Coil Span =%d\n',coilspan);
%Angle of chording---alpha
if(((coilspan/2)-floor(coilspan/2))==0);
    alpha =180/coilspan;fprintf('Alpha in degrees  is : %f',alpha);
    Kp=cosd(alpha/2);
    fprintf('Kp=%f\n',Kp);
else if(((coilspan/2)-floor(coilspan/2))~=0)
    %fprintf('Cs is an odd number hence full pitch winding is used for which Kp=1\n');
    Kp=1;
    end
end
Kd=sind(60/2)/(qs*sind(60/(2*qs)));
fprintf('Kd=%f\n',Kd);
Kwss=Kp*Kd;
fprintf('Kwss=%f\n',Kwss);
%Stator Line Current ---Is
Is=(Kw*1000)/(3*Esa*pf*eff);
fprintf('Stator Line Current in Amperes is %f\n ',Is);
StatorConductors=6*Ts;fprintf('Stator Conductors = %f\n',StatorConductors);
%Stator Slot Pitch Yss
yss = pi*(D*1000)/Ss;
fprintf('Stator slot pitch in mm is %f\n',yss);
if(yss>25);
    fprintf('the stator slot pitch of small induction motor should be less than 25mm.');
   break;
end
currentdensity=4;
AreaofStatorConductor=Is/currentdensity;
fprintf('Area of Stator Conductor  required in mm2=%f\n',AreaofStatorConductor);
Zss=(phases*poles*Ts)/Ss;
statorconductorsperslot=(phases*poles*Ts)/Ss;
fprintf('Stator Conductors per slot=%f\n',ceil(statorconductorsperslot));
Zss*AreaofStatorConductor;
fprintf('Space Required for bare conductor in each slot =%f\n',Zss*AreaofStatorConductor);
%Slot Dimensions------------------------------
spacefactor=0.4;
(Zss*AreaofStatorConductor)/spacefactor;
fprintf('Area of each slot in mm^2 = %f\n',Zss*AreaofStatorConductor/spacefactor);
%Minimum width of stator teeth ---Wts----------------------------
Wts=Fluxperpole*poles*1000/(1.7*Ss*Li);
fprintf('Minimum Width of stator teeth in mm is =%.1f\n',Wts);
 
lip=1;wedge=1.95;height=12.5;
dss=lip+wedge+height;
fprintf('Depth of stator Slot is %0.2f\n',floor(dss));
%Flux Density in Stator Tooth---Fdst-------------------------------------
Fdst=Fluxperpole*poles/(Ss*Li*6*10^(-3));
if(Fdst>=1.7)
    fprintf('Change the Slot Dimensions and start\n');
    break
end
fprintf('Flux Density in Stator Tooth in Wb/m2 =%f\n',Fdst);
%Stator Core
%area of stator core----Acs.Flux  density in stator core is assumed as 1.5wb/m2.
Acs=Fluxperpole/(2*1.5);
dcs=Acs/(Li);
fprintf('Depth of stator Core in mm  =%0.2f\n',dcs*1000);
Bcs=(dcs*1000*1.5)/(ceil(dcs*1000));
fprintf('Value of Stator core flux density in Wb/m2 is %0.2f\n',Bcs);
Do=(D*1000)+(2*dss+2*dcs*1000);
fprintf('Outside diameter including Laminations in mm =%0.2f\n',Do);
lg=0.2+(2*sqrt(D*L));
fprintf(' THE LENGTH OF AIR GAP OF THE INDUCTION MOTOR in mm =%.4f  \n',lg);
Dr = D*1000-(2*lg);
fprintf(' THE OUTER DIAMETER OF ROTOR in mm =%.4f mm \n',Dr);
Lmts=(2*L*1000)+(2.3*tau*1000)+24;
fprintf('Mean length of turns in mm =%0.3f\n',Lmts);
% to find no of rotor slots 
disp('In order to avoid synchronous cusps the difference of stator and rotor slots should not be +1*poles or -1*Poles or a multiple of poles');
Sr=Ss-(9*poles/4);
fprintf('the no of rotor slots is =%d\n',Sr);
Ysr=(3.14*D*1000)/Sr;
fprintf('the rotor slot pitch in mm =%d\n',Ysr);
if(Ysr>25)
fprintf('Modify the value of rotor slots');
end
%to find the rotor bar current
Ib=pf*((phases*qs*Kwss*Is*Ts)/Sr);%Power Factor is 0.71
fprintf('\n the rotor bar current in Amperes=%0.4f\n',Ib);
% to find  rotor bar cross sectional area
Db=6;%Current Density
Ab = Ib /Db;
fprintf('\n THE  AREA OF  each ROTOR BAR in mm2 =%f\n',Ab);
Wsr=6 ;
depth=8.5;
k=Wsr*depth;
if(k<=Ab)
    fprintf('Change the values of the width of rotor bar and depth such that it is slightly greater or equal to area of rotor bar');
    break;
end
%Rotor core
%dcr----Depth of Rotor Core
dcr=dcs;
%Di----Inner Diameter of Rotor Lamination 
Di=Dr-(2*depth)-(2*dcr);
fprintf('Inner Diameter of Rotor Lamination in mm =%f\n',Di);
 %the slot pitch at bottom of rotor bar in mm is  ;
 Ysb=3.14*((1000*D)-(2*depth))/Sr ;
Wt=Ysb-Wsr;
fprintf('\n  the tooth width at root of rotor bar in mm is  =%f\n',Wt);
Fdrt=Fluxperpole/((Sr/poles)*Li*(Wt*0.001));
fprintf('Flux density at root of rotor bar is %f\n',Fdrt);
if(Fdrt>1.7)
    fprintf('Modify the rotor slots and dimensions of rotor bar \n');
end
%copper loss in rotor bars
%A---Allowance for skewing in mm
A=12;
Lb=(L*1000)+A;
fprintf(' the length of rotor bar  in mm is  =%f\n',Lb);
Rb=0.021*Lb*0.001/Ab;
CLRb=Sr*Rb*(Ib^(2));
fprintf('the copper loss in rotor bars in Watts is=%f\n',CLRb);
%End Ring Current
Ie= (Sr*Ib)/(3.14*poles);
fprintf(' the end ring current in Amperes=%f\n',Ie);
%----------------------Area of end ring
De=6; 
Ae = Ie/De;
fprintf(' the area of end rings in mm2 =%f\n',Ae);
%Copper loss in End Rings  
%Diameter of outer end ring and inner ring are Doe and Die.
Doe=Dr-(2*depth);
Die=Doe-(2*10);
Dme=(Doe+Die)/2;
fprintf('Mean Diameter of end ring in mm is %f\n',Dme);
%Resistance of each end ring Re
Re= 0.021*(Dme*0.001*3.14/ Ae );
%Copper loss in two end rings ---Cr
CLER=2*Ie*Ie*Re;
fprintf(' the total copper loss in end rings in watts =%f\n',CLER);
Totalrotorcopperloss=CLRb+CLER;
fprintf('Total Rotot Copper Loss in watts = %f\n',Totalrotorcopperloss);
%--------------------------------------------------------------
slotopening=2;
gaplength=lg;
qi=slotopening/gaplength;
Kcs=1/(1+(5*lg/2));
fprintf('Carter coefficient is %f\n',Kcs);
Kgss=yss/(yss-(Kcs*2));
fprintf('gap contraction factor for stator slots =%0.3f\n',Kgss);
 
%Gap Contraction factor for rotor slots----Kcsr
Kcsr=1/(1+(5*lg/2));
Kgsr=Ysr/(Ysr-Kcsr*2);
fprintf('gap contraction factor for rotor slots=%0.3f\n',Kgsr);
 
Kgs=(Kgss*Kgsr);
fprintf('gap contraction factor for slots=%0.3f\n',Kgs);
%Gap Contraction Factor for ducts---kgd
Kgd=1;
Kg=(Kgs*Kgd);
fprintf('gap contraction factor=%f\n',Kg);
%--------area of air gap---Ag
Ag=pi*D*L/poles;
fprintf('area of air gap in m2 =%f\n',Ag);
Bg=1.36*Bav;
% Length of effective air gap---lge
lge=Kg*lg;
 fprintf('effective length of air gap in mm =%f\n',lge);
ATg=800000*Kg*Bg*0.3*10^(-3);
fprintf('Air gap in Ampere Turns/metre =%f\n',ATg);
%--------stator teeth-------
St=(Ss/poles)*Wts*Li;
fprintf('area of stator teeth per pole in mm2=%f\n',St);
Btssixty=1.36*Fdst;
fprintf('Value of Btssixty is %0.2f\n',floor(Btssixty));
B1=(0.0045-(0.0019*(Btssixty)));
atsst=Btssixty/B1;
fprintf('The ampere turns for stator teeth corresponding to Btssixty are %d\n',atsst);
ATST=atsst*dss*0.001;
fprintf('MMF required for stator teeth are %f\n',ATST);
%------Stator Core---------
%length of magnetic path through stator core ---lcs;
lcs=(pi*(D+(2*dss*0.001)+dcs))/(3*poles);
fprintf('length of stator core in metres =%f\n',lcs);
B2=(0.0045-(0.0019*floor(Bcs)));
ATsc=Bcs/B2;
fprintf('The ampere turns required for Stator core corresp. to Bcs  are %0.2f\n',ATsc);
ATCS=ATsc*lcs;
fprintf('MMF required for stator core =%0.3f\n',ATCS);
%rotor teeth------
Wtsonethird=(pi*(Dr-(4*depth/3))/Sr)-(Wsr);
fprintf('Wtsonethird in mm =%f\n',Wtsonethird);
%Flux density of rotor teeth at one-third height Btr60
Btr60=((1.36*Wt)/Wtsonethird);
fprintf('Flux density in rotor teeth at one third height is %0.2f\n',Btr60);
%Area of rotor teeth at one third height from narrow end.---Atr
Atr=(Sr/poles)*Wtsonethird*0.001*Li;
%fprintf('Atr=%f\n',Atr);
BTrt60=1.36*Btr60;fprintf('BTrt60=%f\n',BTrt60);
B3=(0.0045-(0.0019*(BTrt60)));
ATrt=(BTrt60)/B3;
fprintf('Ampere turns of the rotor teeth corresponding to one third height of rotor teeth is %f\n',ATrt);
ATRT=ATrt*depth*0.001;
fprintf('MMF required for rotor teeth are %0.3f\n',ATRT);
Brc=Bcs;
ATrc=Brc/B2;
fprintf('Ampere turns corresponding to rotor core flux density is %f\n',ATrc);
%lcr----Length of Flux path in Rotor Core
lcr=pi*Di*0.001/(3*poles);
ATRC=ATrc*lcr;
fprintf('MMF for rotor core ATcr=%f\n',ATRC);
AT=ATg+ATST+ATCS+ATRT+ATRC;
fprintf('Total AT = %f\n',AT);
Im=0.427*poles*AT/(Kwss*Ts);
 
fprintf('Im in Amperes =%0.3f\n',Im);
Vsteeth=Ss*dss*0.001*6*Li*0.001;
Wsteeth=Vsteeth*(7.6)*10^(3);
fprintf('Weight of stator teeth in Kg =%f\n',Wsteeth);
Mfdt=MaxFluxDensityinteeth; MaxFluxDensityinteeth=(pi/2)*Fdst;
 %Hystersis loss per cycle is 400J/m3.
ILst=((400*f/7650)+(((pi^(2))*(f^(2))*(Mfdt(2))*(0.0005^(2)))/(6*densitycopper*resistivitycopper)))*Wsteeth;

fprintf('Iron Loss in Stator Teeth = %d\n',ILst);
%Iron Loss in Stator core----Ilsc
Volumestatorcore=3.14*((Do*0.001)-dcs)*dcs*Li;
fprintf('Volume of stator core in kg/m3 = %f\n',Volumestatorcore);
Weightstatorcore=Volumestatorcore*7.65*10^(3);
fprintf('Weight of stator core in kg is %f\n',Weightstatorcore);
ILcore=Weightstatorcore*((400*f/7650)+(((pi^(2))*(f^(2))*(Bcs^(2))*(0.0005^(2)))/(6*densitycopper*resistivitycopper)));
TotalIronLoss=ILst+ILcore;
fprintf('Total Iron Loss in Watts  = %f\n',2*TotalIronLoss);
FandWLoss=4.4*Kw*1000/100;
TotalNoLoadLosses=(2*TotalIronLoss)+FandWLoss;
fprintf('Total No Load Loss=%f\n',TotalNoLoadLosses);
Il=(TotalNoLoadLosses)/(phases*Esa);
fprintf('Loss component of no load current =%f\n',Il);
Io=sqrt((Im*Im)+(Il*Il));
fprintf('No Load Current in Amperes = %d\n',Io);
fprintf('No Load Current as expressed as percentage of Stator Current = %d\n',(Io*100/Is));
%Short Circuit Current
%Stator Slot Leakage
h1=height;
h2=0;
h3=wedge;
h4=lip;
w0=2;w1=(pi*((D*1000)+(2*lip+2*wedge))/Ss)-Wts;w2=11;
%ws---Width at bottom of slot
Ws=14.2;
c=((2*h1)/(3*Ws+3*w2))+(2*h2)/(w1+w2)+(2*h3)/(w1+w0)+(h4/w0);
ssl=4*c*pi*10^(-7);
fprintf('Specific Slot Permeance =%f\n',ssl);
%Rotor Slot Leakage
h5=6;h7=0.3;
h6=1;w01=2;WSr=5.5;
c1=(h5/(3*WSr))+(h7/WSr)+((2*h6)/(WSr+w01))+(h4/w01);
rsp=4*pi*(c1)*10^(-7);
fprintf('Rotor Slot Permeance =%f\n',rsp);
rsps=rsp*(Ss/Sr)*Kwss^2;
Tssp=ssl+rsps;
fprintf('Total Specific Slot Permeance= %f\n',Tssp);
%Slot Leakage Reactance----Xs
Xs=8*pi*f*Ts*Ts*L*Tssp/(poles*qs);
fprintf('Slot Leakage Reactance in ohms =%0.3f\n',Xs);
% Overhang Leakage---X0
X0=8*pi*f*(Ts^(2))*(4*pi*10^(-7)*0.875*tau*tau)/(pi*8*yss*10^(-3));
fprintf('Overhang Leakage Reactance in ohms=%f\n',X0);
%Zigzag Leakage----Xz
Xm=Esa/Im;qr=Sr/poles/phases;
Xz=(5/6)*(Xm/(phases*phases))*((1/(qs*qs))+(1/(qr*qr)));
fprintf('Zigzag Leakage Reactance in ohms =%0.3f\n',Xz);
fprintf('Total Leakage Reactance referred to Stator in ohms = %f\n',(Xs+X0+Xz));
%Resistance of stator winding per phase ----rs
resistivitycopper=0.021;
rs=Lmts*Ts*resistivitycopper/(AreaofStatorConductor*1000);
fprintf('Resistance of stator winding per phase in ohms=%0.2f\n',rs);
Totalstatorcopperloss=phases*Is*Is*rs;
fprintf('Total stator copper loss in watts=%d\n',Totalstatorcopperloss);
rotorcopperlossperphase=Totalrotorcopperloss/phases;
fprintf('Rotor Copper Loss per phase in watts = %0.1f\n',rotorcopperlossperphase);
Rotorresistancereferredtostator=rotorcopperlossperphase/(((Is^(2))*(pf^(2))));
fprintf('Rotor resistance referred to stator in ohms= %0.2f\n',Rotorresistancereferredtostator);
%input----inp
inp=(Kw*1000)+Totalstatorcopperloss+TotalIronLoss+Totalrotorcopperloss+FandWLoss;
fprintf('Input in watts = %d\n',inp);
efficiency=(Kw*1000)/inp;
fprintf('Calculated Efficiency = %0.4f\n',efficiency*100);
slip=1/((Kw*1000/Totalrotorcopperloss)+1);fprintf('Slip is %f\n',slip);
Nr=750*(1-slip);fprintf('Speed of Rotor in RPM is %f\n',Nr);
Zs=sqrt(((rs+Rotorresistancereferredtostator)^(2))+(Xs+X0+Xz)^(2));
Isc=Es/Zs;fprintf('Short Circuit current per phase in Amperes is %0.2f\n',Isc);
fprintf('ratio of short circuit current is %0.2f\n',Isc/Is);
%Total Impadance of rotor at stand still is TIR
TIR=sqrt(((rs+Rotorresistancereferredtostator)^(2))+(Xs+X0+Xz)^(2));
fprintf('Impedance of Rotor  at stand still in ohms is %f\n',TIR);
%No load power factor----Nlpf
Nlpf=(Il/Io);
fprintf('No load power factor is %f\n',Nlpf);
fprintf('No load power factor in degrees is %f\n',acosd(Nlpf));
%short circuit power factor----pfsc
scpf=(Rotorresistancereferredtostator+rs)/sqrt((Rotorresistancereferredtostator+rs)^(2)+(Xs+X0+Xz)^(2));
fprintf('short circuit  power factor is %f\n',scpf);
fprintf('short circuit power factor in degrees is %f\n',acosd(scpf));