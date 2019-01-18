This motor is designed for good power factor and good overall design by taking 
ratio (L/tau) =1.24

%P1= 1 atm
P1=1;
P2=input('Enter the pressure in bar');
V1=input('Enter V1 in m3');
%Gamma for air lies between 1.2 to 1.4
gamma=1.34;
K=((P1/P2)^(1/1.34));
V2=K*V1;
fprintf('V2 in m3/min is %f\n',V2);
%Work done in watts Wwatt
Wwatt=(gamma/(gamma-1))*((P2*V2)-(P1*V1))*((10^(5)))/(60*1000);
fprintf('Work done in kilo watt by compressor is %f\n',Wwatt);
HP=floor(Wwatt)/0.746;
Hp=ceil(HP);
fprintf('Nearest available rating of motor in HP is %f\n',Hp);
Kw=0.746*Hp;
phases=3;
Es=415;
f=50;
poles=4;
% winding in Stator
        Esa=Es 
 %Synchronous speed ---ns
 Ns=(120*f)/poles;
ns=(2*f)/poles;
fprintf('Speed in RPS is = %d\n',ns); kw=0.955; eff=0.925;

pf=0.86;
disp('For 50Hz machines  design the value of Bav lies between 0.3 to 0.6 Wb/m2');
Bav=0.32;%Bav in Wb/m2
ac=16300;%ac is ampere conductors per metre.
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
ratio=1.24;
D=(p1/(pi*ratio/poles))^0.33;
fprintf('The value of D in metres  is %f\n',D);
L=(pi*ratio/poles)*D;
fprintf('The length in metres  =%f\n',L);
%Net Iron Length ---Li
nd=8;Kcd=0.87;wd=8;
Li=0.9*(L-(nd*wd*0.001));
fprintf('Net Iron Length in metres = %f\n',Li);
tau=(pi*D)/poles;
fprintf('tau----Pole Pitch =%f\n',tau);
Fluxperpole =Bav*tau*L;
fprintf('Flux per pole in Wb =%f\n',Fluxperpole);
%Stator Turns per phase ---Ts
Ts=Esa/(4.44*f*kw*Fluxperpole);
fprintf('Ts=%d\n',ceil(Ts));
qs=3;
%(qs<2)is not good for design.
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
is=(22*1000)/(3*Esa*pf*eff);
fprintf('Stator Current per phase in Amperes =%f\n',is);
Is=1.732*is;
fprintf('Stator Line Current in Amperes is %0.3f\n ',Is);
StatorConductors=6*ceil(Ts);fprintf('Stator Conductors = %d\n',StatorConductors);
%Stator Slot Pitch Yss
yss = pi*(D*1000)/Ss;
fprintf('Stator slot pitch in mm is %f\n',yss);
if(yss>25);
    fprintf('the stator slot pitch of induction motor should be less than 25mm.');
    disp('Change the slots per pole per phase');
   break;
end
currentdensity=3.51;
AreaofStatorConductor=is/(currentdensity);
fprintf('Area of Stator Conductor  required in mm2=%0.3f\n',ceil(AreaofStatorConductor));
%Designer specification of area of stator conductor of 7mm2 is 6.25*1.13.
disp('Designer specification of area of stator conductor of 7mm2 is 6.25*1.13');
Zss=(6*Ts)/Ss;
statorconductorsperslot=(6*Ts)/Ss;fprintf('Stator Conductors per slot=%0.1f\n',floor(statorconductorsperslot));
%slot insulation=0.3mm and slack=1mm
slotinsl=0.3;
slack=1;lip=1;wedge=2;
dss=(floor(statorconductorsperslot)*1.13)+lip+(3*slotinsl)+slack+wedge;
fprintf('Total depth of stator slot in mm is = %.1f\n',dss);
WSt=(pi*((D*1000)+(2*(dss)))/Ss)-11;
fprintf('Width at the root of stator teeth  in mm is %f\n',WSt);
%Flux Density at root of stator teeth Fdst
Fdst=Fluxperpole/((Ss/poles)*WSt*Li*0.001);
fprintf('Flux density at root of stator tooth is %f\n',Fdst);
%Stator Core
%area of stator core----Acs.Flux  density in stator core is assumed as 1.5wb/m2.
Acs=Fluxperpole/(2*1.45);
dcs=Acs/(Li);
fprintf('Depth of stator Core in mm  =%0.2f\n',dcs*1000);
Bcs=(dcs*1000*1.5)/(ceil(dcs*1000));
fprintf('Value of Stator core flux density in Wb/m2 is %0.2f\n',Bcs);
Do=(D*1000)+(2*dss)+(2*dcs*1000);
fprintf('Outside diameter including Laminations in mm =%0.2f\n',Do);
Lmts=(2*L*1000)+(2.3*tau*1000)+24;
fprintf('Mean length of turns in mm =%0.3f\n',Lmts);
%Calculation of air gap length and Rotor Design
lg=0.2+D;
fprintf(' THE LENGTH OF AIR GAP OF THE INDUCTION MOTOR in mm =%.4f  \n',lg);
Dr = D*1000-(2*(1.0*lg));
fprintf(' THE OUTER DIAMETER OF ROTOR in mm =%.4f mm \n',Dr);
% to find no of rotor slots 
disp('In order to avoid synchronous cusps the difference of stator and rotor slots should not be +1*poles or -1*Poles or a multiple of poles');
 
 
Sr=Ss+(9*poles/4);
fprintf('the no of rotor slots is =%d\n',Sr);
Ysr=(3.14*D*1000)/Sr;
fprintf('the rotor slot pitch in mm =%0.4f\n',Ysr);
if(Ysr>25)
fprintf('Modify the value of rotor slots');
break;
end
%to find the rotor bar current
Ib=0.86*((6*Is*Kwss*Ts)/Sr);
fprintf('\n the rotor bar current in Amperes  =%.4f\n',Ib);
% to find  rotor bar cross sectional area
Db=4.5;%Current Density
Ab = Ib /Db;
fprintf('\n THE  AREA OF  each ROTOR BAR in mm2 =%f \n',ceil(Ab));
%designer specification to use Deep T rotor bars of ((10*7.9)+(10*1))
disp('designer specification to use Deep T rotor bars of ((10.4*6.7)+(10*2.2))');
depth=19.4;
width=7.1;
Ysb=pi*((D*1000)-(2*depth))/Sr;
Wtr=(pi*((D*1000)-(2*depth))/Sr)-7.1;
fprintf('Width at root of rotor teeth in mm is %f\n',Wtr);
%Flux Density at root of rotor teeth Fdrt
Fdrt=Fluxperpole/((Ss/poles)*Li*Wtr*0.001);
fprintf('Flux density at root of rotor tooth is %f\n',Fdrt);
if(Fdrt>1.7);
    fprintf('Verify the value of width at rotor of rotor teeth and calculate Fdrt');
    break;
end
 
dcr=dcs;
Di=Dr-(2*depth)-(2*dcr);fprintf('Inner Diameter of Rotor Lamination in mm =%f\n',Di);
 %the slot pitch at bottom of rotor bar in mm is  ;
 
 
%copper loss in rotor bars
%A---Allowance for skewing in mm
A=10;
Lb=(L*1000)+A;
fprintf(' the length of rotor bar in mm is  =%f\n',Lb);
Rb=0.021*Lb*0.001/Ab;
CLRb=Sr*Rb*(Ib^(2));
fprintf('the copper loss in rotor bars in Watts is=%f\n',CLRb);
%End Ring Current
Ie= (Sr*Ib)/(3.14*poles);
fprintf(' the end ring current in Amperes =%f\n',Ie);
De=4%End ring current density is 4 A/mm2; 
Ae = Ie/De;
fprintf(' the area of end rings in mm2=%f\n',Ae);
%Copper loss in End Rings  
%Diameter of outer end ring and inner ring are Doe and Die.
Doe=Dr-(2*depth);
Die=Doe-(2*50);
Dme=(Doe+Die)/2;
fprintf('Mean Diameter of end ring in mm is %f\n',Dme);
%Resistance of each end ring Re
Re= 0.021*(Dme*0.001*3.14/ Ae );
%Copper loss in two end rings ---Cr
%Since two ed rings are present,Copper loss n two end rings are 2times.
CLER=2*Ie*Ie*Re;
fprintf(' the total copper loss in end rings in watts =%0.3f\n',CLER);
Totalrotorcopperloss=CLRb+CLER;
fprintf('Total Rotor[ Copper Loss in watts = %0.3f\n',Totalrotorcopperloss);
%--------------------------------------------------------------
slotopening=5;
gaplength=lg;
qi=slotopening/gaplength;
Kcs=1/(1+(5*lg/slotopening));
fprintf('Carter coefficient is %f\n',Kcs);
Kgss=yss/(yss-(Kcs*5));
fprintf('gap contraction factor for stator slots =%0.3f\n',Kgss);
 
%Gap Contraction factor for rotor slots----Kcsr
Kcsr=1/(1+(5*lg/2.2));
Kgsr=Ysr/(Ysr-Kcsr*2.2);
fprintf('gap contraction factor for rotor slots=%0.3f\n',Kgsr);
 
Kgs=(Kgss*Kgsr);
fprintf('gap contraction factor for slots=%0.3f\n',Kgs);
%Gap Contraction Factor for Ducts
nd=8;Kcd=0.87;wd=8;
Kgd=(L*1000)/((L*1000)-(nd*wd*Kcd));
fprintf('Gap Contraction Factor for ducts is %f\n',Kgd);
Kg=Kgs*Kgd;
lge=Kg*lg;Bgs=1.36*Bav;
 fprintf('effective length of air gap in mm =%f\n',lge);
ATg=800000*Kg*Bgs*lg*10^(-3);
fprintf('MMF of Air Gap is %f\n',ATg);
St=(Ss/poles)*WSt*Li;
fprintf('area of teeth per pole in mm2=%f\n',St);
Btssixty=1.36*Fdst;
fprintf('Value of Btssixty is %f\',floor(Btssixty));
B1=(0.0045-(0.0019*(Btssixty)));
atsst=Btssixty/B1;
fprintf('The ampere turns for stator core corresponding to Btssixty are %d\n',atsst);
ATST=atsst*dss*0.001;
fprintf('MMF required for stator teeth are %f\n',ATST);
%------Stator Core---------
%length of magnetic path through stator core ---lcs;
lcs=(pi*(D+(2*dss*0.001)+dcs))/(3*poles);
fprintf('length of stator core in metres =%f\n',lcs);
B2=(0.0045-(0.0019*floor(Bcs)));
ATsc=Bcs/B2;
fprintf('The ampere conductors required for Stator core corresp. to Bcs  are %0.2f\n',ATsc);
ATCS=ATsc*lcs;
fprintf('MMF required for stator core =%0.3f\n',ATCS);
%rotor teeth------
Wtsonethird=(pi*(Dr-(4*depth/3))/Sr)-(Wtr);
fprintf('Wtsonethird in mm =%f\n',Wtsonethird);
%Flux density of rotor teeth at one-third height Btr60
Btr60=((1.36*Wtr)/Wtsonethird);
fprintf('Flux density in rotor teeth at one third height is %0.2f\n',Btr60);
Atr=(Sr/poles)*Wtsonethird*0.001*Li;%Area of rotor teeth at one third height from narrow end.---Atr
 
BTrt60=1.36*Btr60;fprintf('BTrt60=%f\n',BTrt60);
B3=(0.0045-(0.0019*(BTrt60)));
ATrt=(BTrt60)/B3;
fprintf('Ampere turns of the rotor teeth corresponding to Btr60 is %f\n',ATrt);
ATRT=ATrt*depth*0.001;
fprintf('MMF required for rotor teeth are %0.3f\n',ATRT);
Brc=Bcs;
%atcr----Ampere Turns corresponding to Rotor core
ATrc=Brc/B2;
fprintf('Ampere turns corresponding to rotor core flux density is %f\n',ATrc);
%lcr----Length of Flux path in Rotor Core
lcr=pi*Di*0.001/(3*poles);
ATRC=ATrc*lcr;
fprintf('MMF for rotor core ATcr=%f\n',ATRC);
AT=ATg+ATST+ATCS+ATRT+ATRC;
fprintf('Total AT = %f\n',AT);
Im=0.427*poles*AT/(Kwss*Ts);
fprintf('Magnetizing Current in Amperes is %f\n',Im);
Vsteeth=Ss*dss*0.001*WSt*Li*0.001;
Wsteeth=Vsteeth*7.6*10^(3);
fprintf('Weight of stator teeth in Kg =%f\n',Wsteeth);
Mfdt=MaxFluxDensityinteeth;MaxFluxDensityinteeth=(pi/2)*Fdst;
ILst=((400*f/7650)+(((pi^(2))*(f^(2))*(Mfdt^(2))*(0.0005^(2)))/(6*densitycopper*resistivitycopper)))*Wsteeth;
fprintf('Iron Loss in Stator Teeth = %d\n',ILst);
Volumestatorcore=3.14*((Do*0.001)-dcs)*dcs*Li;fprintf('Volume of stator core in Kg/m3 = %f\n',Volumestatorcore);
Weightstatorcore=Volumestatorcore*7.65*10^(3);
ILcore=Weightstatorcore*((400*f/7650)+(((pi^(2))*(f^(2))*(Bcs^(2))*(0.0005^(2)))/(6*densitycopper*resistivitycopper)));
TotalIronLoss=ILst+ILcore;
fprintf('Total Iron Loss in Watts  = %d\n',2*TotalIronLoss);
FandWLoss=1.8*Kw*1000/100;
Totalnoloadloss=(2*TotalIronLoss)+FandWLoss
fprintf('Total No Load Loss in Watts =%d\n',Totalnoloadloss);
%Loss component of No Load current ---Il
Il=(Totalnoloadloss)/(phases*Esa);
fprintf('Loss Component of No Load Current in Amperes=%f\n',Il);
Io=sqrt((Im*Im)+(Il*Il));fprintf('No Load Current in Amperes = %f\n',Io);
fprintf('No Load Current as expressed as percentage of Stator Current = %d\n',(Io*100/Is));
%Short Circuit Current
%Stator Slot Leakage
h1=Zss*1.13;
h2=1;   h3=2;  h4=1; ws1=10.3; w0=5;
c0=4*pi*(10^(-7))*((h1/(3*ws1))+(h2/ws1)+(2*h3/(ws1+w0))+(h4/w0));
fprintf('Specific Slot Permeance of Stator =%f\n',c0);
a1=10.4*7.1;a2=10*2.2;
H1=10.4;H2=9;H3=1;Ws=7.1;W0=2.2;
c1=4*pi*(10^(-7))*(((H1*a1*a1/(3*Ws))+(H2/w0)*((a1*a1)+(a1*a2)+(a2*a2)/3))/((a1+a2)^(2))+(H3/Ws));
fprintf('Specific Slot Permeance of Rotor is %f\n',c1);
c2=c1*Kwss^(2)*Ss/Sr
SlotLeakageReactance=8*pi*f*(Ts^(2))*L*(c0+c2)/(poles*qs);
fprintf('Slot leakage Reactance is %f\n',SlotLeakageReactance);
% Overhang Leakage---X0
X0=8*pi*f*(Ts^(2))*(4*pi*10^(-7)*tau*tau)/(poles*qs*pi*yss*10^(-3));
fprintf('Overhang Leakage Reactance in ohms=%f\n',X0);
%Zigzag Leakage----Xz
Xm=Esa/Im;qr=Sr/poles/phases;
Xz=(5/6)*(Xm/(phases*phases))*((1/(qs*qs))+(1/(qr*qr)));
fprintf('Zigzag Leakage Reactance in ohms =%0.3f\n',Xz);
Xs=SlotLeakageReactance;
fprintf('Total Leakage Reactance referred to Stator in ohms = %f',(Xs+X0+Xz))
%Resistance of stator winding per phase ----rs
resistivitycopper=0.021;
rs=Lmts*Ts*resistivitycopper/(AreaofStatorConductor*1000);
fprintf('Resistance of stator winding per phase in ohms=%0.2f\n',rs);
Totalstatorcopperloss=phases*Is*Is*rs;fprintf('Total stator copper loss in watts=%d\n',Totalstatorcopperloss);
rotorcopperlossperphase=Totalrotorcopperloss/phases;
fprintf('Rotor Copper Loss per phase in watts = %0.1f\n',rotorcopperlossperphase);
Rotorresistancereferredtostator=rotorcopperlossperphase/((is^(2))*(pf^(2)));
fprintf('Rotor resistance referred to stator in ohms= %0.2f\n',Rotorresistancereferredtostator);
 
%input----inp
inp=(Kw*1000)+Totalstatorcopperloss+TotalIronLoss+Totalrotorcopperloss+FandWLoss;
fprintf('Input in watts = %d\n',inp);
efficiency=(Kw*1000)/inp;
fprintf('Efficiency is= %f\n',efficiency*100);
slip=1/((Kw*1000/Totalrotorcopperloss)+1);fprintf('Slip is %f\n',slip);
Nr=Ns*(1-slip);fprintf('Speed of Rotor in RPM is %f\n',Nr);
 
%Total Impadance of rotor at stand still is TIR
TIR=sqrt(((rs+Rotorresistancereferredtostator)^(2))+(Xs+X0+Xz)^(2));fprintf('Rotor Resistance at stand still in ohms is %f\n',TIR);
%No Load power factor calculation---Nlpf
Nlpf=(Il/Io);
fprintf('No load power factor is %f\n',acosd(Nlpf));
%Short circuit power factor ----scpf
scpf=((rs+Rotorresistancereferredtostator)/sqrt(((rs+Rotorresistancereferredtostator)^(2))+(Xs+X0+Xz)^(2)));
fprintf('Short circuit power factor %f\n',acosd(scpf));
Isc=Es/TIR;fprintf('Short Circuit current per phase  in Amperes is %0.2f\n',Isc);
fprintf('ratio of short circuit current to stator per phase current is %0.2f\n',Isc/is);


