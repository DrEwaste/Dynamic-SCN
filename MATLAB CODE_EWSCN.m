function main_ahppcaS1
%PROGRAM FOR AHP+PCA COMBINED METHOD FOR MECHANICAL RECYCLING PLANT SUPPLY CHAIN
xa = 0; xb = 900;
solinit = bvpinit(linspace(xa,xb,3),@ew_init);
options=bvpset('Reltol',1e-3, 'AbsTol',1e-5);
sol = bvp4c(@ew_func,@bc,solinit,options);

xint = linspace(xa,xb,1000);
S_i = deval(sol,xint);

figure(1)
plot(xint,S_i(1,:),'k-.');
grid on
xlabel('time (day)');
ylabel('Volume of CO2 generated (VCO2)_ AHP');
saveas(gcf,'C_VCO2_AHP','jpg')
saveas(gcf,'C_VCO2_AHP','eps')

figure(2)
plot(xint,S_i(2,:),'k-.');
grid on
xlabel('t(day)');
ylabel('Energy Consumed(Ec)_ AHP');
saveas(gcf,'C_Ec_AHP','jpg')
saveas(gcf,'C_Ec_AHP','eps')

figure(3)
plot(xint,S_i(3,:),'k-.');
grid on
xlabel('t(day)');
ylabel('No. of awareness activities (N3)_ AHP');
saveas(gcf,'C_N3_AHP','jpg')
saveas(gcf,'C_N3_AHP','eps')

figure(4)
plot(xint,S_i(4,:),'k-.');
grid on
xlabel('time(day)');
ylabel('No. of product sold (N4)_ AHP');
saveas(gcf,'C_N4_AHP','jpg')
saveas(gcf,'C_N4_AHP','eps')

xpa = 0; xpb = 900;
solinit = bvpinit(linspace(xpa,xpb,3),@ewp_init);
options=bvpset('Reltol',1e-3, 'AbsTol',1e-5);
sol = bvp4c(@ew_funcp,@bc,solinit,options);

xpint = linspace(xpa,xpb,1000);
Sp_i = deval(sol,xpint);

figure(5)
plot(xpint,Sp_i(1,:),'k-.');
grid on
xlabel('time (day)');
ylabel('Volume of CO2 generated (VCO2)_ AHP+PCA');
saveas(gcf,'C_VCO2_PCA','jpg')
saveas(gcf,'C_VCO2_PCA','eps')

figure(6)
plot(xpint,Sp_i(2,:),'k-.');
grid on
xlabel('t(day)');
ylabel('Energy Consumed(Ec)_ AHP+PCA');
saveas(gcf,'C_Ec_PCA','jpg')
saveas(gcf,'C_Ec_PCA','eps')

figure(7)
plot(xpint,Sp_i(3,:),'k-.');
grid on
xlabel('t(day)');
ylabel('No. of awareness activities (N3)_ AHP+PCA');
saveas(gcf,'C_N3_PCA','jpg')
saveas(gcf,'C_N3_PCA','eps')

figure(8)
plot(xpint,Sp_i(4,:),'k-.');
grid on
xlabel('time(day)');
ylabel('No. of product sold (N4)_ AHP+PCA');
saveas(gcf,'C_N4_PCA','jpg')
saveas(gcf,'C_N4_PCA','eps')

figure(9)
plot(xint,S_i(1,:),'k-.');
grid on
hold on
plot(xpint,Sp_i(1,:),'k-');
hold on 
grid on
xlabel('time(day)');
ylabel('Volume of CO2 generated (VCO2)_ Comp');
hold on
legend('VCO2_ AHP','VCO2_ AHP+PCA','Location','best')
saveas(gcf,'Comp_VCO2','jpg')
saveas(gcf,'Comp_VCO2','eps')

figure(10)
plot(xint,S_i(2,:),'k-.');
grid on
hold on
plot(xpint,Sp_i(2,:),'k-');
hold on 
grid on
xlabel('time (day)');
ylabel('Energy Consumed (Ec)_ Comp');
hold on
legend('Ec_ AHP','EC_ AHP+PCA','Location','best')
saveas(gcf,'Comp_Ec','jpg')
saveas(gcf,'Comp_Ec','eps')

figure(11)
plot(xint,S_i(3,:),'k-.');
grid on
hold on
plot(xpint,Sp_i(3,:),'k-');
hold on 
grid on
xlabel('time (day)');
ylabel('No. of Awareness (N3)_ Comp');
hold on
legend('N3_ AHP','N3_ AHP+PCA','Location','best')
saveas(gcf,'Comp_N3','jpg')
saveas(gcf,'Comp_N3','eps')

figure(12)
plot(xint,S_i(4,:),'k-.');
grid on
hold on
plot(xpint,Sp_i(4,:),'k-');
hold on 
grid on
xlabel('time (day)');
ylabel('No. of product sold (N4)_ Comp');
hold on
legend('N4_ AHP','N4_ AHP+PCA','Location','best')
saveas(gcf,'Comp_N4','jpg')
saveas(gcf,'Comp_n4','eps')

function res = bc(ya,yb)
res = [ ya(1)-1.2; ya(2)-0.002; ya(3)-2; ya(4)-0.3; yb(1)-0.84;yb(2)-0.0015;yb(3)-4;yb(4)-0.6];

function dxdt=ew_func(t,x)
%Defining Parameters and Variables of the Model

Vco2=1.2;           %Volume of CO2 generated
Ec=0.002;           %Energy Consumed (GW)
Wp=30000;           %Water Consumed (Litres)
Ww=Wp*0.9;          %Waste Water Produced
N1=16;              %No of labours
N3=2;               %number of awareness activities and markenting per year 
N4=6;               %number of recycled materials 
N5=7;               %number of operations
N7=2;               %number of logistics
N8=2;               %number of waste materials going to TSDF
N9=1;               %number of taxes
f1=10000;           %unit cost of CO2 recovery 
f2=20;              %unit cost of  energy used 
f3=1500;            %unit cost of water used 
f4=60000;           %unit cost of wastewater treatment 
f5=108000;          %salary of one labour  
f6=100000;          %average cost of awareness activity & makreting
f7=4260000;         %unit revenue earned from recycled product sold 
f8=4500;            %unit cost of each operation 
f10=180000;         %unit cost of disposal in TSDF
f11=100000;         %unit cost of Tax

%Defining Weigh Factors derived from AHP

e1=0.41606;
e2=0.45793;
e3=0.12601;
A1=0.1793;
A2=0.0631;
A3=0.0240;
A4=0.1497;
A5=0.0315;
A6=0.0945;
A7=0.1998;
A8=0.1193;
A9=0.0850;
A10=0.0319;
Ah11=0.0219;

%Defining Inter-dependency values derived from AHP

Z=0;
a1=0.25;
a2=0.75;
a12=a1*a2;
ap1=(a1)^0.5;
ap2=(a2)^0.5;
b1=0.16667;
b2=0.83333;
b12=b1*b2;
bp1=(b1)^0.5;
bp2=(b2)^0.5;
c1=0.75;
c2=0.25;
d1=0.33333;
d2=0.66667;
d12=d1*d2;
dp1=(d1)^0.5;
dp2=(d2)^0.5;
alpha1=0.83333;
alpha2=0.16667;
alpha12=alpha1*alpha2;
alphap1=(alpha1)^0.5;
alphap2=(alpha2)^0.5;
beta1=0.5;
beta2=0.5;
beta12=beta1*beta2;
betap1=(beta1)^0.5;
betap2=(beta2)^0.5;
gamma=1;

%defining the scaling factor

q=(1/30000);

%Defining Inter-dependency values of Lagrange multipliers

l1=0.41606;
l2=0.45793;
l3=0.12601;
%l4=0.1497;

%Defining Kappa (ki)values for Unconstrained System in case of a Mechanical Recyling plant

%k1=e1*A1*f1;
%k2=e1*A2*f2;
%k3=e2*A6*f6;
%k4=e3*A7*f7;

k1=(e1*A1*f1-l1*f1);
k2=e1*A2*f2;
k3=(e2*A6-l2)*f6;
k4=(e3*A7-l3)*f7;

%defining the scaling factor

%q=(1/30000);

%Defining Elements of the Matix 'A' i.e. for 
%Unconstrained System for a Mechanical Recyling plant

Ak1=((2*ap1*e3*A8*f8)/((a1+a12*N7+2*ap1*N5)^2))-((e1*A1*f1*a12)/(gamma*((a2+a12*N7+2*ap1*N5)^2)));
Ak2= 2*ap1*e3*A8*f8*((a1+a12*N7+2*ap1*N5)^3);
Ak3=((-2)*ap2*e1*A1*f1*(1/gamma))*((a2+a12*N5+2*ap2*N7)^3);
Ak4=((2*ap2)/(a2+a12*N7+2*ap1*N5))+((a12)/(a1+a12*N7+2*ap1*N5));
Ak5=((1/(a1+a12*N7+2*ap1*N5))*Ak1)+Ak2+Ak3+(e1*A1*f1*gamma*Ak4);
A11=Ak5*q;

At1=e3/(a1+a12*N7+2*ap2*N7);
At2=((2*bp2*A8*f8)/((b2+b12*N4+2*bp2*N5)^2));
At3=((A7*f7*b12)/((b1+b12*N5+2*bp1*N4)^2));
At4=(gamma*e1*A1*f1*a12)/(b2+b12*N4+2*bp2*N5);
A12=((At1*(At2-At3))+At4)*q;
A13=0;
A31=0;
Ae1=(1/(a1+a12*N7+2*ap1*N5));
Ae2=(e1*A2*f2*b12);
Ae3=(2*e3*A8*f8*betap2)/((beta2+beta12*Ec+2*betap2*N5)^2);
A14=(Ae1*(Ae2+Ae3))*q;
A21=A12;
Aq1=e3/(b1+b12*N4+2*bp2*N5);
Aq2=((2*bp2*A8*f8)/((b2+b12*N4+2*bp2*N5)^2))-((A7*f7*b12)/((b1+b12*N5+2*bp1*N4)^2));
Aq3=2*bp2*e3*A8*f8*((b2+b12*N4+2*bp2*N5)^3);
Aq4=(beta1+beta12*N5+2*betap1*Ec);
Aq5=((e3*b12*A8*f8)/((b2+b12*N4+2*bp2*N5)^2))-((e3*A7*f7*2*bp1)/((b1+b12*N5+2*bp1*N4)^2));
Aq6=((2*bp1)/((beta1+beta12*N5+2*betap1*Ec)^3));
Aq7=((e1*A2*f2)*(b1+b12*N5+2*bp1*N4))-((e3*A8*f8)/((beta2+beta12*Ec+2*betap2*N5)^2));
A22=((Aq1*Aq2)+Aq3+(Aq4*Aq5)-(Aq6*Aq7))*q;
As1=(beta1+beta12*N5+2*betap1*Ec);
As2=((alpha12*Ah11*f11)/((alpha2+alpha12*N4+2*alphap2*N9)^2))-((2*alphap1*A7*f7)/((alpha2+alpha12*N9+2*alphap1*N4)^2));
A23=(As1*(e3*As2))*q;
A24=(e3*(((b12*A8*f8)/((b2+b12*N4+2*bp2*N5)^2))-((A7*f7*2*bp1)/((b1+b12*N5+2*bp1*N4)^2))))*q;
A32=A23;
Aw1=(e3/(alpha1+alpha12*N9+2*alphap1*N4));
Aw2=((alpha12*Ah11*f11)/((alpha2+alpha12*N4+2*alphap2*N9)^2))-((2*alphap1*A7*f7)/((alpha2+alpha12*N9+2*alphap1*N4)^2));
Aw3=(-2*alphap1)*((alpha1+alpha12*N9+2*alphap1*N4)^3);
Aw4=((e1*A2*f2)*(b1+b12*N5+2*bp1*N4))-((e3*A8*f8)/((beta2+beta12*Ec+2*betap2*N5)^2));
Aw5=(2*alphap2*e3*A11*f11)*((alpha2+alpha12*N4+2*alphap2*N9)^3);
A33=((Aw1*Aw2)+(Aw3*Aw4)+Aw5)*q;
A34=(e3*(((alpha12*Ah11*f11)/((alpha2+alpha12*N4+2*alphap2*N9)^2))-((A7*f7*2*alphap1)/((alpha1+alpha12*N9+2*alphap1*N4)^2))))*q;
A41=A14;
A42=A24;
A43=A34;
Ar1=e3*(b1+b12*N4+2*bp2*N5);
Ar2=((b12*A8*f8)/((b2+b12*N4+2*bp2*N5)^2))-((2*bp1*A7*f7*b12)/((b1+b12*N5+2*bp1*N4)^2));
Ar3=((-2*e3*betap1)/((b1+b12*N4+2*bp2*N5)^3));
Ar4=((A7*f7)/((b1+b12*N4+2*bp2*N5)^2))-((A8*f8)/((b2+b12*N4+2*bp2*N5)^2));
Ar5=(1/(beta2+beta12*Ec+2*betap2*N5));
Ar6=(e1*b12*A2*f2)+((2*betap2*e3*A8*f8)/((beta2+beta12*Ec+2*betap2*N5)^2));
Ar7=(2*betap2*e3*A8*f8)*((beta2+beta12*Ec+2*betap2*N5)^3);
A44=((Ar1*Ar2)+(Ar3*Ar4)+(Ar5*Ar6)+Ar7)*q;

%Defining Elements of the Matix 'C' i.e. for 
%Unconstrained System for a Mechanical Recyling plant

C11=A11/k1;
C12=A12/k1;
C13=A13/k1;
C14=A14/k1;
C21=A21/k2;
C22=A22/k2;
C23=A23/k2;
C24=A24/k2;
C31=A31/k3;
C32=A32/k3;
C33=A33/k3;
C34=A34/k3;
C41=A41/k4;
C42=A42/k4;
C43=A43/k4;
C44=A44/k4;

dxdt=[x(5);
      x(6);
      x(7);
      x(8);
     (C11)*x(1)+(C12)*x(2)+(C13)*x(3)+(C14)*x(4);
     (C21)*x(1)+(C22)*x(2)+(C23)*x(3)+(C24)*x(4);
     (C32)*x(2)+(C33)*x(3)+(C34)*x(4);
     (C41)*x(1)+(C42)*x(2)+(C33)*x(3)+(C44)*x(4)];

function v=ew_init(x)
        v=[1.2;
           0.002;
           2;
           0.3;
           0;
           0;
           0;
           0]; 
       
       

function dxpdt=ew_funcp(t,xp)
%Defining Parameters and Variables of the Model

Vco2=1.2;           %Volume of CO2 generated
Ec=0.002;           %Energy Consumed (GW)
Wp=30000;           %Water Consumed (Litres)
Ww=Wp*0.9;          %Waste Water Produced
N1=16;              %No of labours
N3=2;               %number of awareness activities and markenting per year 
N4=6;               %number of recycled materials 
N5=7;               %number of operations
N7=2;               %number of logistics
N8=2;               %number of waste materials going to TSDF
N9=1;               %number of taxes
f1=10000;           %unit cost of CO2 recovery 
f2=20;              %unit cost of  energy used 
f3=1500;            %unit cost of water used 
f4=60000;           %unit cost of wastewater treatment 
f5=108000;          %salary of one labour  
f6=100000;          %average cost of awareness activity & makreting
f7=4260000;         %unit revenue earned from recycled product sold 
f8=4500;            %unit cost of each operation 
f10=180000;         %unit cost of disposal in TSDF
f11=100000;         %unit cost of Tax

%Defining Weigh Factors derived from AHP+PCA

e1=0.713;
e2=0.286;
e3=0.001;
A1=0.1142;
A2=0.0083;
A3=0.0636;
A4=0.0954;
A5=0.0200;
A6=0.0602;
A7=0.2521;
A8=0.1135;
A9=0.1047;
A10=0.0726;
Ah11=0.0950;

%Defining Inter-dependency values derived from AHP+PCA

Z=0;
a1=0.17748;
a2=0.25;
a12=a1*a2;
ap1=(a1)^0.5;
ap2=(a2)^0.5;
b1=0.08475;
b2=0.8334;
b12=b1*b2;
bp1=(b1)^0.5;
bp2=(b2)^0.5;
c1=0.6649;
c2=0.4421;
d1=0.0333;
d2=0.667;
d12=d1*d2;
dp1=(d1)^0.5;
dp2=(d2)^0.5;
alpha1=0.834;
alpha2=0.166;
alpha12=alpha1*alpha2;
alphap1=(alpha1)^0.5;
alphap2=(alpha2)^0.5;
beta1=0.3435;
beta2=0.5;
beta12=beta1*beta2;
betap1=(beta1)^0.5;
betap2=(beta2)^0.5;
gamma=1;

%defining the scaling factor

q=(1/30000);

%Defining Inter-dependency values of Lagrange multipliersderived from AHP+PCA

l1=0.6221;
l2=0.2963;
l3=0.0816;

%Defining Kappa (ki)values for Unconstrained System in case of a Mechanical Recyling plant

%k1=e1*A1*f1;
%k2=e1*A2*f2;
%k3=e2*A6*f6;
%k4=e3*A7*f7;

k1=(e1*A1*f1-l1*f1);
k2=e1*A2*f2;
k3=(e2*A6-l2)*f6;
k4=(e3*A7-l3)*f7;

%defining the scaling factor

%q=(1/30000);

%Defining Elements of the Matix 'A' i.e. for 
%Unconstrained System for a Mechanical Recyling plant

Ak1=((2*ap1*e3*A8*f8)/((a1+a12*N7+2*ap1*N5)^2))-((e1*A1*f1*a12)/(gamma*((a2+a12*N7+2*ap1*N5)^2)));
Ak2= 2*ap1*e3*A8*f8*((a1+a12*N7+2*ap1*N5)^3);
Ak3=((-2)*ap2*e1*A1*f1*(1/gamma))*((a2+a12*N5+2*ap2*N7)^3);
Ak4=((2*ap2)/(a2+a12*N7+2*ap1*N5))+((a12)/(a1+a12*N7+2*ap1*N5));
Ak5=((1/(a1+a12*N7+2*ap1*N5))*Ak1)+Ak2+Ak3+(e1*A1*f1*gamma*Ak4);
A11=Ak5*q;

At1=e3/(a1+a12*N7+2*ap2*N7);
At2=((2*bp2*A8*f8)/((b2+b12*N4+2*bp2*N5)^2));
At3=((A7*f7*b12)/((b1+b12*N5+2*bp1*N4)^2));
At4=(gamma*e1*A1*f1*a12)/(b2+b12*N4+2*bp2*N5);
A12=((At1*(At2-At3))+At4)*q;
A13=0;
A31=0;
Ae1=(1/(a1+a12*N7+2*ap1*N5));
Ae2=(e1*A2*f2*b12);
Ae3=(2*e3*A8*f8*betap2)/((beta2+beta12*Ec+2*betap2*N5)^2);
A14=(Ae1*(Ae2+Ae3))*q;
A21=A12;
Aq1=e3/(b1+b12*N4+2*bp2*N5);
Aq2=((2*bp2*A8*f8)/((b2+b12*N4+2*bp2*N5)^2))-((A7*f7*b12)/((b1+b12*N5+2*bp1*N4)^2));
Aq3=2*bp2*e3*A8*f8*((b2+b12*N4+2*bp2*N5)^3);
Aq4=(beta1+beta12*N5+2*betap1*Ec);
Aq5=((e3*b12*A8*f8)/((b2+b12*N4+2*bp2*N5)^2))-((e3*A7*f7*2*bp1)/((b1+b12*N5+2*bp1*N4)^2));
Aq6=((2*bp1)/((beta1+beta12*N5+2*betap1*Ec)^3));
Aq7=((e1*A2*f2)*(b1+b12*N5+2*bp1*N4))-((e3*A8*f8)/((beta2+beta12*Ec+2*betap2*N5)^2));
A22=((Aq1*Aq2)+Aq3+(Aq4*Aq5)-(Aq6*Aq7))*q;
As1=(beta1+beta12*N5+2*betap1*Ec);
As2=((alpha12*Ah11*f11)/((alpha2+alpha12*N4+2*alphap2*N9)^2))-((2*alphap1*A7*f7)/((alpha2+alpha12*N9+2*alphap1*N4)^2));
A23=(As1*(e3*As2))*q;
A24=(e3*(((b12*A8*f8)/((b2+b12*N4+2*bp2*N5)^2))-((A7*f7*2*bp1)/((b1+b12*N5+2*bp1*N4)^2))))*q;
A32=A23;
Aw1=(e3/(alpha1+alpha12*N9+2*alphap1*N4));
Aw2=((alpha12*Ah11*f11)/((alpha2+alpha12*N4+2*alphap2*N9)^2))-((2*alphap1*A7*f7)/((alpha2+alpha12*N9+2*alphap1*N4)^2));
Aw3=(-2*alphap1)*((alpha1+alpha12*N9+2*alphap1*N4)^3);
Aw4=((e1*A2*f2)*(b1+b12*N5+2*bp1*N4))-((e3*A8*f8)/((beta2+beta12*Ec+2*betap2*N5)^2));
Aw5=(2*alphap2*e3*A11*f11)*((alpha2+alpha12*N4+2*alphap2*N9)^3);
A33=((Aw1*Aw2)+(Aw3*Aw4)+Aw5)*q;
A34=(e3*(((alpha12*Ah11*f11)/((alpha2+alpha12*N4+2*alphap2*N9)^2))-((A7*f7*2*alphap1)/((alpha1+alpha12*N9+2*alphap1*N4)^2))))*q;
A41=A14;
A42=A24;
A43=A34;
Ar1=e3*(b1+b12*N4+2*bp2*N5);
Ar2=((b12*A8*f8)/((b2+b12*N4+2*bp2*N5)^2))-((2*bp1*A7*f7*b12)/((b1+b12*N5+2*bp1*N4)^2));
Ar3=((-2*e3*betap1)/((b1+b12*N4+2*bp2*N5)^3));
Ar4=((A7*f7)/((b1+b12*N4+2*bp2*N5)^2))-((A8*f8)/((b2+b12*N4+2*bp2*N5)^2));
Ar5=(1/(beta2+beta12*Ec+2*betap2*N5));
Ar6=(e1*b12*A2*f2)+((2*betap2*e3*A8*f8)/((beta2+beta12*Ec+2*betap2*N5)^2));
Ar7=(2*betap2*e3*A8*f8)*((beta2+beta12*Ec+2*betap2*N5)^3);
A44=((Ar1*Ar2)+(Ar3*Ar4)+(Ar5*Ar6)+Ar7)*q;

%Defining Elements of the Matix 'C' i.e. for 
%Unconstrained System for a Mechanical Recyling plant

C11=A11/k1;
C12=A12/k1;
C13=A13/k1;
C14=A14/k1;
C21=A21/k2;
C22=A22/k2;
C23=A23/k2;
C24=A24/k2;
C31=A31/k3;
C32=A32/k3;
C33=A33/k3;
C34=A34/k3;
C41=A41/k4;
C42=A42/k4;
C43=A43/k4;
C44=A44/k4;

dxpdt=[xp(5);
      xp(6);
      xp(7);
      xp(8);
     (C11)*xp(1)+(C12)*xp(2)+(C13)*xp(3)+(C14)*xp(4);
     (C21)*xp(1)+(C22)*xp(2)+(C23)*xp(3)+(C24)*xp(4);
     (C32)*xp(2)+(C33)*xp(3)+(C34)*xp(4);
     (C41)*xp(1)+(C42)*xp(2)+(C33)*xp(3)+(C44)*xp(4)];

function v=ewp_init(xp)
        v=[1.2;
           0.002;
           2;
           0.3;
           0;
           0;
           0;
           0]; 
 