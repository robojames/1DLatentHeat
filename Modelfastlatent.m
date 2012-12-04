%% 1-D multi layer conduction transient control volume code
% with emphasis on thermoelectric heating and cooling 
% Written by Logan Compton
% 2-16-2012

% Code pulls existing electrical current data from experimental results and
% models the temperature of the top node of the thermal electric module
% which consists of layers of ceramic alumnia, alumnium, bismuth telluride,
% alumnium and ceramic alumnia. Applys peltier heating and cooling at
% alumnium and bisumth telluride interfaces along with joule heating within
% the bismuth telluride piece. Constant boundary condition on the bottom end
% node to represent cold block and heat transfer coefficient on the top
% node to represent natural convection. 
 

clear all
clc
tic
%% Build model by adding thermodynamic properties into array

%% layer properties 
k_SENSOR = 92.45;
cp_SENSOR = 450;
rho_SENSOR = 8049;

% Layering Scheme: ceramic, copper, bismuth telluride, copper, ceramic, copper, glass,
% bismuth telluride, copper, water, Type E themoelectric

%k=[35,223,1.48,223,35,1.4,223,1.2,223,(21.9+16.7),223,0.69]; %(W/mK)Thermal conductivity of layers (1,2,3..) 
k=[35,223,1.48,223,35,1.4,223,1.2,223, (21.9+16.7), 0.69]; %(W/mK)Thermal conductivity of layers (1,2,3..) 

%rho=[3750,8800,7700,8800,3750,2225,8800,7700,8800,8730,8800,1000];  % (kg/m^3)density of layers (1,2,3..)
rho=[3750,8800,7700,8800,3750,2225,8800,7700,8800, 8730,1000];  % (kg/m^3)density of layers (1,2,3..)

%cp=[775,420,122,420,775,835,420,122,420,(4.07e-1+3.984e-1),420,4218.3]; % (J/kg*K)heat capacity of layers (1,2,3..) 499
cp=[775,420,122,420,775,835,420,122,420, (4.07e-1+3.984e-1),4218.3]; % (J/kg*K)heat capacity of layers (1,2,3..) 499

len=[0.0006858,0.0004318,0.00127,0.0004318,0.0006858,0.0001524,0.0001524,500e-6,1e-6,25e-6,250e-6]; %(m)length of layers (1,2,3..)  %(m)length of layers (1,2,3..) 
%nodes=[round(len(1)/cvl),round(len(2)/cvl),round(len(3)/cv2),round(len(4)/cvl),round(len(5)/cvl),round(len(6)/cv3),round(len(7)/cv4)]; %nodes of layers (1,2,3..) 

 nodes=[500,600,1200,600,500,500,500,800,300,400,600];%mesh independent volume 


%% Experimental Model Properties
  Jo=xlsread('1Kpersec','E1:E1200'); % retreives current profile 
  t=xlsread('1Kpersec','B1:B1200'); % retreives time profile 
   T_init=xlsread('Tinit','A1:A7004'); % retreives intial temperature distrubtion within the model 
  
  %% Specify Time step
  Dt=0.005; % constant time step [s]
  
 %%%%%%%%%%%%%%%%%%%%%calculation insures even time step%%%%%%%%%%%%%%%%%%%%%%
  t1=diff(t);% data time step [s]
  ts=round(t1/Dt);  % number of timing loops per data time step
  n_time=sum(ts);   % total number of timing loops per data time step. 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 %% Specify Initial Temperature
  %T_init = 288; % initial temperature */
  J=(Jo./(2.17e-6))';%* electric current density in [A/m2] 
 rho_e = [(2*8.82e-6),1.68e-8];%/* electric resistivity in [ohm-m] */
 %rho_e = [(1.026e-3),1.68e-8];
 alpha =0.378*1.92e-4; %/* Seebeck coefficient [V/K] */
  
%% Internal Heat generation & Boundary Conditions 
spa=[0,0,0,0,0,0,0,0,0,0,0,0];
sca=[0,0,0,0,0,0,0,0,0,0,0,0];
layerbc=[0,1,-1,0,0,0,0,0,0,0,0];% applies peltier cooling/heating at layer interfaces -1=cooling 1=heating 
control = [0,1]; %/* control should be 0 if a temperature BC is used and 1 otherwise */
T_bc =[251,0]; %Temperature Boundary Conditions 
q = [0,0]; % conduction B.C. 
h = [0,0.5]; % heat transfer coefficient for convection B.C. 
T_inf = [0,288]; % ambient air temperature for convection B.C. 

%% Natural convection information 
L=(1.95e-6)/(0.0056); %charcteristic length (A/P)
Beta=1/T_inf(2); %
Nav=14.54e-6; 
Naalpha=20.62e-6;
Nak=24e-3; 

%% Latent Heat 
Tm=273-23; %Freezing Temperature K
Lm=0; % Latent heat of melting J/Kg;
cs=211.0; % specific heat of ice (J/KgK); 
cl=211.0; % specific heat of water (J/Kg/K); 
flag=0; 
count=0; 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%Control Volume Creation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% Set up a numbering scheme to make node indexing easier */
n(1)=nodes(1);
for i=2:length(nodes)
    n(i)=n(i-1)+nodes(i); 
end
% Set up a length scheme to make position indexing easier */
l(1)=len(1); 
for i=2:length(nodes)
    l(i)=l(i-1)+len(i); 
end
% control volume
cv(1)=2.0*(len(1)/(2.0*nodes(1) + 1));%/* control volume for first layer */
cv(length(nodes))=2.0*(len(length(nodes))/(2.0*nodes(length(nodes)) + 1)); %/* control volume for last layer */
for i=2:length(nodes)-1
    cv(i)=len(i)/nodes(i);  %/* control volume for ith layer */
end

%% assigning node position and intial temperature for each node. 
x(1,1) = 0;
T(1,1) = T_init(1);
x(1,n(length(nodes)) + 1) = l(length(nodes)-1) + len(length(nodes));
T(1,n(length(nodes)) + 1) = T_init(length(T_init));
for i=2:nodes(1)
		x(1,i) = x(1,i-1) + cv(1);
        T(1,i) = T_init(i);
end
     
    for j=2:length(nodes)
	        for i = 1:nodes(j) 
		x(1,n(j-1) + i) = l(j-1) + cv(j)/2.0 + cv(j)*(i-1);
        T(1,n(j-1) + i) = T_init(n(j-1) + i);
        
            end
    end
    
    %% Pre Allocating memory

    AP=zeros(1,(n(length(nodes))+1)); 
    AE=zeros(1,(n(length(nodes))+1)); 
    AW=zeros(1,(n(length(nodes))+1)); 
    Ao=zeros(1,(n(length(nodes))+1)); 
    sp=zeros(1,(n(length(nodes))+1)); 
    sc=zeros(1,(n(length(nodes))+1)); 
    b=zeros(1,(n(length(nodes))+1)); 
    P=zeros(1,(n(length(nodes))+1)); 
    Q=zeros(1,(n(length(nodes))+1)); 
    f=zeros(1,length(nodes)); 
    k_eq=zeros(1,length(nodes)); 
    Tlastnode=zeros(n_time,1); 
    
    
    Tlastnode(1,:)=T_init(length(T_init)); 
    
    %% Start of timing loop
	z=1; % z denotes which timing loop (ts)
    s=1; % s denotes which step inside current timing loop (ts)
  for time = 1:n_time

if s==ts(z) 
    z=z+1; % if s is equal to the total number of steps inside timing loop  
    s=0;  % then progress towards next timing loop and reset s
end
clc


%disp('percent complete')
if (mod((z/length(t)),2)==0)
    disp(z/length(t)*100) 
end
s=s+1;


%% Natural Heat Transfer Coefficient 
Ra=abs(real((9.81*Beta*(T(1,n(length(nodes))+1)-T_inf(2))*L^3)/(Nav*Naalpha))); %Rayleigh number
Nul=0.27*Ra^(1/4); %Nusselt number
Nah=(Nul*Nak)/L; % heat transfer coefficient for natural convection
h=[0,Nah]; 
JH=J(z)*J(z).*rho_e;
sca=[0,JH(2),JH(1),JH(2),0,0,0,0,0,0,0,0];

%% Layering 
sc(1) = sca(1);
sp(1) = spa(1);
AE(1) = control(1)*k(1)/cv(1);
AW(1) = 0.0;
Ao(1) = control(1)*rho(1)*cp(1)*cv(1)*0.5/Dt;
AP(1) = AE(1) + AW(1) +  Ao(1) + h(1) - sp(1)*cv(1)*0.5*control(1) + 1.0*(1.0 - control(1));
b(1) = T_bc(1) + q(1) + h(1)*T_inf(1) + Ao(1)*T(1,1) + sc(1)*cv(1)*0.5*control(1);
      
     
        for o=2:nodes(1)
               sc(o) = sca(1);
		       sp(o) = spa(1);
			   AE(o) = k(1)/cv(1);
			   AW(o) = k(1)/cv(1);
			   Ao(o) = rho(1)*cp(1)*cv(1)/Dt;
			   AP(o) = AE(o) + AW(o) + Ao(o) - sp(o)*cv(1);
			   b(o) = Ao(o)*T(1,o) + sc(o)*cv(1);
        end
           
   for p=1:length(nodes)-1
       for i = 2:nodes(p+1)
                         
           sc(n(p) + i) = sca(p);
               sp(n(p) + i) = spa(p);
			   AE(n(p) + i) = k(p+1)/cv(p+1);
			   AW(n(p) + i) = k(p+1)/cv(p+1);
			   Ao(n(p) + i) = rho(p+1)*cp(p+1)*cv(p+1)/Dt;
			   AP(n(p) + i) = AE(n(p) + i)+ AW(n(p) + i) + Ao(n(p) + i) - sp(n(p) + i)*cv(p+1);
			   b(n(p) + i) = Ao(n(p) + i)*T(1,n(p) + i) + sc(n(p) + i)*cv(p+1);
       end
            
         f(p) = (l(p) - x(n(p)))/(x(n(p)+1) - x(n(p)));
         k_eq(p) = 1.0/((1 - f(p))/k(p) + f(p)/k(p+1));
         
         AE(n(p)) = k_eq(p)/(0.5*(cv(p) + cv(p+1)));
         AP(n(p)) = AE(n(p)) + AW(n(p)) + Ao(n(p)) - sp(n(p))*cv(p);
       
       
       sc(n(p) + 1) = sca(p);
	   sp(n(p) + 1) = spa(p);
	   AE(n(p) + 1) = k(p+1)/cv(p+1);
	   AW(n(p) + 1) = k_eq(p)/(0.5*(cv(p) + cv(p+1)));
	   Ao(n(p) + 1) = rho(p+1)*cp(p+1)*cv(p+1)/Dt;
	   AP(n(p) + 1) = AE(n(p) + 1) + AW(n(p) + 1) + Ao(n(p) + 1) - sp(n(p) + 1)*cv(p + 1);
	   b(n(p) + 1) = Ao(n(p) + 1)*T(1,n(p) + 1) + sc(n(p) + 1)*cv(p + 1); 
   
   end
   
 %/* take care of first and last boundary node coefficients,
            %based on type of B.C. */
sc(1) = sca(1);
sp(1) = spa(1);
AE(1) = control(1)*k(1)/cv(1);
AW(1) = 0.0;
Ao(1) = control(1)*rho(1)*cp(1)*cv(1)*0.5/Dt;
AP(1) = AE(1) + AW(1) +  Ao(1) + h(1) - sp(1)*cv(1)*0.5*control(1) + 1.0*(1.0 - control(1));
b(1) = T_bc(1) + q(1) + h(1)*T_inf(1) + Ao(1)*T(1,1) + sc(1)*cv(1)*0.5*control(1);
      
 
j=(n(length(nodes))+1); m=length(nodes); 
sc(j) = sca(m);
sp(j) = spa(m);
 AW(j) = control(2)*k(m)/cv(m);
 AE(j) = 0.0;
 Ao(j) = control(2)*rho(m)*cp(m)*cv(m)*0.5/Dt;
 AP(j) = AE(j) + AW(j) +  Ao(j) + h(2) - sp(j)*cv(m)*0.5*control(2) + 1.0*(1.0 - control(2));
 b(j) = T_bc(2) + q(2) + h(2)*T_inf(2) + Ao(j)*T(1,j) + sc(j)*cv(length(nodes))*0.5*control(2);
 
 
%%   applies peliter cooling to nodes around the interface layer
 for a=1:length(layerbc)
  if layerbc(a) ~=0

      qv1(time)=(alpha*J(z)*T(1,n(a)))/(cv(a));

sc(n(a)+(1)) = qv1(time)*layerbc(a); % Heating/Cooling  at the direct interface
b(n(a) + 1) = Ao(n(a) + 1)*T(1,n(a) + 1) + sc(n(a) + 1)*cv(a+1); 

sc(n(a)-(1)) = qv1(time)*layerbc(a); % Heating/Cooling  at the direct interface
b(n(a) - 1) = Ao(n(a) - 1)*T(1,n(a) - 1) + sc(n(a) - 1)*cv(a - 1);

sc(n(a)) = qv1(time)*layerbc(a); % Heating/Cooling  at the direct interface
b(n(a)) = Ao(n(a))*T(1,n(a)) + sc(n(a))*cv(a);

  end
 end

 %% Latent Heat
 p1=10; %node indexing of water only
 p2=10;
 Lim=0.2; %limits
 Mp=100; % Skewness of relaxation profile, Mp=1: Linear Line, Mp=1000: square wave
 for i = 2:nodes(p1+1)
 Twat(i)=T(1,n(p1) + i);% Temperature of Water at each node
 end
 for kkk=nodes(p1+1)+1:nodes(p2+1)
 Twat(kkk)=T(1,n(p2) + kkk);
 end
 Twater=mean(Twat); % average temperature of Water
 
 if Twater<=Tm+Lim && Twater>=Tm-Lim && flag==0
 
 YL=atan(Mp*-Lim); 
 YU=atan(Mp*Lim)-YL; 
 fr=abs(atan(((Tm-Twater)*Mp)-YL)*(1/YU)); %relaxation factor  
 dfdT(time)=(Mp/((4*Mp*Twater^2)-(8*Tm*Twater)+(4*Mp*Tm^2)+1)); %derivative of relaxation factor relative to Temperature 
 M=rho(10)*(cs+(fr*(cl-cs))+(dfdT(time)*(((cs-cl)*Tm)+Lm)));
     if flag==0
     Lat=Lm*rho(10); 
     Aop1=M*cv(p1+1)/Dt;
     end
 if flag==1
     cp(10)=2110; %specific heat of ice
     rho(10)=919; %density of ice
     k(10)=2.2; %thermal conductivity of ice
     Lat=0;
      
 end
 flag=1; 
 
 
 for i = 2:nodes(p1+1)
      
       sc(n(p1) + i) =Lat; 
       sp(n(p1) + i) =0;
               AE(n(p1) + i) = k(p1+1)/cv(p1+1);
			   AW(n(p1) + i) = k(p1+1)/cv(p1+1);
			   Ao(n(p1) + i) = Aop1; 
               AP(n(p1) + i) = AE(n(p1) + i)+ AW(n(p1) + i) + Ao(n(p1) + i) - sp(n(p1) + i)*cv(p1+1);
			   b(n(p1) + i) = Ao(n(p1) + i)*T(1,n(p1) + i) + sc(n(p1) + i)*cv(p1+1);
 end
%  for i = 2:nodes(p2+1)
%        sc(n(p2) + i) =Lat;
%        sp(n(p2) + i) =0;
%                AE(n(p2) + i) = k(p2+1)/cv(p2+1);
% 			   AW(n(p2) + i) = k(p2+1)/cv(p2+1);
% 			   Ao(n(p2) + i) = Aop1;
% 			   AP(n(p2) + i) = AE(n(p2) + i)+ AW(n(p2) + i) + Ao(n(p2) + i) - sp(n(p2) + i)*cv(p2+1);
% 			   b(n(p2) + i) = Ao(n(p2) + i)*T(1,n(p2) + i) + sc(n(p2) + i)*cv(p2+1);
%  end
 end

 if flag==1
     cp(10)=2110; %specific heat of ice
     rho(10)=919; %density of ice
     k(10)=2.2; %thermal conductivity of ice
 end
   %/* Solve for the temperature field using TDMA.  Move in positive x dir. */

		 P(1) = AE(1)/AP(1);
		 Q(1) = b(1)/AP(1);

		 for i = 2:n(length(nodes))+1 
			P(i) = AE(i)/(AP(i) - AW(i)*P(i-1));
			Q(i) = (b(i) + AW(i)*Q(i-1))/(AP(i) - AW(i)*P(i-1));
         end

		 T(1, n(length(nodes)) + 1) = Q(n(length(nodes)) + 1);

		 for i = n(length(nodes)):-1:1
			T(1,i) = P(i)*T(1,i+1) + Q(i);
         end

%Tlastnode(time+1,:)=T(1,n(length(nodes))+1); 
Tlastnode(time+1,:)=T(1,5700); 
      
  end
  % creates time array for plotting 
  ttime(1)=0;
  for q=2:n_time+1
      ttime(q)=ttime(q-1)+Dt;
  end
  
 
  figure(1)
  hold on
  plot(ttime,Tlastnode)
  xlabel('time (s)')
  ylabel('Temperature of last node') 
  hold off
 
toc
MM=[ttime',Tlastnode]; 
xlswrite('Data_legitbaseline',MM)