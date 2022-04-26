clear
%%% Set Params %%%%%%
T = 100;
meshs= 101;
mesht= 100001;
ds=1/(meshs-1);
dt=T/(mesht-1);
s  = 0:ds:1;
t  = 0:dt:T;

S = 0.0001;  %sigma
B = 1;  %beta

a = .1;
P0 = exp(-((-(1/2)+s).^2/a^2))/(a.*sqrt(pi).*erf(1/(2*a)));
P0 = P0/trapz(s,P0);
%a = 30;
%P0 =(2*a*(1+cos(a*s.^2)))./(3*a+cos(a)*sin(a));
%P0 = P0/trapz(s,P0);
%P0 = ones(1,meshs); %initial distribution % Uniform
%P0 = [1.5*ones(1,(meshx-1)/2) 1 .5*ones(1,(meshx-1)/2)];
sbar0 = trapz(s,s.*P0); %initial xbar

A1 = S*dt/ds^2;
A2 = B*dt;

P  = [P0 ; zeros(mesht-1,meshs)];
sbar = [sbar0 zeros(1,mesht-1)];
norm = [1 zeros(1,mesht-1)];

 for i=2:mesht
     P(i,2:meshs-1)=P(i-1,2:meshs-1)+A1*(P(i-1,3:meshs)-2*P(i-1,2:meshs-1)+P(i-1,1:meshs-2))+...
         A2*(s(2:meshs-1)-sbar(i-1)).*P(i-1,2:meshs-1);
     P(i,1)=P(i-1,1)+2*A1*(P(i-1,2)-P(i-1,1))+ A2*(-sbar(i-1)).*P(i-1,1);
     P(i,meshs)=P(i-1,meshs)+2*A1*(P(i-1,meshs-1)-P(i-1,meshs))+ ...
         A2*(1-sbar(i-1)).*P(i-1,meshs);

     norm(i) = trapz(s,P(i,:)); %Compute total mass.
     P(i,:) = P(i,:)/norm(i); %Ensure it equals 1.

     sbar(i) = trapz(s,s.*P(i,:));
 end
 figure
 surf(s(1:10:meshs),t(1:1000:mesht),P(1:1000:mesht,1:10:meshs))