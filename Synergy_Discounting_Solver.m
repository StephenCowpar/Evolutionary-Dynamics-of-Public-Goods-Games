clear
%%% Set Params %%%%%%
sd = -0.05; %epsilon- small pos num for syn, small neg num for disc.
r = 7;%Multiplication Factor.
G = 6;%Group Size.
T = 1000;%max Time.
meshs= 1001;
mesht= 100001;
%%%%%%%%%%%%%%%%%%%%%%

ds=1/(meshs-1);
dt=T/(mesht-1);
s  = 0:ds:1;
t  = 0:dt:T;

 a = .1;
 %P0 =(2*a*(1+cos(a*s.^2)))./(3*a+cos(a)*sin(a));
 P0 = (2* a *s + 2*(1 - a)*(1-s));
 P0 = P0/trapz(s,P0);
%P0 = ones(1,meshs); %initial distribution.

sbar0 = trapz(s,s.*P0); %initial sbar.
s2bar0 = trapz(s,s.^(2).*P0); %initial s^2 bar.
P  = [P0 ; zeros(mesht-1,meshs)]; %Initialize solution.
sbar = [sbar0 zeros(1,mesht-1)]; %Initialize mean.
s2bar = [s2bar0 zeros(1,mesht-1)]; %Initialize s^2 mean.
norm = [trapz(s,P0) zeros(1,mesht-1)]; %Initialize norm tracking.

S = sd; 
D = dt;
B = r/G;
v = ones(1,meshs); %useful.

 for i=2:mesht
     P(i,:)= P(i-1,:)+D*P(i-1,:).*(((s-(v.*sbar(i-1))).*(B-1+B*(G-1)*S.*v.*sbar(i-1)))+...
         (B*(S/2).*(s.^2-(v.*s2bar(i-1)))));

     norm(i) = trapz(s,P(i,:)); %Compute total mass.
     P(i,:) = P(i,:)/norm(i); %Ensure it equals 1.

     sbar(i) = trapz(s,s.*P(i,:)); %Compute the averages for next step.
     s2bar(i) = trapz(s,s.^(2).*P(i,:));
 end

var = s2bar-sbar.^2;
s3bar = zeros(1,mesht);
s4bar = zeros(1,mesht);
s5bar = zeros(1,mesht);
for i=1:mesht 
    s3bar(i) = trapz(s,s.^(3).*P(i,:));
    s4bar(i) = trapz(s,s.^(4).*P(i,:));
    s5bar(i) = trapz(s,s.^(5).*P(i,:));
end
sig = s3bar -3*sbar.*var -sbar.^3;

rcrit = (2*G*var)./(2*(1+G*S*sbar).*var+S.*sig);

figure
surf(s(1:10:meshs),t(1:(mesht-1)/100:mesht),P(1:(mesht-1)/100:mesht,1:10:meshs))
figure
hold on
plot(t,sbar,'b')
plot(t,var,'k')
plot(t,sig,'r')
hold off
figure
plot(t,rcrit,'m',t,r*ones(1,mesht),'y')
hold off
