close all; clear all;
M = dlmread('CRYSTALS_Si_Palik.csv.txt');
x = M(:,1);
lambda = 1000.*x; %wavelength in nanometers
length = size(x);
l= 10; %number of layers of hbn

n0 = ones(size(x));
n1 = (2.6 - 1.3i) .* ones(size(x)); %ref index of graphene
n2 = 1.85 .* ones(size(x)); %ref index of hbn
n3 = sqrt( 1 + 0.6961663*power(x,2)./(power(x,2)-power(0.0684043,2)) +(0.4079426*power(x,2))./(power(x,2)-power(0.1162414,2)) + (0.8974794*power(x,2))./(power(x,2)-power(9.896161,2)));  %refractive index of sio2
n4 = M(:,2) + 1i.*M(:,3); %ref index of silicon


r(:,1) = (n0-n1)./(n0+n1);
r(:,2) = (n1-n2)./(n1+n2);
r(:,3) = (n2-n3)./(n2+n3);
r(:,4) = (n3-n4)./(n3+n4);

d2 = 0.023;
%0.0 + l*0.0004; %size of hbn layer in micron + first term contaminant
d1 = 0.00034; %size of graphene monolayer micron
d3 = 0.290; %size of sio2 layer in microns

p(:,1) = exp(-2i*(2*pi*n1*d1)./x);  %graphene phase
p(:,2) = exp(-2i*(2*pi*n2*d2)./x);  %hbn phase
p(:,3) = exp(-2i*(2*pi*n3*d3)./x);  %Sio2 phase

g(:,4) = r(:,4);

for j= 1:3
    g(:,4-j) = (r(:,4-j) + g(:,5-j).* p(:,4-j))./ (1 + r(:,4-j).*g(:,5-j).*p(:,4-j));
end

I = (abs(g(:,1))).^2;


rn(:,1) = (n0 - n2)./(n0 + n2);
rn(:,2) = r(:,3);
rn(:,3) = r(:,4);
h(:,3) = rn(:,3);
q(:,2) = p(:,3);
q(:,1) = p(:,2);


for j=1:2   
    h(:,3-j) = (rn(:,3-j) + h(:,4-j).*q(:,3-j))./(1 + rn(:,3-j).*h(:,4-j).*q(:,3-j));
end
K = abs(h(:,1)).^2;


rnnn(:,1) = r(:,1);
rnnn(:,2) = (n1-n2)./(n2+n1);
rnnn(:,3) = r(:,4);
hh(:,3) = rnnn(:,3);
qq(:,2) = p(:,3);
qq(:,1) = p(:,1);


for j=1:2   
    hh(:,3-j) = (rnnn(:,3-j) + hh(:,4-j).*qq(:,3-j))./(1 + rnnn(:,3-j).*hh(:,4-j).*qq(:,3-j));
end
M = abs(hh(:,1)).^2;


rnn(:,1) = (n0-n3)./(n0+n3);
rnn(:,2) = r(:,4);
t = p(:,3);
f(:,2) = rnn(:,2);

f(:,1) = (rnn(:,1) + f(:,2).*t)./(1 + rnn(:,1).*f(:,2).*t);
L = abs(f(:,1)).^2;


C1 = (K - I)./K; %h+s - g+h+s
C2 = (L - K)./L; % s - h+s
C3 = (L - I)./L; %s - g+h+s
C4 = (M - I)./M; %g+s - g+h+s
C5 = (M - I - K + L)./M; %g+s - g+h+s - h+s + s
C6 = (L - M)./M; %s - g+s

%close all
%hold  on
figure
plot(lambda,C6,'r')
xlim([410 740]) 
% figure
% plot(lambda,C2)
% xlim([410 740]) 
% figure
% plot(lambda,C3,'g')
% figure
% plot(lambda,C4,'y')
% figuren4
% 
% plot(lambda,C5,'r')
% xlim([410 740]) 