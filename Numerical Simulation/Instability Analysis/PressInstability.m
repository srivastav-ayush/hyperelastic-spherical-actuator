%clc; clear all; close all; 

A=  10.5;

B=  0.09;

Im = 100.8; 

syms L

i = 1;

for Lb=1:1e-3:3.8

A1=0.3;  B1=0.5; 

Lb;

La=((1+(B1./A1).^3.*(Lb.^3-1))).^(1./3);

fun = @(L)(2.*(((((2.*A)./(Im-(1./(L.^4) + 2.*(L.^2)))).*(L.^2-(1./(L.^4))))-((B./sqrt(L.^4 + 2./(L.^4))).*((1./(L.^2))-L.^4)))./(L.*(1-L.^3))));


q_AG = integral(fun,La,Lb);

store(i)=q_AG;
store1(i) =Lb;
store2(i)=La;
i=i+1;
end


%%%%%%%%%%%%%%%% FEM data %%%%%%%%%%%%%%%%
Pressure =[0.02
0.04
0.06
0.08
0.1
0.12
0.14
0.16
0.18];

La1=[1.0252
1.053363333
1.085233333
1.121913333
1.16512
1.217766667
1.285596667
1.383126667
1.5887]; 

Lb1=[1.00555
1.012008
1.019636
1.028836
1.040238
1.05495
1.07519
1.106718
1.181696]; 


figure(1)
subplot(1,2,1)
plot(store2,store)
hold on
plot(La1,Pressure,'o')
subplot(1,2,2)
plot(store1,store)
hold on
plot(Lb1,Pressure,'o')







