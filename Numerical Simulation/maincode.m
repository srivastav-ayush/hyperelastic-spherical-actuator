% clc; close all; 
%%%%%%%%%%%%%%%%%%% Geometry constraints %%%%%%%%%%%%%%%%%%%%%%%
syms r 
A=0.3;        %Initial Inner Radius(A)
a=0.333668;   %Deformed Inner Radius(a)
q=0.609254;   %Deformed Outer Radius(q)
p=0.1;        %Internal Pressure(p)
omeg=80e-4;       %Rotational Velocity(omeg)
row=956;      %Density of Material(row)
R=(A^3+r^3-a^3)^(1/3);      %Incompressible constraint - Volume constant 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% Hyper elastic theory - matrix and invariants %%%%%%%%%%
L = r/R;
Lz=L;
F = [1/(L*Lz) 0 0;0 L 0;0 0 Lz];
B = F*(F.');
B1 = inv(B);
I1 = (B(1,1) + B(2,2) + B(3,3)); 
I2 = ((1/B(1,1)) + (1/B(2,2)) + (1/B(3,3)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%% Augmented Gent Model Coeeficient %%%%%%%%%%%%%%%%
% Im = 84; 
% A = 10.5;  B = 0.09; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%% Hyper elastic theory - calculations %%%%%%%%%%%%%%
syms b 
W_1 = -7/(54*(I1/81 - 28/27)); 
W_2 = 9/(200*I2^(1/2)); 
Im1 = [1 0 0;0 1 0;0 0 1];

str = -p*Im1 + 2*W_1*B - 2*W_2*B1;
sig_rr = str(1,1);
sig_tt = str(2,2);
sig_zz = str(3,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%% Differentiate sig_rr %%%%%%%%%%%%%%%%%%%%%% 
str_dif = (2*(sig_tt-sig_rr))/r;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%% Definite Integration %%%%%%%%%%%%%%%%%%%%%%
sigma_rr1 =  int(str_dif,r,a,b);  % lower limit == a; higher limit ==r;
sigma_rr = sigma_rr1;

RR = a:0.001:q;          
sigma_rr = subs(sigma_rr,b,RR); 
sigma_rr = -p + double(sigma_rr);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%% Rotational Calculations %%%%%%%%%%%%%%%%%%%%
ome_dif = pi*row*omeg*omeg*r;
sigma_rr2 = int(ome_dif,r,a,b);
RR1 = q:-0.001:a;     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%% Figure Plots %%%%%%%%%%%%%%%%%%%%%%%%%%
sigma_rr3 = subs(sigma_rr2,b,RR1);
sigma_rr3 = double(sigma_rr3);
sigma_rr4 = sigma_rr+sigma_rr3; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NDr = RR./q;
subplot(1,2,1)
plot(NDr,sigma_rr4,'b')
xlabel('$\frac{\textit{r}}{\textit{b}}$', 'FontName', 'Times New Roman', ...
       'FontSize',14,'Color','k', 'Interpreter', 'LaTeX')
ylabel('Radial Stress, \sigma_{r} [MPa]', 'FontName', 'Times New Roman', ...
       'FontSize',12,'Color','k')
hold on 
%%%%%%%%%%%%%%%%%%%%
sigma_tt = sig_tt-sig_rr; 
sigma_tt = sigma_rr + subs(sigma_tt,r,RR);
sigma_tt = double(sigma_tt)+sigma_rr3;

subplot(1,2,2)
plot(NDr,sigma_tt,'b')
xlabel('$\frac{\textit{r}}{\textit{b}}$', 'FontName', 'Times New Roman', ...
       'FontSize',14,'Color','k', 'Interpreter', 'LaTeX')
ylabel('Hoop Stress, \sigma_{\theta} [MPa]', 'FontName', 'Times New Roman', ...
       'FontSize',12,'Color','k')
hold on
%%%%%%%%%%%%%%%%%%%%
sigma_zz = sig_zz-sig_rr; 
sigma_zz = sigma_rr + subs(sigma_zz,r,RR);
sigma_zz = double(sigma_zz)+sigma_rr3;

% figure(3)
% plot(NDr,sigma_zz,'b')
% hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







%%%%%%%%%%%%%%%%%%%%%%%% Abaqus FEM Data %%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Deformed Radius 1 %%%%%
rad1 = [0.35476
0.36313
0.371618
0.380212
0.388904
0.397686
0.40655
0.41549
0.4245
0.433574
0.442708
0.451897
0.461137
0.470425
0.479757
0.48913
0.498542
0.50799
0.517472
0.526986
0.53653
0.546102
0.555701
0.565325
0.574974
0.584645
0.594339];
rad1 = rad1./q;
rad1 = transpose(rad1); 
%%%%% Radial Stress %%%%%
S_11 = [-0.0730792
-0.0649829
-0.0586172
-0.0528884
-0.047719
-0.0430428
-0.0388029
-0.0349503
-0.0314428
-0.0282434
-0.0253203
-0.0226453
-0.020194
-0.0179446
-0.0158781
-0.0139775
-0.012228
-0.010616
-0.00912984
-0.00775879
-0.00649324
-0.00532462
-0.0042452
-0.00324801
-0.0023268
-0.0014756
-0.000690276]; 
S_11 = transpose(S_11);
%%%%% Deformed Radius 1 %%%%%
rad2 = [0.338414
0.346517
0.35476
0.36313
0.371618
0.380212
0.388904
0.397686
0.40655
0.41549
0.4245
0.433574
0.442708
0.451897
0.461137
0.470425
0.479757
0.48913
0.498542
0.50799
0.517472
0.526986
0.53653
0.546102
0.555701
0.565325
0.574974
0.584645
0.594339
0.604053];
rad2 = rad2./q;
rad2 = transpose(rad2);
%%%%% Hoop Stress %%%%%
S_33 = [0.103021
0.0959775
0.0897061
0.0840939
0.0790502
0.0744977
0.070372
0.0666182
0.0631896
0.0600462
0.0571539
0.054483
0.0520079
0.0497061
0.0475584
0.0455476
0.0436589
0.041879
0.0401964
0.0386008
0.037083
0.0356349
0.0342495
0.03292
0.0316409
0.0304069
0.0292134
0.0280561
0.0269314
0.0258324]; 
S_33 = transpose(S_33); 

% S_22 = [0.124901
% 0.117404
% 0.110727
% 0.104767
% 0.0994264
% 0.0946253
% 0.0902948
% 0.0863772
% 0.0828228
% 0.079589
% 0.0766391
% 0.0739414
% 0.0714681
% 0.0691954
% 0.0671018
% 0.0651692
% 0.063381
% 0.061723
% 0.0601821
% 0.0587477
% 0.0574053]; 
% S_22 = transpose(S_22); 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%% Plotting Commands %%%%%%%%%%%%%%%%%%%%%%%
% subplot(1,2,1)
% plot(rad1,S_11,'o')
% hold on 
% 
% legend('Analytical, {\it{p}} = 0.1 MPa, {\omega} = 80 rpm','FEM, {\it{p}} = 0.1 MPa, {\omega} = 80 rpm')
% 
% subplot(1,2,2)
% plot(rad2,S_33,'o')
% hold on 
% 
% legend('Analytical, {\it{p}} = 0.1 MPa, {\omega} = 80 rpm','FEM, {\it{p}} = 0.1 MPa, {\omega} = 80 rpm')

% figure(3)
% plot(rad,S_22,'o')
% hold on 
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%