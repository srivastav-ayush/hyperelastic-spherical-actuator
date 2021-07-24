%%%%% Internal Pressure and Deformed Inner Radius %%%%%
P = [Removed for Confidentiality];

Def_R = [Removed for Confidentiality];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%% Exponential Curve Fitting %%%%%%%%%%%%%%%
a = Removed for Confidentiality;
b = Removed for Confidentiality;
c = Removed for Confidentiality;
d = Removed for Confidentiality;
x = [Removed for Confidentiality];
y = a*exp(b*x) + c*exp(d*x);
plot(x,y)
hold on
plot(P,Def_R,'o')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
