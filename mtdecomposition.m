function [ M, Miso, Mdc, Mclvd, Mdcdiag, Mclvddiag, isoper, dcper, clvdper, Mo, Mw, phi2, delta2, lambda2 ] = mtdecomposition( Cvo,angs )
%This funtion will make a moment tensor decompostion by using the formulas in
%the article of M.L. Jost and R. B. Herrmann A Student's Guide to and Review of Moment Tensors 1989 and
%Václav Vavry?uk  Focal mechanisms in anisotropic media 2005.
% Link of the article1: ftp://ftp.ingv.it/pub/elisa.tinti/articoli_source/Jost%201989.pdf
% Link of the article2: https://academic.oup.com/gji/article/161/2/334/557793/Focal-mechanisms-in-anisotropic-media

%Inputs of this function are Cvo=Elasticity Tensor with Voigt Notation in 
%PGa and angs vector which has [strike dip rake] angles.

%This function will give you the Moment Tensor (M), ISO part of the moment 
%tensor (Miso), DC part of the Moment Tensor (Mdc), CLVD part of the Moment 
%Tensor (Mclvd), Diagonalize version of the DC part of the Moment Tensor 
%(Mdcdiag), Diagonalize version of the CLVD part of the Moment Tensor 
%(Mclvddiag), ISO percentage in the Moment Tensor (isoper), DC percentage 
%in the moment tensor (dcper), CLVD percentage of the Moment Tensor 
%(clvdper), Seismic Moment (Mo), Moment Magnitude (Mw), New Strike (phi2), 
%New Dip (delta2) and New Rake  (lambda2) by using the DC part of the 
%Moment Tensor.

%Deniz Ertuncay, 2016
%Bogazýcý Unýversýty Kandilli Observatoy and Earthquake Research Instýtute
%Department of Geophysics

delta=angs(1,2)*pi/180; %Dip
phi=angs(1,1)*pi/180; %Strike
lambda=angs(1,3)*pi/180; %Rake

u=[cos(lambda)*cos(phi)+cos(delta)*sin(lambda)*sin(phi) cos(lambda)*sin(phi)-cos(delta)*sin(lambda)*cos(phi) -sin(delta)*sin(lambda)]; %slip vector
v=[-sin(delta)*sin(phi) sin(delta)*cos(phi) -cos(delta)]; %fault normal
d=[u(1,1)*v(1,1) u(1,2)*v(1,2) u(1,3)*v(1,3) u(1,2)*v(1,3)+u(1,3)*v(1,2) u(1,1)*v(1,3)+u(1,3)*v(1,1) u(1,1)*v(1,2)+u(1,2)*v(1,1)];

m=Cvo*d';
M=[m(1) m(6) m(5);m(6) m(2) m(4);m(5) m(4) m(3)];

[V,D]=eig(M);
%Finding the position of max,min and other values of eigenvalues of M.
%w1,w2 and w3 are for determining eigenvectors locations on V.
for i=1:3
        if D(i,i) == max([D(1,1) D(2,2) D(3,3)]);
            posmax=D(i,i) ; %Max eigenvalue is denoted by posmax
            w1=[i];
        elseif D(i,i) == min([D(1,1) D(2,2) D(3,3)]);
            posmin=D(i,i) ; %Min eigenvalue is denoted by posmin
            w2=[i];
        else
            posmed=D(i,i) ; %2nd eigenvalue is denoted by posmed
            w3=[i];
        end
end

%Finding first and second biggest eigenvalues in absulut sense. We can find
%Scalar Moment Magnitude by using them
if posmax > abs(posmin) && posmax > abs(posmed) && posmed > abs(posmin)
    absposmax=abs(posmax);
    absposmed=abs(posmed);
    absposmin=abs(posmin);
elseif posmax > abs(posmin) && posmax > abs(posmed) && posmed < abs(posmin)
    absposmax=abs(posmax);
    absposmed=abs(posmin);
    absposmin=abs(posmed);
elseif posmed > abs(posmax) && posmed > abs(posmin) && posmax > abs(posmin)
    absposmax=abs(posmed);
    absposmed=abs(posmax);
    absposmin=abs(posmin);
elseif posmed > abs(posmin) && posmed > abs(posmax) && posmin > abs(posmax)
    absposmax=abs(posmed);
    absposmed=abs(posmin);
    absposmin=abs(posmax);
elseif posmin > abs(posmax) && posmin > abs(posmed) && posmax > abs(posmed)
    absposmax=abs(posmin);
    absposmed=abs(posmax);
    absposmin=abs(posmed);
else
    absposmax=abs(posmin);
    absposmed=abs(posmed);
    absposmin=abs(posmin);
end

M0=0.5*(absposmax+absposmed); %Slacar Moment Magnitude
Mw=(2/3)*log(M0)-10.7;

D1=posmax; %1st eigenvalue
D2=posmed; %2nd eigenvalue
D3=posmin; %3rd eigenvalue

V1=[V(1,w1) V(2,w1) V(3,w1)]; %1st eigenvector
V2=[V(1,w3) V(2,w3) V(3,w3)]; %2nd eigenvector
V3=[V(1,w2) V(2,w2) V(3,w2)]; %3rd eigenvector
trace=[(D(1,1)+D(2,2)+D(3,3))./3 0 0;0 (D(1,1)+D(2,2)+D(3,3))./3 0;0 0 (D(1,1)+D(2,2)+D(3,3))./3];
Mstar=M-trace;
mstar=eig(Mstar);
[~,idx]=sort(abs(mstar)); %Sorting the mstar vector. 3rd number is the maximum
result=mstar(idx); %biggest value (absulute sense). 1st number is the minimum
mstar=result; %smallest number.
%Double Couple - CLVD - ISO Decomposition 
Mdcc1=trace; %Isotropic Part



[Vvav,Dvav]=eig(Mstar);
%Finding the position of max,min and other values of eigenvalues of M.
%w1,w2 and w3 are for determining eigenvectors locations on V.
for i=1:3
        if Dvav(i,i) == max([Dvav(1,1) Dvav(2,2) Dvav(3,3)]);
            posmax2=Dvav(i,i) ; %Max eigenvalue is denoted by posmax
            wvav1=[i];
        elseif Dvav(i,i) == min([Dvav(1,1) Dvav(2,2) Dvav(3,3)]);
            posmin2=Dvav(i,i) ; %Min eigenvalue is denoted by posmin
            wvav2=[i];
        else
            posmed2=Dvav(i,i) ; %2nd eigenvalue is denoted by posmed
            wvav3=[i];
        end
end
Vvav1=[Vvav(1,wvav1) Vvav(2,wvav1) Vvav(3,wvav1)]; %1st eigenvector
Vvav2=[Vvav(1,wvav3) Vvav(2,wvav3) Vvav(3,wvav3)]; %2nd eigenvector
Vvav3=[Vvav(1,wvav2) Vvav(2,wvav2) Vvav(3,wvav2)]; %3rd eigenvector



%Vavrycuk style Eigenvalues
if posmax2 > abs(posmin2) && posmax2 > abs(posmed2) && posmed2 > abs(posmin2)
    absposmax2=abs(posmax2);
    absposmed2=abs(posmed2);
    absposmin2=abs(posmin2);
elseif posmax2 > abs(posmin2) && posmax2 > abs(posmed2) && posmed2 < abs(posmin2)
    absposmax2=abs(posmax2);
    absposmed2=abs(posmin2);
    absposmin2=abs(posmed2);
elseif posmed2 > abs(posmax2) && posmed2 > abs(posmin2) && posmax2 > abs(posmin2)
    absposmax2=abs(posmed2);
    absposmed2=abs(posmax2);
    absposmin2=abs(posmin2);
elseif posmed2 > abs(posmin2) && posmed2 > abs(posmax2) && posmin2 > abs(posmax2)
    absposmax2=abs(posmed2);
    absposmed2=abs(posmin2);
    absposmin2=abs(posmax2);
elseif posmin2 > abs(posmax2) && posmin2 > abs(posmed2) && posmax2 > abs(posmed2)
    absposmax2=abs(posmin2);
    absposmed2=abs(posmax2);
    absposmin2=abs(posmed2);
else
    absposmax2=abs(posmin2);
    absposmed2=abs(posmed2);
    absposmin2=abs(posmin2);
end

eps=-mstar(1)/abs(mstar(3));

%   Herrmann's Way
% F=-(mstar(1)/mstar(3));
% if abs(posmax) >= abs(posmed) && abs(posmed) >= abs(posmin)
%     Mdcc2=mstar(3).*(1-2*F).*((V1'*V1)-(V2'*V2)); %Major Double Couple
%     Mdcc3=mstar(3).*F.*(2.*(V1'*V1)-(V2'*V2)-(V3'*V3)); %CLVD
% elseif abs(posmax) >= abs(posmin) && abs(posmin) >= abs(posmed)
%     Mdcc2=mstar(3).*(1-2*F).*((V1'*V1)-(V3'*V3)); %Major Double Couple
%     Mdcc3=mstar(3).*F.*(2.*(V1'*V1)-(V2'*V2)-(V3'*V3)); %CLVD
% elseif abs(posmed) >= abs(posmax) && abs(posmax) >= abs(posmin)
%     Mdcc2=mstar(3).*(1-2*F).*((V2'*V2)-(V1'*V1)); %Major Double Couple
%     Mdcc3=mstar(3).*F.*(2.*(V2'*V2)-(V3'*V3)-(V1'*V1)); %CLVD
% elseif abs(posmed) >= abs(posmin) && abs(posmin) >= abs(posmax)
%     Mdcc2=mstar(3).*(1-2*F).*((V2'*V2)-(V3'*V3)); %Major Double Couple
%     Mdcc3=mstar(3).*F.*(2.*(V2'*V2)-(V3'*V3)-(V1'*V1)); %CLVD
% elseif abs(posmin) >= abs(posmax) && abs(posmax) >= abs(posmed)
%     Mdcc2=mstar(3).*(1-2*F).*((V3'*V3)-(V1'*V1)); %Major Double Couple
%     Mdcc3=mstar(3).*F.*(2.*(V3'*V3)-(V2'*V2)-(V1'*V1)); %CLVD
% else
%     Mdcc2=mstar(3).*(1-2*F).*((V3'*V3)-(V2'*V2)); %Major Double Couple
%     Mdcc3=mstar(3).*F.*(2.*(V3'*V3)-(V1'*V1)-(V2'*V2)); %CLVD
% end

%   Vavrycuk's Way
if abs(posmax2) >= abs(posmed2) && abs(posmed2) >= abs(posmin2)
    Mdcc2=((1-2*abs(eps))*absposmax2).*((Vvav1'*Vvav1)-(Vvav2'*Vvav2)); %Major Double Couple
    Mdcc3=(abs(eps)*absposmax2).*(2.*(Vvav1'*Vvav1)-(Vvav2'*Vvav2)-(Vvav3'*Vvav3)); %CLVD
elseif abs(posmax2) >= abs(posmin2) && abs(posmin2) >= abs(posmed2)
    Mdcc2=((1-2*abs(eps))*absposmax2).*((Vvav1'*Vvav1)-(Vvav3'*Vvav3)); %Major Double Couple
    Mdcc3=(abs(eps)*absposmax2).*(2.*(Vvav1'*Vvav1)-(Vvav2'*Vvav2)-(Vvav3'*Vvav3)); %CLVD
elseif abs(posmed2) >= abs(posmax2) && abs(posmax2) >= abs(posmin2)
    Mdcc2=((1-2*abs(eps))*absposmax2).*((Vvav2'*Vvav2)-(Vvav1'*Vvav1)); %Major Double Couple
    Mdcc3=(abs(eps)*absposmax2).*(2.*(Vvav2'*Vvav2)-(Vvav3'*Vvav3)-(Vvav1'*Vvav1)); %CLVD
elseif abs(posmed2) >= abs(posmin2) && abs(posmin2) >= abs(posmax2)
    Mdcc2=((1-2*abs(eps))*absposmax2).*((Vvav2'*Vvav2)-(Vvav3'*Vvav3)); %Major Double Couple
    Mdcc3=(abs(eps)*absposmax2).*(2.*(Vvav2'*Vvav2)-(Vvav3'*Vvav3)-(Vvav1'*Vvav1)); %CLVD
elseif abs(posmin2) >= abs(posmax2) && abs(posmax2) >= abs(posmed2)
    Mdcc2=((1-2*abs(eps))*absposmax2).*((Vvav3'*Vvav3)-(Vvav1'*Vvav1)); %Major Double Couple
    Mdcc3=(abs(eps)*absposmax2).*(2.*(Vvav3'*Vvav3)-(Vvav2'*Vvav2)-(Vvav1'*Vvav1)); %CLVD
else
    Mdcc2=((1-2*abs(eps))*abs(mstar(3))).*((Vvav3'*Vvav3)-(Vvav2'*Vvav2)); %Major Double Couple
    Mdcc3=(abs(eps)*abs(mstar(3))).*(2.*(Vvav3'*Vvav3)-(Vvav1'*Vvav1)-(Vvav2'*Vvav2)); %CLVD
end


% Diagonal Matrix Form
if abs(posmax2) >= abs(posmed2) && abs(posmed2) >= abs(posmin2)
    Mdccc2=((1-2*abs(eps))*absposmax2).*[-1 0 0;0 0 0;0 0 1]; %Major Double Couple
    Mdccc3=(abs(eps)*absposmax2).*[-1 0 0;0 -1 0;0 0 2]; %CLVD
elseif abs(posmax2) >= abs(posmin2) && abs(posmin2) >= abs(posmed2)
    Mdccc2=((1-2*abs(eps))*absposmax2).*[-1 0 0;0 0 0;0 0 1]; %Major Double Couple
    Mdccc3=(abs(eps)*absposmax2).*[-1 0 0;0 -1 0;0 0 2]; %CLVD
elseif abs(posmed2) >= abs(posmax2) && abs(posmax2) >= abs(posmin2)
    Mdccc2=((1-2*abs(eps))*absposmax2).*[-1 0 0;0 0 0;0 0 1]; %Major Double Couple
    Mdccc3=(abs(eps)*absposmax2).*[-1 0 0;0 -1 0;0 0 2]; %CLVD
elseif abs(posmed2) >= abs(posmin2) && abs(posmin2) >= abs(posmax2)
    Mdccc2=((1-2*abs(eps))*absposmax2).*[-1 0 0;0 0 0;0 0 1]; %Major Double Couple
    Mdccc3=(abs(eps)*absposmax2).*[-1 0 0;0 -1 0;0 0 2]; %CLVD
elseif abs(posmin2) >= abs(posmax2) && abs(posmax2) >= abs(posmed2)
    Mdccc2=((1-2*abs(eps))*absposmax2).*[-1 0 0;0 0 0;0 0 1]; %Major Double Couple
    Mdccc3=(abs(eps)*absposmax2).*[-1 0 0;0 -1 0;0 0 2]; %CLVD
else
    Mdccc2=((1-2*abs(eps))*abs(mstar(3))).*[-1 0 0;0 0 0;0 0 1]; %Major Double Couple
    Mdccc3=(abs(eps)*abs(mstar(3))).*[-1 0 0;0 -1 0;0 0 2]; %CLVD
end


isoper=((1/3)*((D(1,1)+D(2,2)+D(3,3))/abs(absposmax2)))*100; %ISO percentage
clvdper=2*eps*(100-abs(isoper)); %CLVD percentage
dcper=100-abs(isoper)-abs(clvdper); %DC percentage
Mson=Mdcc1+Mdcc2+Mdcc3;


[V22,D22]=eig(Mdcc2);
for i=1:3
        if D22(i,i) == max([D22(1,1) D22(2,2) D22(3,3)]);
            posmax2=D22(i,i) ; %Max eigenvalue is denoted by posmax
            w11=[i];
        elseif D22(i,i) == min([D22(1,1) D22(2,2) D22(3,3)]);
            posmin2=D22(i,i) ; %Min eigenvalue is denoted by posmin
            w22=[i];
        else
            posmed2=D22(i,i) ; %2nd eigenvalue is denoted by posmed
            w33=[i];
        end
end
V221=[V22(1,w11) V22(2,w11) V22(3,w11)]; %1st eigenvector
V222=[V22(1,w33) V22(2,w33) V22(3,w33)]; %2nd eigenvector
V223=[V22(1,w22) V22(2,w22) V22(3,w22)]; %3rd eigenvector

uu1=(V221+V223)./norm(V221+V223);
uu2=(V221-V223)./norm(V221-V223);
uu3=(-V221-V223)./norm(-V221-V223);
uu4=(-V221+V223)./norm(-V221+V223);
uu=[uu1; uu2; uu3; uu4];
vv1=(V221-V223)./norm(V221-V223);
vv2=(V221+V223)./norm(V221+V223);
vv3=(-V221+V223)./norm(-V221+V223);
vv4=(-V221-V223)./norm(-V221-V223);
vv=[vv1; vv2; vv3; vv4];
for i=1:4
    dp(i)=dot(vv(i,:),v); %We want to find the closest vv to 1. 
    dp2(i)=dot(uu(i,:),u); %We want to find the closest uu to 1. 
end
J=find(dp==max(dp));
J2=find(dp2==max(dp2));
v2=vv(J(1),:);
u2=uu(J2(1),:);

delta2=acosd(-v2(3)); %Dip
% if delta2 > 90
%     delta2=180-delta2;
% end
%phi2=acosd(v2(2)/sind(delta2)); 
phi2=asind(-v2(1)/sind(delta2)); %Strike
% if phi2 < 0
%     phi2=360+phi2;
% end
lambda2=asind(-u2(3)/sind(delta2)); %Slip(rake)

%Naming the variables for the output
M;Miso=Mdcc1;Mdc=Mdcc2;Mclvd=Mdcc3;Mdcdiag=Mdccc2;Mclvddiag=Mdccc3;
isoper;dcper;clvdper;Mo=M0;Mw=Mw;phi2;delta2;lambda2;
end