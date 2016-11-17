printf('Taking Input data For Y,Z Graph:-\n');
data = load('yzdata.txt');
Y_fig = data(:,1);
Z_fig= data(:,2);
Z_fig = -1*Z_fig;

figure;
plot(Z_fig,Y_fig,'r-','Markersize',10); xlabel('Z-Axis'); ylabel('Y-Axis'); hold on; %fprintf('Tell the point where you want to break the beam:-'); %cross = load('yzcross.txt'); %Y_cross = cross(:,1); %Z_cross = cross(:,2); %p = cross(:,3); %m1 = length(Y_cross); %plot(Z_cross,Y_cross,'b-','Markersize',10);

%%second Approch ...here first approach is not working fprintf('\nTaking Length Input values for Y,Z:-\n'); length = load('length.txt'); Y = length(:,1); Z = length(:,2); yc = length(:,3); zc = length(:,4); %y1 = length(1,1); %y2 = length(2,1); %y3 = length(3,1); E = length(:,5); %z1 = length(1,2); %z2 = length(2,2); %z3 = length(3,2); A = Y.*Z;
%A1 = y1*z1; %A2 = y2*z2; %A3 = y3*z3;

fprintf('Calculating Netural Point:-\n');
y_center = (sum(yc.*A.*E))/sum(E'*A);
z_center = (sum(zc.*A.*E))/sum(E'*A);

fprintf('Y_center:-');
disp(y_center);
fprintf('Z_center:-');
disp(z_center);

plot(zc,yc,'bx','Markersize',10);
plot(z_center,y_center,'bo');
legend('Boundry of Beam','Centriod of a particular Section','Neutral point');

% Here We find Inertia
%first Of all I_YY fprintf('Calculation of I_yy;-\n'); I_yy1 = ((zc-z_center).^2).*A; I_yy2 = (Y.*(Z.^3))/12; I_yy = I_yy1 + I_yy2; I_yy = E.*I_yy.*10^-8; fprintf('I_yy Part Wise:-\n'); disp(I_yy);
tot_I_yy = sum(I_yy);
fprintf('Total I_yy');
disp(tot_I_yy);

%Calculation Of I_zz; fprintf('Calculation of I_zz:-\n'); I_zz1 = ((yc-y_center).^2).*A; I_zz2 = (Z.*(Y.^3))/12; I_zz = I_zz1 + I_zz2; I_zz = E.*I_zz.*10^-8; fprintf('I_zz Part Wise:-\n'); disp(I_zz);
tot_I_zz = sum(I_zz);
fprintf('Total I_zz');
disp(tot_I_zz);

%Calculation Of I_yz; fprintf('Calculation of I_yz;-\n'); I_yz1 = ((yc-y_center).*(zc-z_center)).*A; I_yz2 = 0; I_yz = I_yz1 + I_yz2; I_yz = E.*I_yz.*10^-8; fprintf('I_yz Part Wise:-\n'); disp(I_yz);
tot_I_yz = sum(I_yz);
fprintf('Total I_yz');
disp(tot_I_yz);

fprintf('Inertia Matrix:-\n');
I = [tot_I_yy tot_I_yz; tot_I_yz tot_I_zz;];

disp(I);

fprintf('Principal Inertia Matrix:-\n');
[V,D] = eig(I);
disp(V);
disp(D);

B = [I(2,1) I(1,1); -I(2,2) -I(1,2);];
prompt = 'Enter Value of My:-';
My = input(prompt);
prompt = 'Enter Value of Mz:-';
Mz = input(prompt);

M = [My ; Mz;];
C = (1/det(I)).*B*M*10^-2;

disp(C);
syms y z;
sigmaXX = -((y)*C(1)*10^-2 + (z)*C(2)*10^-2);
sigmaXX = eval(sigmaXX);
fprintf('Value of SigmaaXX:-\n');
disp(sigmaXX);
fprintf('Do You want to find sigmaaXX at any Particular Value of (Y,Z)??[Y/N]\n');
prompt = '[Y/N]:-';
check = input(prompt);

if(check == Y) prompt = 'Enter Value for y:-'; y_test = input(prompt); prompt ='Enter Value for z:-'; z_test = input(prompt);
    y= y_test-y_center;
    z = z_test-z_center;
    sigmaXX_test = subs(sigmaXX);
    sigmaXX_test = eval(sigmaXX_test);
    fprintf('Value for SigmaXX:-');
    disp(sigmaXX_test);
    syms y z;
else
    pause;
end
