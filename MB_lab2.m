%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Montasir Bechir Nasir
%10081129
%V.01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('goph_547_lab_2_data_w2016_rev.mat');


base_point=grav_survey_data(:,1);
NumRows=size(grav_survey_data,1);
day=grav_survey_data(:,2);
total_time=grav_survey_data(:,3);
total_time=grav_survey_data(:,4);
x=grav_survey_data(:,5);
y=grav_survey_data(:,6);
z=grav_survey_data(:,7);
g=grav_survey_data(:,8);
dg_tide=grav_survey_data(:,9);

%Part 1
%create plots of relative gravity vs total time
g_case2=g;
g_case2(1)=0;
figure;
plot(total_time,g_case2,'r');
hold on;
%superimpose plot of tidal variation vs total time on the previous figure 
plot(total_time,dg_tide,'b');
title('Relative gravity and Tidal Variations Versus Time Graph');
xlabel('time(s)');
ylabel('Relative gravity(MicroGal)');
legend('relative gravity','tidal gravity variations');
%create a bestfit line for relative gravity readings 
g_polyfit=polyfit(total_time,g_case2,20);


%Part 2
%create a contour plot of the terrain survey data using Xt, Yt, Zt 
figure;
contourf(Xt,Yt,Zt);
colorbar;
title('Terrain Contour Plot of the Survey Data');
xlabel('Easting (km)');
ylabel('Northing (km)');
z_datum = mean(mean(Zt));
Zt_range = range(Zt);

%Part 3 
%create a sorted matrix containing x,y,z data 
x_sort=[x,y,z];
%sort through x, and then y and save the indeces 
[x_sort,ind_sort]=sortrows(x_sort,[1,2]);

%Part 4
%extract only the unique indeces from the looped input data 
[x_sort,ind_uniq]=unique(x_sort,'rows');

%Part 5
%save the unique sorted values 
Xg=x_sort(:,1);
Yg=x_sort(:,2);
Zg=x_sort(:,3);
Nx=sqrt(2601);
Ny=Nx;

%Part 6
% Reshape matrices
Xg = reshape(Xg, Nx, Ny);
Yg = reshape(Yg, Nx, Ny);
Zg = reshape(Zg, Nx, Ny);
%do not forget to check meshgrid case

%Part 7
%initialize corrected data vector
g_corr=g;
%saving the raw data
g_raw=g_corr;
%overwrite g_raw to obtain absolute gravity effect readings for all points
for i = 2:length(g_raw)
    g_raw(i) = g_raw(1) + g_raw(i);
end

%overwrite g_raw to sort and cont contain unique values
g_raw=g_raw(ind_sort);
g_raw=g_raw(ind_uniq);
%create a grid of g_raw
g_raw=reshape(g_raw,Nx,Ny);

%create a contourplot of the absolute gravity
figure;
contourf(Xg,Yg, g_raw);
colorbar;
title('Absolute Gravity Contour Plot');
xlabel('Easting(km)');
ylabel('Northing(km)');

%part 8
%compute gt using the IGF formula for ground surveys
gt=9.7803278*(1+(.0053024*((sin(49.1286))^2))-(.0000058*((sin(2*49.1286))^2)));
%unit conversion to microgal
gt=gt*(1*10^8);
g_corr(1)=g_corr(1)-gt;
leftOver=g_corr(1);
%1gal=1cm/s^2

%Part 9
%create a conturplot with normal gravity correction
g_norm=g_corr;
for i=2:length(g_norm)
    g_norm(i)=g_norm(i)+g_norm(1);
end

g_norm=g_norm(ind_sort);
g_norm=g_norm(ind_uniq);
g_norm=reshape(g_norm,Nx,Ny);
figure;
contourf(Xg,Yg, g_norm);
colorbar;
title('Normal Correction Contour Plot');
xlabel('Easting(km)');
ylabel('Northing(km)');

%DRIFT CORRECTIONS
%Part 10

for i=1:length(g_corr)
    g_corr(i)=g_corr(i)-dg_tide(i);
end

%Part 11
%tidal drift corrections
g_tide=g_corr;
for i=2:length(g_tide)
    g_tide(i)=g_tide(i)+g_tide(1);
end

g_tide=g_tide(ind_sort);
g_tide=g_tide(ind_uniq);
g_tide=reshape(g_tide,Nx,Ny);

figure;
contourf(Xg,Yg, g_tide);
colorbar;
title('Tidal Drift Corrected Contour plot');
xlabel('Easting(km)');
ylabel('Northing(km)');

%Part 12
%Instrument drift corrections 
h=1;
for i=1:length(base_point)
    if base_point(i)==2 %end of current loop
        drift=g_corr(i)
        %compute dt
        t_2=total_time(i);
        t_1=total_time(h);
        dt=t_2-t_1;
        %compute drift rate for this loop
        drift_rate=drift/dt;
        %correct for instrument drift only 
        g_corr(h+1:i-1)=g_corr(h+1:i-1)-drift_rate*(total_time(h+1:i-1)- total_time(h));
        
        %subtract this from all future points
        g_corr(i:length(base_point))= g_corr(i:length(base_point))-drift;
    end
    if i>1 & grav_survey_data(i,1)==1
        h=i;
        drift=g_corr(h);
        g_corr(h:length(base_point))=g_corr(h:length(base_point))-drift;
    end
end

IsTrue=0;
for i=1:length(base_point)
    if base_point(i)==1 | base_point(i)==2
        if g_corr(i)==0
            IsTrue=IsTrue;
        end
        if g_corr(i)~0
            IsTrue=IsTrue+1;
        end
    end
        
end

for i=2:length(g_corr)
    g_corr(i)=g_corr(i)+g_corr(1);
end

%Time series plot of gcor against time
figure;
plot(total_time,g_corr);
title('Time Series Plot of Drift Corrected Gravity vs Time');
xlabel('total time(s)');
ylabel('(microgal)');


%Part 14

g_corr=g_corr(ind_sort);
g_corr=g_corr(ind_uniq);
g_corr=reshape(g_corr,Nx,Ny);

g_drift=g_corr;

figure;
contourf(Xg,Yg, g_drift);
colorbar;
title('Drift Corrected Gravity Contour Plot');
xlabel('Easting(km)');
ylabel('Northing(km)');

%Elevation and terrain correction
%Part 15
G=6.674*10^-11;
Me=5.972*10^24;
Re=6371000;
dZ=Zg-z_datum;

%dg_FA=-(G*Me*dZ)/(
dg_FA=size(Zg);
for i=1:length(Zg)
    for j=1:length(Zg)
        dg_FA(i,j) = -.3086*dZ(i,j)*1000; %~in mgal/m
    end
end

%convert dg_FA from m/s^2 to microgal
dg_FA=dg_FA.*10^8;

for i=1:length(dg_FA)
    for j=1:length(dg_FA)
        g_corr(i,j)=g_corr(i,j)+dg_FA(i,j);
    end
end

%Part 16
%create a free air correction contour plot
g_FA=g_corr;
figure;
contourf(Xg,Yg, g_FA);
colorbar;
title('Free Air Corrected Gravity Contour Plot');
xlabel('Easting(km)');
ylabel('Northing(km)');

%Part 17
%Bouguer plate correction for surface
rhoB=2.65*(100^3)*(1/1000);%g/cc
dg_BP=size(dZ);
for i=1:length(dZ)
    for j=1:length(dZ)
        dg_BP(i,j)=2*pi*G*rhoB.*dZ(i,j);
    end
end

%unit conversion from cm/s^2 to microgal
dg_BP=dg_BP*(10^8);

for i=1:length(dg_BP)
    for j=1:length(dg_BP)
        g_corr(i,j)=g_corr(i,j)+dg_BP(i,j);
    end
end

g_elev=g_corr;
figure;
contourf(Xg,Yg, g_elev);
colorbar;
title('Bouguer Plate Corrected Gravity Contour Plot');
xlabel('Easting(km)');
ylabel('Northing(km)');

%Part 19
%Terrain correction
dg_terr=zeros(Nx,Ny);
rho_B=2650;
G=6.67*10^-11;

%get x and y coordinates from Xg and Yg and save them as a row vector xi
%using proper units
for j=1:Nx,
    for i=1:Ny,
        xi = [Xg(i,j)*1000, Yg(i,j)*1000,z_datum];
        for k = 1:151,
            for n=1:151,
            %save the coordinates of the local centres of mass relative to
            %the datum
            xm=[Xt(k,n)*1000, Yt(k,n)*1000, .5*(Zt(k,n)-z_datum)];
            dA= 1000*1000;
            %compute the mass of terrain above and below datum
            dm=rho_B*((Zt(k,n)-z_datum))*dA;
            %calculate the gravity effect using function from last lab 
            dg_terr(i,j)=dg_terr(i,j)+abs(grav_eff_point(xi,xm,dm,G));
            end
        end
    end
end

%changing units from m/s^2 to microgal
dg_terr= dg_terr*10^8; %to microGal

%increment dg_terr by the absolute gravity effect
g_corr=g_corr + dg_terr;
g_terr=g_corr;

%create contour plot for g_terr
figure;
contourf(Xg,Yg,g_terr);
xlabel('X [km]', 'fontweight', 'bold');
ylabel('Y [km]', 'fontweight', 'bold');
title('Terrain Correction Contour Plot');

%set colourbar limits
g_terr_min=min(min(g_terr));
g_terr_max=max(max(g_terr));

%add a colourbar
h_c=colorbar;
ylabel(h_c, 'g [microGal]', 'fontweight', 'bold');
caxis([g_terr_min, g_terr_max]);

%Part 21
dg_rgnl=mean(mean(g_corr));
for i=1:Nx
    for j=1:Ny
        g_corr(i,j)=g_corr(i,j)-dg_rgnl;
    end
end

g_anom=g_corr;

%create contour plot for g_anom
figure;
contourf(Xg,Yg,g_anom);
xlabel('X [km]', 'fontweight', 'bold');
ylabel('Y [km]', 'fontweight', 'bold');
title('Geoid Effects Corrected Contour Plot');

%set colourbar limits
g_anom_min=min(min(g_anom));
g_anom_max=max(max(g_anom));

%add a colourbar
h_c=colorbar;
ylabel(h_c, 'g [microGal]', 'fontweight', 'bold');
caxis([g_anom_min, g_anom_max]);

%Part 22
%create contour plots of dgdx and dgdy
dgdy=zeros(size(g_anom));
dgdx=zeros(size(g_anom));
%start at 2 to avoid index errors
for i=2:Nx
    for j=2:Ny
        dgdx(i,j)=((g_anom(i,j)-g_anom(i,j-1))/1000);
        dgdy(i,j)=((g_anom(i,j)-g_anom(i-1,j))/1000);
    end
end

%create a contour plot of dgdx
figure;
contourf(Xg,Yg, dgdx);
h_c = colorbar;
ylabel(h_c,'\partialg/\partialx [\muGal/m]','fontweight','bold');
title('dg/dx Contour Plot');
xlabel('Easting(km)');
ylabel('Northing [km]','fontweight','bold');

%create a contour plot of dgdy
figure;
contourf(Xg,Yg, dgdy);
h_c = colorbar;
ylabel(h_c,'\partialg/\partialy [\muGal/m]','fontweight','bold');
title('dg/dy Contour Plot');
xlabel('Easting(km)');
ylabel('Northing(km)');

%compute second derivative of dg/dz
dgdz_second=zeros(size(g_anom));
%make sure index doesn't crash
for i =2:Nx-1
    for j=2:Ny-1
        dgdz_second(i,j)=(1/(1000^4))*((2*((1000^2)+(1000^2))*g_anom(i,j))-(((1000^2)*(g_anom(i+1,j)+g_anom(i-1,j)))+((1000^2)*(g_anom(i,j+1)+g_anom(i,j-1)))));
    end
end

%plot contour of second derivative of dg/dz
figure;
contourf(Xg,Yg, dgdz_second);
h_c = colorbar;
ylabel(h_c,'[\muGal/m]','fontweight','bold');
title('Second Partical Derivative of Gravity with Respect to Depth Contour Plot');
xlabel('Easting(km)');
ylabel('Northing(km)');



