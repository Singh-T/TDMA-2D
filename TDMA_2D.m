% 2 D TDMA %
% Heat transfer in a 2D plate with thickness 1cm%
% Top boundary maintainsed at 100 degree Celcius.
% Constant heat flux of 500 kW/m^2 applied at west boundary
% other two boundaries are insulated. See ﻿Versteeg and Malalasekara
% Example 7.2 for problem discription.


%%%% Variables definition %%%%%
L=0.3;
nx=6; % number of points in x;
B=0.4;
ny=8; % number of points in y;
k=1000;
q=500;

%%%%%
%Considering equal dx,dy and cell area%%
dx=L/nx;
dy=B/ny;
A=0.01*0.1;
a=(k*A)/dx;

aw=zeros(nx,ny)+a;
ae=zeros(nx,ny)+a;
an=zeros(nx,ny)+a;
as=zeros(nx,ny)+a;
ap=zeros(nx,ny);
Sp=zeros(nx,ny);
Su=zeros(nx,ny);

T=zeros(nx,ny);

x = linspace(0,0.3,nx);
y = linspace(0,0.4,ny);
[X,Y] = meshgrid(x,y);

%   Discretized equation : aPTP = aWTW + aETE + aSTS + aNTN
%   Interior points : ap=aw+ae+as+an

for i = 2:nx-1
    for j = 2:nx
        ap(i,j)=aw(i,j)+ae(i,j)+as(i,j)+an(i,j);
    end
end

%Boundaries: ap=aw+ae+as+an-Sp

for i = 1:nx
    if i==1 %left boundary
        for j = 1:ny
            if j==1
                aw(i,j)=0; as(i,j)=0;
                Su(i,j)=500;
                Sp(i,j)=0;
                
            elseif j==ny 
                aw(i,j)=0; an(i,j)=0;
                Su(i,j)=2500;
                Sp(i,j)=-20;
                
            else 
                aw(i,j)=0;
                Su(i,j)=500;
                Sp(i,j)=0;
            end
        end
    elseif i==nx %right boundary
        for j = 1:ny
            if j==1
                ae(i,j)=0; as(i,j)=0;
                Su(i,j)=0;
                Sp(i,j)=0;
                
            elseif j==ny 
                an(i,j)=0; ae(i,j)=0;
                Su(i,j)=2000;
                Sp(i,j)=-20;
                
            else 
                Su(i,j)=0; ae(i,j)=0;
                Sp(i,j)=0;
            end
        end
    else
        for j = 1:ny
            if j==1
                as(i,j)=0;
                Su(i,j)=0;
                Sp(i,j)=0;
                
            elseif j==ny  
                an(i,j)=0;
                Su(i,j)=2000;
                Sp(i,j)=-20;
            end
        end
    end
    for j=1:ny
        ap(i,j)=aw(i,j)+ae(i,j)+as(i,j)+an(i,j)-Sp(i,j);
    end
end

%% TDMA:

alpha=zeros(nx,ny);
beta=zeros(nx,ny);
D=zeros(nx,ny);
C=zeros(nx,ny);

A=zeros(nx,ny);
Cdash=zeros(nx,ny);

run=1;
conv_residual=0.1; % set desired convergence value
max_res=10;
max_iter=500;

while(run<max_iter && max_res>(conv_residual))
    for i = 1:nx
         for j = 1:ny
            alpha(i,j)=an(i,j);
            beta(i,j)=as(i,j); 
            D(i,j)=ap(i,j);
            if i==1
                C(i,j)=((ae(i,j)*T(i+1,j))+Su(i,j));
            elseif i==nx
                C(i,j)=((aw(i,j)*T(i-1,j))+Su(i,j)); 
            else
                C(i,j)=((aw(i,j)*T(i-1,j))+(ae(i,j)*T(i+1,j))+Su(i,j));
            end
         end
        
         for j=1:ny
            if j==1 
                A(i,j) = (alpha(i,j)/(D(i,j)));
                Cdash(i,j) = ((0+ C(i,j))/(D(i,j)));
            else
                A(i,j) = (alpha(i,j)/(D(i,j)-(beta(i,j)*A(i,j-1))));
                Cdash(i,j) = (((beta(i,j)*Cdash(i,j-1))+ C(i,j))/(D(i,j)-(beta(i,j)*A(i,j-1))));
            end
         end
         t=T;
         for j=ny:-1:1
             if j==ny
                 T(i,j)= (0+Cdash(i,j));
             else
                 T(i,j)= ((A(i,j)*T(i,j+1))+Cdash(i,j));
             end
         end
    end
    run=run+1;
    res=abs(t-T);
    max_res=max(res, [], 'all');
    
    fprintf('Iteration.... %i', run);
    fprintf('\n');
    pause(0.05);                        % Pauses computation for contour update. Value in seconds.
    contourf(X,Y,T','ShowText','on')
    colorbar; title('Temperature distribution (K)','FontSize',22);
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
end

if max_res<(conv_residual)
    fprintf('Solution converged');
else
    fprintf('Max. iterations reached');
end
