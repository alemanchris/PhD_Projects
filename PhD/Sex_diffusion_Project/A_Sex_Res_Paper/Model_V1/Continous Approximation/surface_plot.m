%% graph of the surface of the excess demands
%%
%clear all
%load('equilibrium.mat')
%%
clear all
%[X,Y] = meshgrid(linspace(0,0.04,4),linspace(0.0001,0.4,4));
[X,Y] = meshgrid(linspace(0,0.04,4),linspace(0.0001,0.8,4));
Z_ass = zeros(size(X));
Z_sex = zeros(size(X));
for i=1:size(X,1)
    for j=1:size(Y,1)
        tic
        disp([i,j])
        excess=clearing([X(i,j),Y(i,j)]);
        Z_ass(i,j)=excess(1);
        Z_sex(i,j)=excess(2);
        toc
    end
end
% Graph Excess sex
%%
figure(1)
surf(X,Y,Z_ass)
ylabel('Price')
xlabel('Rate')
title('Excess Assets') 

%
figure(2)
surf(X,Y,Z_sex)
ylabel('Price')
xlabel('Rate')
title('Excess Sex') 
%%
figure(3)
%CO(:,:,1) = zeros(4); % red
CO(:,:,1) = zeros(4); % red
CO(:,:,2) = ones(4).*linspace(0.5,0.6,4); % green
CO(:,:,3) = ones(4).*linspace(0,0.1,4); % blue
CO2(:,:,1) = ones(4).*linspace(0.4,1,4); % green
CO2(:,:,2) = zeros(4); % green
CO2(:,:,3) = zeros(4); % blue
surf(X,Y,Z_ass,CO); hold on
surf(X,Y,Z_sex,CO2)
%surf(X,Y,Z_sex,'FaceAlpha',0.5)
ylabel('Price')
xlabel('Rate')
title('General Equilibrium: Sex and Assets market')  
%}
% rate [0,0.014] 0.0122
% Price [0.24,0.4] 0.1466

%%
for k=1:3
  saveas(figure(k),sprintf('FIG_EQUILIBIUM%d.png',k)); % will create FIG1, FIG2,...
end
%%
filename = 'equilibrium.mat';
save(filename)