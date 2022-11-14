close all; clear; clc;
% This script was written by Anthony Owusu-Mensah
% current as at November 2019
% This script uses finite difference method to implement the class project(EEN582-Course project.pdf)
% EEN582-Course project.pdf is attached to this repository

lX = 3; % Length in the x-direction
lY = lX; % Lenght in y-direction
tend = 25; % Simulation duration
N = 10; % Number of Nodes
dx = lX/(N - 1); % deltaX
dy = dx; % deltaY
h = dx*dx;
dt = 0.001; % time step
t = 0:dt:tend;
Nt = length(t);
r = dt/h;

%% Initial condition
% Q(x,y,t) is a function of x,y and t
Q = -80*ones(N,N,Nt);
%% No flux Boundary Conditions
Q(:,end,:) = Q(:,end-1,:); % Top 
Q(:,1,:) = Q(:,2,:); % Bottom
Q(1,:,:) = Q(2,:,:); % Left 
Q(end,:,:) = (4/3)*Q(end-1,:,:); % Rigth 

%% Dirichlet Boundary conditions
Q(end-1:end, end-1:end,:) = 100; %Top 4 Rigth corners
Q(end-1:end, 1:2 ,:) = -100; %Bottom 4 Rigth corners
%%

%% Parameters for convergence
iter = 0;
%% err = 5;
errorTolerance  = 1e-6;
numberOfIterations = 200;
%% Phi changes in space and time

while iter < numberOfIterations

    iter = iter + 1;
    Qold = Q;
    %% Solve the interior Nodes (i,j,k are indices for x,y and t respectively)
    for k = 1:Nt-1   
        for i = 2:N-1
            for j = 2:N-1
                Q(i,j,k+1)= r*(Q(i+1,j,k)+ Q(i-1,j,k)+ Q(i,j+1,k)+ Q(i,j-1,k))...
                    -(4*r + dt-1)*Q(i,j,k);
            end
        end
    end
    
    % Test for convergence
     sumOld = sum((Qold.*Qold),'all') ; sumNew = sum((Q.*Q), 'all');   
     error_ = abs(sumNew - sumOld); 
     if(error_ <= errorTolerance)
        break;
     end

 %% Apply boundary conditions
 %% No flux Boundary Conditions (Neumann Boundary Condition)
    Q(:,end,:) = Q(:,end-1,:); % Top 
    Q(:,1,:) = Q(:,2,:); % Bottom
    Q(1,:,:) = Q(2,:,:); % Left 
    Q(end,:,:) = (4/3)*Q(end-1,:,:); % Rigth 

    %% Dirichlet boundary conditions
    %% Top Rigth corners
    Q(end-1:end, end-1:end,:) = 100;  
    %% Bottom Rigth corners
    Q(end-1:end, 1:2,:) = -100;
end

%% This section of the code is for Visualization 
figure('Color','w')
numFigs = [0,0.3333,0.5,2*dx, 3*dx,1.5,5*dx,2,3,5,10,25]; %% Timepoints of interest
[~,m] = size(numFigs);
for xx = 1: m
    ab_ = find(t >= numFigs(xx)); % The time array is of type float and will assume 0.5 != 0.500.
    idxT = ab_(1); %% Index corresponding to the time
    name = ['Time - ',num2str(numFigs(xx)),'s'];
    subplot(4,3,xx)
% %     f = figure('Name',name,'NumberTitle','off');
% %     f.Color = [1 1 1];
        X = Q(:,:, idxT);
% %         imagesc((X))
        imagesc(rot90(X)); % Rotate to match the cartesian coordinates
        colormap(jet(260))
        shading interp
        title(name);
% %         if xx == m
        c = colorbar;
        c.Location = 'eastoutside';
        c.Label.String = 'Electrical Potential (mV)';
        caxis([-100, 100])
%         end
end

