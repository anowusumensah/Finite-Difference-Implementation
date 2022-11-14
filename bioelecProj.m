close all; clear; clc;

lX = 3; % Length in the x-direction
lY = lX; % Lenght in y-direction
tend = 25; % Simulation end time
N = 10; % Number of Nodes
dx = lX/(N - 1); % deltaX
dy = dx; % deltaY
h = dx*dx;
dt = 0.001; % time step
t = 0:dt:tend;
Nt = length(t);
r = dt/h;
%% Initial condition
Q = -80*ones(N,N,Nt);

%% Parameters for convergence
iter = 0;
err = 5;
erroTolerance  = 1e-6;
numberOfIterations = 200;
%% Phi changes in space and time

%% Boundary Conditions
Q(:,end,:) = Q(:,end-1,:); % Top 
Q(:,1,:) = Q(:,2,:); % Bottom
Q(1,:,:) = Q(2,:,:); % Left 
Q(end,:,:) = (4/3)*Q(end-1,:,:); % Rigth 

%% Top Rigth corners
Q(end-1:end, end-1:end,:) = 100;  
%% Bottom Rigth corners
Q(end-1:end, 1:2 ,:) = -100;
%%
% % while err > tol
while iter < numberOfIterations
    iter = iter + 1;
    Qold = Q;
    %% Solve the interior Nodes
    for k = 1:Nt-1   
        for i = 2:N-1
            for j = 2:N-1
                Q(i,j,k+1)= r*(Q(i+1,j,k)+ Q(i-1,j,k)+ Q(i,j+1,k)+ Q(i,j-1,k))...
                    -(4*r + dt-1)*Q(i,j,k);
            end
        end
    end

     sumold = sum((Qold.*Qold),'all') ; aa = sum((Q.*Q), 'all');   
     err = abs(aa-sumold); 
     if(err <= erroTolerance)
      break;
     end

 %% Apply boundary conditions
 %% No Boundary Conditions
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
    ab_ = find(t >= numFigs(xx)); %% Index corresponding to the time
    idxT = ab_(1);
    name = ['Time - ',num2str(numFigs(xx)),'s'];
    subplot(4,3,xx)
% %     f = figure('Name',name,'NumberTitle','off');
% %     f.Color = [1 1 1];
        X = Q(:,:, idxT);
% %         imagesc((X))
        imagesc(rot90(X)); % So it matches the number line
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

