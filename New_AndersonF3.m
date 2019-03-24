function [YY,Z] = New_AndersonF3(ss,n,c,influence,P,trials)

%
% Pre:
%   s:              vector containing values for desired fractional 
%                   exponents
%   n:              number of time steps (n-2 is actual number of steps)
%   c:              c-value for random parameter interval length
%   influence:      how many neighbor interactions to allow (note that if
%                   this parameter is left empty the function will choose
%                   100 by default)
%   P:              chooses the particular delta vector (options are 1-6)
%                   see below for more details on the choices
%   trials:         number of trials to perform for averaging results
%
% Post:             
%   YY:             Stores the average (over all trials) of each s-value
%                   in a particular column (increasing left to right)
%   Z:              Stores each indiviual result for deeper analysis (it
%                   will store the s-values for each trial from left to 
%                   right in a multi-dimensional array)
%
% Note that this functions saves the data in particular directories
%

if nargin <= 3
    influence = 200;
end

if nargin <=4
    P = 2; % Chooses a default vector
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                      %%
%% Code for one-dimensional Anderson localization simulations via the   %%
%% spectral approach.                                                   %%
%%                                                                      %%
%% This code allows for the consideration of fractional diffusion for   %%
%% all fractions 0 < s < 2. The code employs different methods based on %%
%% the choise of s.                                                     %%
%%                                                                      %%
%% All methods herein are based on the Bessel-type approach.            %%
%%                                                                      %%
%% All necessary parameters that may be changed are outlined below.     %%
%%                                                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start of Code                                                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clc
rng('Shuffle')

max = length(ss); % Determines the number of s-values to run
YY = zeros(n-2,max); % Empty vector for storing final result
Z = zeros(n-2,max,trials); % Empty vector for storing individual results

%% The following "parfor" loop while parallelize the process of
%% running multiple trials

parfor iii = 1:trials
    
    %% Initialize relevant vector and values
    
    N = 2*influence*n+1; % Lattice size
    Y = zeros(n-2,max);  % Solution vector for given trial
    Omega = c.*(rand(N,1)-0.5*ones(N,1)); % Vector of random values
    NN = 0.5*(N-1)+1; % Midpoint of lattice (for starting point)
    
    %% The following loops runs computations for each s-value
    
    for ll = 1:max
    
        s = ss(ll); % Extracts s-value from appropriate position
    
        if s ~= 1 % Checks s-value
    
            A = zeros(N,3); % Placeholder vector for the simulation
            A(NN,1) = 1; % Builds delta_0
    
            % Weights for computing action of fractional Laplacian
            
            ks = (4^s)*gamma(0.5+s)/((sqrt(pi))*abs(gamma(-s))); % Constant
            K_s = zeros(influence-1,1); % Vector to store computed weights
            kg = gamma(1-s)/gamma(2+s); 
            for j=1:(influence-1)
                K_s(j) = ks*kg; 
                kg = ((j-s)/(j+1-s))*kg; 
            end
            
            p = 0; % Initializes distance calculation
            
            % Loop to compute first action
    
            for k=NN-influence:NN+influence
                
                f1 = 0; % Placeholder
                
                for j=1:(influence-1)
                    
                    % Primary action of fractional Laplacian at node
                    
                    f1 = f1 + (2*A(k,1)-A(k-j,1)-A(k+j,1))*K_s(j);
                    
                end     
                
                A(k,2) = f1; % Stores final value
                
            end
            
            A(:,2) = A(:,2)./norm(A(:,2)); % Normalizes vector
    
            % Loop to compute remaining actions
            
            for i=2:(n-1)
                
                % Inner loop computes each new "action"
                
                for k=NN-influence*i:NN+influence*i
                    
                f1 = 0; % Placeholder
                
                    for j=1:(influence-1)
                        
                        % Primary action of fractional Laplacian at node
                        
                        f1 = f1 + (2*A(k,2)-A(k-j,2)-A(k+j,2))*K_s(j);
                        
                    end  
                    
                    A(k,3) = f1;
                    
                end
                
                A(:,3) = A(:,3) + A(:,2).*Omega; % Adds randomness
                
                % The following two lines perform Gram-Schmidt
        
                A(:,3) = A(:,3) - ((A(:,3))'*A(:,2)).*A(:,2) - ((A(:,3))'*A(:,1)).*A(:,1);
                A(:,3) = A(:,3)./norm(A(:,3));
                
                % The next set of code computes the distance value based
                % on the P-value submitted
                
                if P == 1
                    
                    p = p + (A(NN-1,3))^2;
                    
                elseif P == 2
                    
                    p = p + 0.25*((A(NN-1,3)-A(NN+1,3)+A(NN-2,3)-A(NN+2,3))^2);
                    
                elseif P == 3
                    
                    pp = 0;
                    for ii = 1:influence-1
                        pp = pp + ((-1)^(ii-1))*(A(NN-ii,3)) + ((-1)^(ii))*(A(NN+ii,3));
                    end
                    p = p + 0.5*(1/((influence-1)))*(pp^2);
                    
                elseif P == 4
                    
                    pp = 0;
                    for ii = 1:10
                        pp = pp + ((-1)^(ii-1))*(A(NN-ii,3)) + ((-1)^(ii))*(A(NN+ii,3));
                    end
                    p = p + 0.5*(1/10)*(pp^2);
                    
                elseif P == 5
                    
                    pp = 0;
                    for ii = 1:50
                        pp = pp + ((-1)^(ii-1))*(A(NN-ii,3)) + ((-1)^(ii))*(A(NN+ii,3));
                    end
                    p = p + 0.5*(1/50)*(pp^2);
                    
                elseif P == 6
                    
                    pp = 0;
                    for ii = 1:100
                        pp = pp + ((-1)^(ii-1))*(A(NN-ii,3)) + ((-1)^(ii))*(A(NN+ii,3));
                    end
                    p = p + 0.5*(1/100)*(pp^2);
                    
                end
                
                % Stops computations when "localization" occurs
                
                if p >= 1
                    disp(sprintf('Numerical localization has occurred after %1.0f iterations. Possibly vary the influence.',i))
                    break
                end
        
                Y(i-1,ll) = sqrt(1-p); % Stores computed distance

                A(:,1) = A(:,2); % Updates vector
                A(:,2) = A(:,3); % Updates vector
                
            end
    
        end
        
        % Due to the sensitivity of this code, we will approximate the s=1
        % case via a simple averaging
        
        if s == 1
            
            p = 0;
            
           %% Generation of second vector under Laplacian
            A(NN-1,2) = -1/sqrt(2);
            A(NN+1,2) = -1/sqrt(2);
            
           %% The following loop computes the remaining iterations under the 
           %% action of the random Schrodinger operator (standard Laplacian)
    
            for i = 2:(n-1)
                for k = NN-i:NN+i
                    A(k,3) = -A(k-1,2) + 2*A(k,2) - A(k+1,2); %Laplacian action
                end
    
                A(:,3) = A(:,3) + A(:,2).*Omega; %Addition of random component
    
               %% The following performs the orthogonalization procedure and 
               %% normalizes the resulting vector
        
                A(:,3) = A(:,3) - ((A(:,3))'*A(:,2)).*A(:,2) - ((A(:,3))'*A(:,1)).*A(:,1);
                A(:,3) = A(:,3)./norm(A(:,3));
                
                % The next set of code computes the distance value based
                % on the P-value submitted
                
                if P == 1
                    
                    p = p + (A(NN-1,3))^2;
                    
                elseif P == 2
                    
                    p = p + 0.25*((A(NN-1,3)-A(NN+1,3)+A(NN-2,3)-A(NN+2,3))^2);
                    
                elseif P == 3
                    
                    pp = 0;
                        for ii = 1:influence-1
                            pp = pp + ((-1)^(ii-1))*(A(NN-ii,3)) + ((-1)^(ii))*(A(NN+ii,3));
                        end
                    p = p + 0.5*(1/((influence-1)))*(pp^2);
                    
                elseif P == 4
                    
                    pp = 0;
                    for ii = 1:10
                        pp = pp + ((-1)^(ii-1))*(A(NN-ii,3)) + ((-1)^(ii))*(A(NN+ii,3));
                    end
                    p = p + 0.5*(1/10)*(pp^2);
                    
                elseif P == 5
                    
                    pp = 0;
                    for ii = 1:50
                        pp = pp + ((-1)^(ii-1))*(A(NN-ii,3)) + ((-1)^(ii))*(A(NN+ii,3));
                    end
                    p = p + 0.5*(1/50)*(pp^2);
                    
                elseif P == 6
                    
                    pp = 0;
                    for ii = 1:100
                        pp = pp + ((-1)^(ii-1))*(A(NN-ii,3)) + ((-1)^(ii))*(A(NN+ii,3));
                    end
                    p = p + 0.5*(1/100)*(pp^2);
                    
                end
                
                % Stops computations when "localization" occurs
                
                if p >= 1
                    disp(sprintf('Numerical localization has occurred after %1.0f iterations. Possibly vary the influence.',i))
                    break
                end
        
                Y(i-1,ll) = sqrt(1-p); % Stores computed distance

                A(:,1) = A(:,2); % Updates vector
                A(:,2) = A(:,3); % Updates vector
                    
            end
            
        end
    
    end
    
    Z(:,:,iii) = Y; % Stores all of the data
    
    YY = YY + Y; % Aids in storing averages
    
end

YY = (1/trials)*YY; % Stores the average of the runs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save Data                                                            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Creates string for c value with appropriate length

s = ss;

cc = floor(c);
cString1 = num2str(cc);
lc1 = length(cString1);
ccc = c - cc;
cString2 = num2str(ccc);
lc2 = length(cString2);
cString = [cString1 '_' cString2(lc1+2:lc2)];

%% Making Directories

if exist('Fractional_Anderson_Data', 'dir') == 0 
    mkdir('Fractional_Anderson_Data');
else
    warning('off','MATLAB:MKDIR:DirectoryExists');
end
if exist(['Fractional_Anderson_Data_s_' num2str(s(1)) '-' num2str(s(max))], 'dir') == 0 
    mkdir('Fractional_Anderson_Data',['Fractional_Anderson_Data_s_' num2str(s(1)) '-' num2str(s(max))]);
else
    warning('off','MATLAB:MKDIR:DirectoryExists');
end
ddir = fullfile('Fractional_Anderson_Data',['Fractional_Anderson_Data_s_' num2str(s(1)) '-' num2str(s(max))]);
if exist(['Fractional_Anderson_Data_s_' num2str(s(1)) '-' num2str(s(max)) '_c_' cString], 'dir') == 0
    mkdir(ddir,['Fractional_Anderson_Data_s_' num2str(s(1)) '-' num2str(s(max)) '_c_' cString]);
else
    warning('off','MATLAB:MKDIR:DirectoryExists');
end

%% Saving Data

FolderDestination = ['Fractional_Anderson_Data_s_' num2str(s(1)) '-' num2str(s(max)) '_c_' cString];
data = ['Frac_Data_s_' num2str(s(1)) '-' num2str(s(max)) '_c_' cString '_n_' num2str(n) '_infl_' num2str(influence) '_trials_' num2str(trials) '_' datestr(datetime('now')) '.mat'];
matfile = fullfile(ddir,FolderDestination,data);

save(matfile);
