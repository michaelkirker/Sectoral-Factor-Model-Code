%% ===================================================================== %%
%%      TRADABLE & NONTRADABLE DYNAMIC FACTOR MODEL OF CORE INFLATION    %%
%% ===================================================================== %%
%
% This function estimates the a dynamic factor model using Bayesian techniques. Factors are
% identified as either tradable or nontradable factors. For further details see RBNZ discussion
% paper DP2010/13 "What drives core inflation? A dynamic factor model analysis of tradable and
% nontradable prices".
%
% An additional adjustment has been made to the model to allow for the
% effects of GST, and other one
% off tax changes that would usually bias the core inflation estimate. This is detailed in a technical
% appendix located with the model files.



function [results] = DFMcore(data,numfact,fac_loadings,reps,burn,keyindx,pbar,restriction)
%   REQUIRED INPUTS:
%
%   [DATA]: A TxNN matrix of input data for the model (where T is the number of periods, and NN is
%   the number of series in the panel). Data must be converted to a APC measure and standardised
%   before entering the model.
%
%   [NUMFACT]: A 1x2 vector. The first element of NUMFACT denoted the number of tradable factor in
%   the model. The second element of NUMFACT denotes the number of nontradable factors in the model.
%
%   [fac_loadings]: A NNx2 matrix. Denotes which series are to be loaded onto the tradable and nontradable
%   factors. For each ith element in the first (second) column, the value is 1 is the ith series in
%   the panel is loaded onto the tradable (nontradable) factor, and 0 otherwise.
%
%   [REPS]: The total number of repetitions/draws of the Gibbs sampling algorithm to be performed.
%
%   [BURN]: The number of draws to be burnt before saving the results.
%
%   [KEYINDX]: 3x1 vector identifying which colums PCPIS (CPI inflation), PTR (Tradable CPI), and
%   PNT (Non-tradable CPI).
%
%   [PBAR]: 1 = display progress bar, 0 dont display
%
%   [RESTRICTION]: 1 = apply TR/NT identification restrictions, 0 = estimate simple DFM model with
%   sum(NUMFACT) factors . Will also need to ensure all fac_loadings are 1.
%
%
%
%   OUTPUTS:
%
%   [RESULTS]: A structure containing the results from the saved draws
%


%=========================================================================%
%                   STEP 1 - INITIAL DATA SETUP ETC
%=========================================================================%

KK = sum(numfact);  % total number of factors in the DFM

[T,NN] = size(data);% T = number of periods, NN = number of series

Flags = 2; % number of lags in the factor transition equation

Ilags = 1; % number of lags in the idiosyncratic error terms

Mlags = max([Flags Ilags]); % Max lags

T = T-Mlags;


%-------------------------------------------------------------------------%
%                   1a - process factor loading fac_loadings
%-------------------------------------------------------------------------%

% Extract tradable and non-tradable series indicies
tronlyindx = fac_loadings(:,1)==1 & fac_loadings(:,2)==0;
ntonlyindx = fac_loadings(:,1)==0 & fac_loadings(:,2)==1;


% dummy matrix (NNxKK) controling which series are loaded onto each factor
dummy_mat = [repmat(fac_loadings(:,1),1,numfact(1)) repmat(fac_loadings(:,2),1,numfact(2))];

%-------------------------------------------------------------------------%
%           1b - Compute principle component model on data
%-------------------------------------------------------------------------%
if restriction == 1
    % Get principle component factors for tradable series
    pmattr = extract(data(:,tronlyindx(4:end)),numfact(1)); % PC factor using tr data (dropping first 3 aggregate series)
    
    % Get principle component factors for non-tradable series
    pmatnt = extract(data(:,ntonlyindx(4:end)),numfact(2)); % PC factor using nt data
    
    %put PC factors intro matrix form
    pmat = [pmattr pmatnt];
    
else
    % No identifying restriction
    pmat = extract(data,sum(numfact));
end

% Calculate factor loadings
F0 = dummy_mat.*(pmat(Mlags+1:end,:)\data(Mlags+1:end,:))';

% To ensure the algorithm can initiate successfully, we need the tradable and non-tradable inflation
% series to be positively loaded onto each factor. If not, we simply invert each that particular
% factor and factor loading.
if restriction == 1
    for ii = 1:numfact(1)  % for each tradable factor
        if F0(keyindx(2),ii)<0       % where the factor loading for tradable inflation is negative...
            pmat(:,ii) = (-1).*pmat(:,ii);
            F0(:,ii)   = (-1).*F0(:,ii);
        end
    end
    
    for ii = numfact(1)+1:KK   % for each nontradable factor
        if F0(keyindx(3),ii)<0           % where the factor loading for nontradable inflation is negative...
            pmat(:,ii) = (-1).*pmat(:,ii);
            F0(:,ii)   = (-1).*F0(:,ii);
        end
    end
else
    
    for ii = 1:sum(numfact)  % for each factor
        if F0(keyindx(1),ii)<0       % where the factor loading for CPI inflation is negative...
            pmat(:,ii) = (-1).*pmat(:,ii);
            F0(:,ii)   = (-1).*F0(:,ii);
        end
    end
    
    
end



verr = std(data(Mlags+1:end,:)-pmat(Mlags+1:end,:)*F0').^2;

scaler = diag(median(verr)./(pmat(Mlags+1:end,:)'*pmat(Mlags+1:end,:)));


N0=diag(scaler,0);%eye(KK)*scaler; % Prior variance of the factor loadings

pmat0=pmat;

beta0=[pmat(Mlags,:) pmat(Mlags-1,:)];  %state vector S[t-1|t-1]                                                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ns=size(beta0,2); % number of state variables
P00=eye(ns)*0.1;  %P[t-1|t-1]                                                                                                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rmat=ones(NN,1); %arbitrary starting value for the variance of the idiosyncratic component
alpha=zeros(NN,1); %arbitrary starting value for the alpha


stdscale=NaN(KK,1);

for ii = 1:KK
    Y=pmat(Mlags+1:end,ii);
    
    X = lag0multi(pmat(:,ii),Mlags);
    
    B=X(Mlags+1:end,:)\Y;
    
    err = Y-X(Mlags+1:end,:)*B;
    stdscale(ii)=std(err);
    
end;


%-------------------------------------------------------------------------%
%           1d - Pre-allocate matrices to store results in
%-------------------------------------------------------------------------%

results = struct(); % This structure will hold all the final outputs of the estimation

results.floads=NaN(NN,reps-burn,KK);                              % output/results matrix storing factor loading draws
results.NCpersistence=NaN(NN,reps-burn);                          % output/results matrix storing the estimates of alpha (AR coefficent of non-core inflation)
results.factors=NaN(T,reps-burn,KK);                              % factor estimates
results.transpara=NaN(2,reps-burn,KK);                            % transitory equation parameters: F = beta1*F{-1} + beta2*F{-2}
results.numfact=numfact;                                          % number of tradable and nontradble factors
results.reps = reps;                                              % Total number of reps in estimation
results.burn = burn;                                              % Number of burn in draws
results.varexpln = [];
results.dummymatrix=dummy_mat;                                    % dummy matrix controling the factor loadings

results.logPost = NaN(reps-burn,1);                               % Log posterior of the model
results.logL = NaN(reps-burn,1);                                  % Log likelihood of the model
results.logLprior = NaN(reps-burn,1);                             % prior likelihood of the model

% Pre-allocate matrices not saved down, but still used for calculations

outH= NaN(NN-3,ns,reps-burn);
outR = NaN(NN-3,NN-3,reps-burn);
outF = NaN(ns,ns,reps-burn);
outQ = NaN(ns,ns,reps-burn);
outlik = NaN(reps-burn,1);



if pbar
    hh   = waitbar(0,'Bayesian estimation starting...');
    set(hh,'Name','Bayesian DFM estimation progress')
end

mm=1; % running index of the draws to save

for m=1:reps
    
    
    %=========================================================================%
    %                   STEP 2 - SAMPLE FACTOR LOADINGS
    %=========================================================================%
    
    fload=NaN(NN,KK);      % stores the factor loading estimates for this iteration
    E=NaN(T+1,NN);           % Matrix to store "non-core" error terms in for this iteration
    
    for i=1:NN % for each series in the panel
        
        y=data(:,i);
        
        if m>1
            x = [pmat0(1:Mlags,:); pmat];
        else
            x=pmat;
        end
        %-------------------------------------------------------------------------%
        %                 2a - Remove serial correlation
        %-------------------------------------------------------------------------%
        
        % Rearrange model to the following form:
        %
        %    ystar = coeff'*xstar;
        %
        % Where:
        %   ystar = pi-alpha*pi(-1)
        %   xstar = [Z factor]-alpha*[Z(-1) factor(-1)]
        %   coeff = [thetaZ fload]
        
        
        % Transform regression variables to remove serial correlation
        ystar = y(Mlags+1:end,:)-y(Mlags:end-1,:)*alpha(i);
        xstar = x(Mlags+1:end,:)-x(Mlags:end-1,:)*alpha(i);
        
        %         % Truncate transformed variables to remove missing observation
        %         ystar = ystar(2:end);
        %         xstar = xstar(2:end,:);
        
        
        %-------------------------------------------------------------------------%
        %                 2b - Sample factor loadings from posterior
        %-------------------------------------------------------------------------%
        
        % factor loadings are normally distributed and can sampled from the
        % posterior distribution:
        %
        %   fload ~ N(M,V)
        %
        % M and V will be a weighted average between prior and MLE estimate
        % (based on their relative variances
        
        % Posterior mean
        M = (inv(N0)+(1/rmat(i))*(xstar'*xstar))\(N0\F0(i,:)' + (1/rmat(i))*xstar'*ystar);
        % Posterior variance
        V=inv(inv(N0)+(1/rmat(i))*(xstar'*xstar)); % equivalent to V=rmat(i)*inv(xstar'*xstar);
        
        
        if restriction == 1
            % Draw new factor loadings for series i
            if i == keyindx(2); % If series i is the tradable inflation series impose factor onto the tradble factor(s) loading is positive
                chck=-1;
                while chck<0;
                    % draw new candidate draw sample
                    floadi=dummy_mat(i,:)'.*(M+(randn(1,KK)*chol(V))');
                    % check to see if factor loading(s) is/are positive
                    if min(floadi(1:numfact(1)))>0
                        chck=1;
                    end
                end
            elseif i ==keyindx(3); % If series i is the nontradable inflation series impose factor onto the nontradble factor(s) loading is positive
                chck=-1;
                while chck<0;
                    % draw new candidate draw sample
                    floadi=dummy_mat(i,:)'.*(M+(randn(1,KK)*chol(V))');
                    % check to see if factor loading(s) is/are positive
                    if min(floadi(numfact(1)+1:end))>0
                        chck=1;
                    end
                end
            else % for all other series in the panel (i does not equal tradable or nontradble series)
                floadi=dummy_mat(i,:)'.*(M+(randn(1,KK)*chol(V))');
                % i.e. we do not impose any sign restrictions
            end
            
            
        else
            
            if i == keyindx(1);
                chck=-1;
                while chck<0
                    floadi=dummy_mat(i,:)'.*(M+(randn(1,KK)*chol(V))');
                    % check to see if factor loading(s) is/are positive
                    if min(floadi)>0
                        chck=1;
                    end
                end
            else
                floadi=dummy_mat(i,:)'.*(M+(randn(1,KK)*chol(V))');
            end
            
            
        end
        
        
        % compute serially correlated residuals for the next step
        E(:,i)=y(Mlags:end,:)-x(Mlags:end,:)*floadi;      % (TxN matrix storing all "noncore" error terms)
        
        % Store accepted factor loading draws
        fload(i,:)= floadi';
        
        
    end
    
    %=========================================================================%
    %       STEP 3 - SAMPLE THE SERIAL CORRELATION COEFFICIENT (ALPHA)
    %=========================================================================%
    
    alpha=NaN(NN,1);      % will store alpha estimates for this iteration inside this matrix
    Epsilon=NaN(T,NN);  % iid residual estimates
    
    %
    % Model:
    %
    % y = alpha*y(-1) + Epsilon
    %
    % Where: y is the E errors from step 2.
    
    for i=1:NN % for each series
        
        y=E(:,i);       % = noncore error term
        x=lag0(y,1);    % = lag of noncore error term
        
        %truncate data to account for lost observation
        y=y(2:end);
        x=x(2:end);
        
        % Estimate of Alpha can be sampled from the posterior:
        %
        % alpha_i ~ N(M,V)
        %
        % Because we have no prior on alpha, M and V are equal to the MLE
        % estimates
        
        %posterior mean
        M=(x'*x)\(x'*y);       % equivalent to M=inv(x'*x)*(x'*y);
        
        
        
        %posterior variance
        V=rmat(i)/(x'*x);      % equivalent to V=rmat(i)*inv(x'*x);
        
        %sample but ensure |alphai|<1 (i.e. non core is a
        %stationary process):
        
        chck=-1;
        while chck<0;
            alphai=M+(randn(1,1)*sqrt(V));
            if abs(alphai)<=1
                chck=1;
            end
        end
        
        %save estimate into alpha matrix
        
        alpha(i) = alphai;
        
        %save iid residuals
        Epsilon(:,i) = y-x*alphai;
    end
    
    
    %% STEP 4 SAMPLE VARIANCE R
    %
    % The variance of the white noise shock to non core inflation
    
    rmat = NaN(NN,1); % store esimtates i rmat
    
    for i=1:NN % for each series
        rmati= IG(0,0,Epsilon(:,i)); % sample variance from the inverse gamma distribution
        rmat(i) = rmati;
    end
    
    
    %% STEP 5 SAMPLE THE COEFFICIENTS OF THE TRANSITION EQUATION
    %
    % F_i,t = rho1*F_i,t-1 + rho2*F_i,t-2 + err
    
    
    rhomat=NaN(2,KK); % store transitory estimates here
    
    for k = 1:KK; % for each factor
        
        y=pmat(:,k); % equals that factors estimate
        
        x=[lag0(y,1) lag0(y,2)];
        
        % truncate data
        y=y(3:end,:);
        x=x(3:end,:);
        
        % rho coefficients can be sampled from normal distribution:
        %
        % rho ~ N(M,V)
        
        
        %posterior mean
        M = ( x'*x)\(x'*y); %M=inv(x'*x)*(x'*y);  (1/(stdscale(k)^2)).*
        
        %posterior variance
        V=inv( (1/(stdscale(k)^2)).* x'*x);  %variance is normalised to 1 for identification of scale  
        
        %sample new estimates, but only accept if coefficients imply stationarity
        chck=-1;
        while chck<0;
            rho1=M+(randn(1,2)*cholx(V))';
            rho1x=[rho1(1:2)';[1 0]]; %companion form
            
            ee=max(abs(eig(rho1x)))>1; %all eigenvalues less than unity for stationarity
            if ee==0
                chck=1;
            end
            
        end
        
        rhomat(:,k)=rho1; % store accepted coefficients
    end
    
    
    
    
    %% STEP 6 SAMPLE FACTORS
    
    %step 6A - form matrices of the state space
    
    % State space representation:
    %       Y-Z=H*factors+e
    %       factors=MU+F*Factors(-1)+v
    %       e~N(0,R)
    %       v~N(0,Q)
    
    H=zeros(NN-3,KK*2);
    H(:,1:KK)=fload(4:end,:);
    H(:,KK+1:2*KK)=fload(4:end,:).*repmat(-alpha(4:end),1,size(fload,2));
    
    %matrix R
    R=diag(rmat(4:end));
    
    % F has the following form (example is for 2 factors)
    %
    % F = [ rho1(1) 0       rho1(2) 0;
    %       0       rho2(1) 0       rho2(2);
    %       1       0       0       0;
    %       0       1       0       0];
    
    F1 = diag(rhomat(1,:));
    F2 = diag(rhomat(2,:));
    F3 = eye(KK);
    F4 = zeros(KK);
    
    F = [F1 F2;
        F3 F4];
    
    %matrix Q
    Q=zeros(size(F,1),size(F,1));
    Q(1:KK,1:KK)=eye(KK);   %normalised
    
    
    datastar=( (data(Mlags+1:end,4:end)) - (data(Mlags:end-1,4:end).*repmat(alpha(4:end)',size(data(Mlags+1:end,4:end),1),1))   );  %remove serial correlation
    %datastar(1,:)=datastar(2,:);
    
    
    %step 6b - Carter and Kohn algorithm to draw the factor
    [beta2,lik]=CarterKohn(T,ns,F,Q,R,H,datastar,beta0,P00,KK);
    
    
    
    pmat=beta2(:,1:KK);   %update the factors
    
    %% IF WE HAVE PASSED THE BURN IN PERIOD, SAVE DOWN THE RESULTS FOR THIS ITERATION
    if m>burn
        
        results.floads(:,mm,:)=fload;                                       % save factor loading estimates
        results.NCpersistence(:,mm)=alpha;                                  % save estimates of alpha (AR(1) coefficient for noncore error term
        results.factors(:,mm,:)=pmat;                                       % save factor estimates
        results.transpara(:,mm,:)=rhomat;                                   % save transatory parameter estimates (rhos)
        
        % compute variance decomposition
        outv=NaN(NN,KK+1); % variance explained will be stored in this matrix
        
        for i=1:NN
            Ei=E(2:end,i);
            A0=chol(cov([pmat Ei]));
            A01=A0./repmat(diag(A0),1,size(A0,2));
            temp=[pmat Ei]/(A01);
            
            %total variance = variance of each times factor loadings squared (factor loading for error term = 1)
            tot=(diag(cov(temp))'.*[fload(i,:).^2 1]);
            
            outv(i,:)=tot./sum(tot,2); % proportion of total variance explained by each element.
        end
        
        results.varexpln(:,mm,:)=outv;
        outH(:,:,mm)=H;
        outR(:,:,mm)=R;
        
        outF(:,:,mm)=F;
        outQ(:,:,mm)=Q;
        outlik(mm)=lik;
        
        
        
        mm=mm+1;
    end
    
    
    %% Show progress of estimation in the waitbar
    if pbar
        waitbar(m/reps,hh,sprintf('%.3f percent complete',m/reps*100));
    end
    
    
end

if pbar
    close(hh); % Close the waitbar
end

%% START OF OTHER CALCULATIONS
%%

% BDIC
Hm=squeeze(mean(outH,3));
Rm=squeeze(mean(outR,3));

Fm=squeeze(mean(outF,3));
Qm=squeeze(mean(outQ,3));

alpham=squeeze(mean(results.NCpersistence,2)); % alpham=squeeze(mean(outalpha,3));




datastarm=( (data(:,4:end)) -(lag0(data(:,4:end),1).*repmat(alpham(4:end)', size(data(:,4:end),1) ,1)));  %remove serial correlation


datastarm(1,:) = datastarm(2,:);

[~,likm]=CarterKohn(T,ns,Fm,Qm,Rm,Hm,datastarm,beta0,P00,KK);


params=mean(-2*outlik)-(-2*likm);
results.bdic = mean(-2*outlik)+params;


