function [driver, SCL, MSE] = sparsEDA(signalIn,sr,graphics,epsilon,Kmax,dmin,rho)
    % Authors, F. Hernando-Gallego, D. Luengo-Garcia and A. Artes-Rodriguez
    % Dpt. of Signal Theory and Communications
    % Universidad Carlos III de Madrid
    % ---------------------------------------------------------------------
    % Abstract. This method provide the low component and the driver
    % response of a galvanic skin response signal longer than 70 seconds.
    % ---------------------------------------------------------------------
    % signal  [vector]: galvanic skin response
    % sr         [int]: sample rate
    % graphics   [0/1]: show graphics (1) or not (0)
    % epsilon [double]: step remainder
    % maxIters   [int]: maximum number of LARS iterations
    % dmin    [double]: maximum distance between sparse reactions
    % th      [double]: minimun threshold of sparse reactions 
    %
    % OUTPUTS
    % driver  [vector]: driver responses, tonico componet
    % SCL     [vector]: low component
    % MSE     [vector]: reminder of the signal fitting

    
    %% 0. Exceptions
    if (length(signalIn)/sr < 80)
        display('Signal not enought large. May be longer than 80 senconds')
    end
    
    if (sum(isnan(signalIn))>0)
        display('Signal contains NaN')
    end
    %% 1.PREPROCESSING

    
    % addition start and end to process the signal completly
    signalAdd = ...
        [signalIn(1)*ones(20*sr,1); signalIn; signalIn(end)*ones(60*sr,1)];
    
    if (sr>8)
        [P,Q] = rat(8/sr);
        signalAdd = resample(signalAdd,P,Q);
        Nss = floor(8*length(signalIn)/sr);
        sr = 8;
    else
        Nss = length(signalIn);
    end
    
    Ns  = length(signalAdd);
    b0 = 0;
     
    pointerS = (20*sr) + 1;
    pointerE = pointerS + Nss;
    signalRs = signalAdd(pointerS:pointerE);
    
    %% Overlap save
    % Loop Analysis
    durationR = 70;  % Slots 70 seconds (60s + 10s for deconvolution)
    Lreg = 20*sr*3;  % 60 seconds
    L = 10*sr;       % pulse
    N = durationR*sr;% 70 seconds
    T = 6;           % CSL additional columns
    
    % Toeplitz Matrix
    % SCR
    Rzeros = zeros(N+L,Lreg);
    srF = sr*[0.5 0.75 1 1.25 1.5]; % Multianalysis
    for j=1:length(srF)
        t_rf = 0:(1/srF(j)):10; % 10 seconds
        taus = [0.5, 2, 1];
        rf_biexp = exp(-t_rf/taus(2)) - exp(-t_rf/taus(1));   
        rf_est = taus(3)*rf_biexp;
        rf_est = rf_est/sqrt((sum(rf_est.^2)));
        Rzeros(1:N,1 + (j-1)*Lreg:j*Lreg) =...
            toeplitz([rf_est zeros(1,N-length(rf_est))],zeros(1,Lreg));
    end
    R = Rzeros(1:N,1:5*Lreg);
    % SCL
    R = [zeros(N,T) R];
    R(1:Lreg,1) = (+1)*linspace(0,1,Lreg)';
    R(1:Lreg,2) = (-1)*linspace(0,1,Lreg)';
    R((Lreg/3)+1:Lreg,3) = (+1)*linspace(0,2/3,(2*Lreg)/3)';
    R((Lreg/3)+1:Lreg,4) = (-1)*linspace(0,2/3,(2*Lreg)/3)';
    R((2*Lreg/3)+1:Lreg,5) = (+1)*linspace(0,1/3,Lreg/3)';
    R((2*Lreg/3)+1:Lreg,6) = (-1)*linspace(0,1/3,Lreg/3)';    
    Cte = sum(R(:,1).^2);
    R(:,1:6) = R(:,1:6)/sqrt(Cte);

    %% Loop
    cutS = 1;
    cutE = N;
    slcAux   = zeros(1,Ns);
    drverAux = zeros(1,Ns);
    resAux   = zeros(1,Ns);
    aux      = 0;

    while (cutE<Ns)
        aux = aux + 1;
        signalCut = signalAdd(cutS:cutE);
        
        if(b0==0)
            b0 = signalCut(1);
        end
        
        signalCutIn = signalCut - b0;
        [beta,~, activationHist,~,~,~] = lasso(R,signalCutIn,sr,Kmax,epsilon);
        % display(activationHist)
        signalEst = R*beta + b0;
        
        remAout = (signalCut - signalEst).^2;
        res2 = sum(remAout(20*sr+1:(40*sr)));
        res3 = sum(remAout(40*sr+1:(60*sr)));
        
        jump = 1;
        if (res2 < 1)
            jump = 2;
            if (res3 < 1)
                jump = 3;
            end
        end     
        
        SCL = R(:,1:6)*beta(1:6) + b0;
        SCRline = beta(7:end);
        SCRaux = zeros(Lreg,5);
        SCRaux(:) = SCRline;   
        driver = sum(SCRaux');
        
        b0 = R(jump*20*sr,1:6)*beta(1:6) + b0;
        
        drverAux(cutS:cutS+(jump*20*sr)-1) = driver(1:jump*sr*20)';
        slcAux(cutS:cutS+(jump*20*sr)-1) = SCL(1:jump*sr*20)';
        resAux(cutS:cutS+(jump*20*sr)-1) = remAout(1:jump*sr*20)';
        cutS = cutS + (jump)*20*sr;
        cutE = cutS + N -1;
    end

    SCRaux    = drverAux(pointerS:pointerE);
    SCL       = slcAux(pointerS:pointerE);
    MSE       = resAux(pointerS:pointerE);
    
    %% POST-PROCESSING: REMOVING UNNECESSARY SCR COMPONENTS

    
    ind = find(SCRaux>0);
    scr_temp = SCRaux(ind);
    [scr_ord, ind2] = sort(scr_temp, 'descend');
    scr_fin = scr_ord(1);
    ind_fin = ind(ind2(1));
    for i = 2:length(ind2)
        if all(abs(ind(ind2(i))-ind_fin)>=dmin)
            scr_fin = [scr_fin scr_ord(i)];
            ind_fin = [ind_fin ind(ind2(i))];
        end
    end
    driver = zeros(size(SCRaux));
    driver(ind_fin) = scr_fin;

    scr_max = scr_fin(1);
    threshold = rho*scr_max;
    driver(driver < threshold) = 0;
    
        %%
    if (graphics)
        f1 = figure('Position',[-1900 8 1890 990]);
        t = 0:(1/sr):(length(signalRs)/sr)-(1/sr);
        subplot(311)
        plot(t,signalRs,'k','Linewidth',1.5)
        legend('signal','Location','NorthEastOutside')
        subplot(312)
        plot(t,signalRs,'k','Linewidth',1.5)
        hold on
        plot(t,SCL)
        legend('signal','SCL','Location','NorthEastOutside')
        subplot(313)
        plot(t,driver,'-o','Linewidth',1.5)
        legend('SCR','Location','NorthEastOutside')
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [beta, numIters, activationHist, duals, lambda, res] = lasso(R, s, sr,maxIters,epsilon)

    N = length(s);
    W = size(R,2);

    OptTol = -10;
    solFreq = 0;
    resStop2 = .0005;
    lambdaStop = 0;

    % Global variables for linsolve function
    global opts opts_tr zeroTol
    opts.UT = true; 
    opts_tr.UT = true; opts_tr.TRANSA = true;
    zeroTol = 1e-5;

    x = zeros(W,1);
    x_old = zeros(W,1);
    iter = 0;

    % Correlation
    c = R'*s;

    % Maximun correlation value
    lambda = max(c);

    if (lambda < 0)
        error('y is not expressible as a non-negative linear combination of the columns of X');
    end

    % argmax de correlacion
    newIndices = find(abs(c-lambda) < zeroTol)';    

    collinearIndices = [];
    beta = [];
    duals = [];
    res = s; % res en principio la señal

    % Check stopping conditions 
    if ((lambdaStop > 0) & (lambda < lambdaStop)) | ((epsilon > 0) & (norm(res) < epsilon))
        activationHist = [];
        numIters = 0;
        % return;
    end

    % Initialize Cholesky factor of A_I 
    R_I = [];
    activeSet = [];
    for j = 1:length(newIndices)
        iter = iter+1;
        [R_I, flag] = updateChol(R_I, N, W, R, 1, activeSet, newIndices(j));
        activeSet = [activeSet newIndices(j)];
    end
    activationHist = activeSet;


    %%%%%% LOOP %%%%%%
    done = 0;
    while  ~done

        if (length(activationHist) == 4)
            lambda = max(c);
            newIndices = find(abs(c-lambda) < zeroTol)';    
            % R_I = [];
            activeSet = [];
            for j = 1:length(newIndices)
                iter = iter+1;
                [R_I, flag] = updateChol(R_I, N, W, R, 1, activeSet, newIndices(j));
                activeSet = [activeSet newIndices(j)];
            end
            activationHist = [activationHist activeSet];
        else
            lambda = c(activeSet(1));
        end

        % Compute Lars direction - Equiangular vector
        dx = zeros(W,1);
        % Solve the equation (A_I'*A_I)dx_I = sgn(corr_I)
        z = linsolve(R_I,sign(c(activeSet)),opts_tr);
        dx(activeSet) = linsolve(R_I,z,opts);

        v = R(:,activeSet)*dx(activeSet);
        ATv = R'*v;

        % For  Lasso, Find first active vector to violate sign constraint
        gammaI = Inf;
        removeIndices = [];

        % Find first inacti ve vector to enter the active set
        inactiveSet = 1:W;
        inactiveSet(activeSet) = 0;
        inactiveSet(collinearIndices) = 0;
        inactiveSet = find(inactiveSet > 0);

        if (length(inactiveSet) == 0)
            gammaIc = 1;
            newIndices = [];
        else
            epsilon = 1e-12; 
            gammaArr = (lambda - c(inactiveSet))./(1 - ATv(inactiveSet) + epsilon);
            gammaArr(gammaArr < zeroTol) = Inf;

            [gammaIc, Imin] = min(gammaArr);
            newIndices = inactiveSet(find(abs(gammaArr - gammaIc) < zeroTol));
         end


        gammaMin = min(gammaIc,gammaI);

        % Compute the next Lars step
        x = x + gammaMin*dx;
        res = res - gammaMin*v;
        c = c - gammaMin*ATv;

        % Check stopping condition
        if ((lambda - gammaMin) < OptTol) | ((lambdaStop > 0) & (lambda <= lambdaStop)) | ((epsilon > 0) & (norm(res) <= epsilon))
            newIndices = [];
            removeIndices = [];
            done = 1;

            if ((lambda - gammaMin) < OptTol)
                display(num2str((lambda - gammaMin)))
            end
            if ((lambdaStop > 0) & (lambda <= lambdaStop))
            end
            if ((epsilon > 0) & (norm(res) <= epsilon))
            end

        end

        if (norm(res(1:sr*20)) <= resStop2)
            done = 1;
            if (norm(res(sr*20:sr*40)) <= resStop2)
                done = 1;
                if (norm(res(sr*40:sr*60)) <= resStop2)
                    done = 1;
                end
            end
        end

        % Add new indices to active set
        if (gammaIc <= gammaI) && (length(newIndices) > 0)
            for j = 1:length(newIndices)
                iter = iter+1;
                % Update the Cholesky factorization of A_I
                [R_I, flag] = updateChol(R_I, N, W, R, 1, activeSet, newIndices(j));
                % Check for collinearity
                if (flag)
                    collinearIndices = [collinearIndices newIndices(j)];
                else
                    activeSet = [activeSet newIndices(j)];
                    activationHist = [activationHist newIndices(j)];
                end
            end
        end

        % Remove violating indices from active set
        if (gammaI <= gammaIc)
            for j = 1:length(removeIndices)
                iter = iter+1;
                col = find(activeSet == removeIndices(j));

                % Downdate the Cholesky factorization of A_I
                R_I = downdateChol(R_I,col);
                activeSet = [activeSet(1:col-1), activeSet(col+1:length(activeSet))];

                % Reset collinear set
                collinearIndices = [];
            end

            x(removeIndices) = 0;  % To avoid numerical errors
            activationHist = [activationHist -removeIndices];
        end

        if iter >= maxIters
            done = 1;
        end

        if (find(x<0)>0)
            x = x_old;
            done = 1;
        else
            x_old = x;
        end

        if done | ((solFreq > 0) & (~mod(iter,solFreq)))
            beta = [beta x];
            duals = [duals v];
        end
end

numIters = iter;
clear opts opts_tr zeroTol
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [R_I, flag] = updateChol(R_I, n, N, R, explicitA, activeSet, newIndex)
% updateChol: Updates the Cholesky factor R of the matrix 
% X(:,activeSet)'*X(:,activeSet) by adding X(:,newIndex)
% If the candidate column is in the span of the existing 
% active set, R is not updated, and flag is set to 1.

global opts_tr zeroTol
flag = 0;


newVec = R(:,newIndex);

if (length(activeSet) == 0)
    R_I = sqrt(sum(newVec.^2));
else
    if (explicitA)       
        p = linsolve(R_I,R(:,activeSet)'*R(:,newIndex),opts_tr);
    else
        AnewVec = feval(R,2,n,length(activeSet),newVec,activeSet,N);
        p = linsolve(R_I,AnewVec,opts_tr);
    end
    q = sum(newVec.^2) - sum(p.^2);
    if (q <= zeroTol) % Collinear vector
        flag = 1;
    else
        R_I = [R_I p; zeros(1, size(R_I,2)) sqrt(q)];
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R = downdateChol(R, j)
% downdateChol: `Downdates' the cholesky factor R by removing the 
% column indexed by j.

% Remove the j-th column
R(:,j) = [];
[m,n] = size(R);

% R now has nonzeros below the diagonal in columns j through n.
% We use plane rotations to zero the 'violating nonzeros'.
for k = j:n
    p = k:k+1;
    [G,R(p,k)] = planerot(R(p,k));
    if k < n
        R(p,k+1:n) = G*R(p,k+1:n);
    end
end

% Remove last row of zeros from R
R = R(1:n,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
