% Huabin Zhang - huabinz@connect.hku.hk, 20250806
% QUASS process for zspec
%   ref: Sun Zhe-2011. MRM 86(2): 765-776.

function [zspecQUASSList, R1rhoList] = QUASSprocess(zspecList, offs_ppm, SATpara, R1w_vec, db0map_vec, displayflag)
% INPUT: 
%   zspecList: [nf,nzspec]
%   offs_ppm: [1, nf]
%   SATpara: [3, nzspec] or [3, 1], for each zspec: [B1(uT), Ts(s), Td(s)]
%   R1w_vec: [1, nzspec], in Hz
%   db0map_vec: [1, nzspec], in ppm

    B0 = 3; % main field

    zspecListReform = reshape(zspecList,size(zspecList,1),[]);
    [nf, nzspec] = size(zspecListReform);

    if ~exist('db0map_vec','var') || isempty(db0map_vec)
        db0map_vec = 0;
    end
    if ~exist('displayflag','var') || isempty(displayflag)
        displayflag = 0;
    end
	
    if size(SATpara,2) == 1
        SATpara = repmat(SATpara,1,nzspec);
    end
    if isscalar(R1w_vec)
        R1w_vec = repmat(R1w_vec,1,nzspec);
    end
    if isscalar(db0map_vec)
        db0map_vec = repmat(db0map_vec,1,nzspec);
    end

    R1rhoList = zeros(nf, nzspec);
    zspecQUASSList = zeros(nf, nzspec); % zspec after QUASS correction
    multistart_N = 1;

    tic
    fprintf('\t processing QUASS, ');
	backNum = 0;
    for idxZ = 1:nzspec
        if displayflag == 1
            fprintf(repmat('\b',1,backNum));
            backNum = fprintf(' %d/%d',idxZ,nzspec);
        end
        zspec = zspecListReform(:,idxZ);
        R1w = R1w_vec(idxZ);
        B1 = SATpara(1,idxZ);
        Ts = SATpara(2,idxZ);
        Td = SATpara(3,idxZ);
        
        % R1rho fitting boundaries (Hz)
        R1rho_iv = 3;
        R1rho_lb = 0.01;
        R1rho_ub = 30;
        
        for idxoffs = 1:nf
            % QUASS model, steady-state: M_sat / M_0
            theta_fit = @(offs_ppm) atan(B1./B0./( offs_ppm-db0map_vec(idxZ) )  );
            zspecModel_QUASS = @(R1rho, offs_ppm) ...
                (...
                    (1-exp(-R1w*Td)) * exp(-R1rho*Ts) + R1w/R1rho*cos(theta_fit(offs_ppm)).^2*(1-exp(-R1rho*Ts))...
                    ) ./ (1-exp(-R1w*(Ts+Td)));

            OPTIONS = optimoptions('lsqcurvefit','Algorithm','Levenberg-Marquardt','Display','off');      
            problem = createOptimProblem('lsqcurvefit','x0',R1rho_iv,'objective',zspecModel_QUASS,...
        'lb',R1rho_lb,'ub',R1rho_ub,'xdata',offs_ppm(idxoffs),'ydata',zspec(idxoffs),'options',OPTIONS);
            ms = MultiStart('Display','off');
            [par,~] = run(ms,problem,multistart_N);
            R1rhoList(idxoffs, idxZ) = par(1);
        end
        
        theta_List = atan(B1./B0./ (offs_ppm-db0map_vec(idxZ) ) ); 
        zspecQUASSList(:, idxZ) = (R1w./reshape(R1rhoList(:, idxZ),[],1)) .* cos(theta_List).^2;
    end
    fprintf(repmat('\b',1,backNum));

    zspecQUASSList = reshape(zspecQUASSList,size(zspecList));
    R1rhoList = reshape(R1rhoList,size(zspecList));
    fprintf('elapsed time: %.4f s\n', toc);
end