function [rBest,azBest,azPerf,elPerf,hist] = ...
    optRandConfig(egain,nTrials,lamEval,azTrue,elTrue,abcdBlock,...
    defdB,azGrid,elGrid,truePol,SNRdB,...
    rFix,azFix,nRand,xGrid,yGrid,iAllow,azAllow,verbose)
% [rBest,azPerf,azPerf,elPerf,hist] = ...
%    optRandConfig(egain,nTrials,lamEval,azTrue,elTrue,abcdBlock,...
%    defdB,azGrid,elGrid,truePol,SNRdB...
%    rFix,azFix,nRand,xGrid,yGrid,iAllow,azAllow,verbose);
% Optimizes a random set of planar array configurations that appear on 
% a fixed grid.  One specifies a set of fixed desired locations, and a set
% of allowed positions, with allowable azimuth rotations.  The given
% number of trials are run, with worst-case threshold performance evaluated
% over the given directions and wavelengths.  The best configuration
% is saved, along with its worst-case specifications.
% Returns an allowed planar array of apertures.
% egain       -- Handle for an element gain function @g such that
%                G = g(i,Azvals,Elvals,lam,ipol) must return the 
%                amplitude gains of element i in matrix G with same dims
%                as Azvals, Elvals, when these are in degrees.
%                ipol = 1 --> vertical (Z Electric Field), 2 --> Horizontal
%                polarization (+ electric field in cross(z,d) direction 
%                where d is the direction of propagation
% nTrials     -- number of random trials to use in the optimization
% lamEval     -- 1 x nLam wavelengths at which to evaluate each array
% azTrue      -- 1 x nTrue different azimuths to the true
%                directions to be evaluated for each wavelength, deg.
% elTrue      -- 1 x nTrue different elevations to true directions
% abcdBlock   -- 3 x 4 x m corners of rectangles, in order around edges
% defdB       -- the maximum allowable threshold for any config to be viable
%                If no threshold, then performance at this level is analyzed
% azGrid      -- nEl x nAz values of azimuth to search for sidelobes, deg.
%                azimuth constant down columns, elevation along rows
% elGrid      -- nEl x nAz values of elevation to search for sidelobes, deg.
% SNRdB       -- 1 x nsnr the SNR values for analysis.
%                SNR is with respect to a the whole array 1 snapshot,
%                i.e., SNRpower = qx * norm(A)^2 / sig^2, where qx is
%                the signal m.s. at a gain-of-1 antenna, norm(A) is the 
%                norm of the vector of gains of all apertures, 
%                and sig = noise rms
% rFix        -- 3 x nfix positions of fixed apertures
% azFix       -- 1 x nfix azimuth orientations of the fixed apertures, deg
% nRand       -- number of random apertures to place
% xGrid,yGrid -- nx x ny X and Y positions of the aperture location grid 
%                points. xGrid(i,:) is constant, yGrid(:,i) is constant
% iAllow      -- nx x ny logical array of which grid points are allowed
% azAllow     -- nx x ny x 2 allowed interval of azimuths for each node
%                [azAllow(i,j,1) azAllow(i,j,2)] is the interval in deg
%                Not necessarily in any particular branch cut.  E.G.,
%                suppose 10-20 deg is not allowed.  Put
%                azAllow(i,j,1:2) = [20 10+360];  Default is [0 360];
% verbose     -- logical, if present, yes/no for intermediate messages
% rBest       -- 3 x nfix+nrand coordinates of the best array apertures
% azBest      -- 1 x nfix+nrand azimuth rotations of the best array 
%                apertures in deg, in (-180,180]
% azPerf      -- nLam x 2 for each wavelength, the best worst-case azimuth
%                performance;  Worst case threshold in dB or defdB, 
%                whichever greater. and rms error at that value 
%                at defdB, whichever is greater
% elPerf      -- nLam x 2 same as azPerf but for elevation
% hist        -- cell array of processing messages of interest

nLam = length(lamEval);
M = 5;    % five snapshot analysis
nTrue = length(azTrue);

verbose = exist('verbose','var') && verbose;

%--------------------------------------------------------------------------
% Loop initialization
%--------------------------------------------------------------------------

iHist = 0;
hist = cell(0);

iDefSNR = find(SNRdB <= defdB,1,'last');    % default SNR index

azPerf = inf(nLam,2);
elPerf = inf(nLam,2);

%--------------------------------------------------------------------------
% Turn blockage rectangles into a useful structure for evaluation.
%--------------------------------------------------------------------------

recBlock = rec2struct(abcdBlock);

%--------------------------------------------------------------------------
% Loop over each independent trial.
%--------------------------------------------------------------------------

iHist = iHist + 1;
hist{iHist} = sprintf('Starting random optimiz. at %s',datestr(now));

iBest = 0;
rBest = [];
azBest = [];

if verbose
    fprintf('\n%s...',hist{iHist});
end

euler = zeros(3,nRand + size(rFix,2));

for iTry = 1:nTrials
    
    if (verbose && rem(iTry,25) == 0)
        fprintf('\niTry = %d...',iTry);
    end
    
    %----------------------------------------------------------------------
    % obtain the random configuration
    %----------------------------------------------------------------------
    
    [r,az] = randConfig(rFix,azFix,nRand,xGrid,yGrid,iAllow,azAllow);
    
    euler(1,:) = az;
    
    %----------------------------------------------------------------------
    % Loop over wavelengths.
    %----------------------------------------------------------------------
    
    azWorst = zeros(nLam,2);
    elWorst = azWorst;
    
    for iLam = 1:nLam
        
        lamChk = lamEval(iLam);
        
        %------------------------------------------------------------------
        % Loop over look directions.
        %------------------------------------------------------------------
        
        for iTrue = 1:nTrue
            
            azChk = azTrue(iTrue);
            elChk = elTrue(iTrue);
            dir = [cosd(azChk)*cosd(elChk); sind(azChk)*cosd(elChk); ...
                sind(elChk)];

            %--------------------------------------------------------------
            % Eliminate apertures that are blocked in this look direction.
            %--------------------------------------------------------------
            
            canSee = seePastBlock(r,dir,recBlock);
            rUse = r(:,canSee);
            eulerUse = euler(:,canSee);
            
            cont = any(canSee);
            if (~cont), break, end;    % nobody can see anything

            %--------------------------------------------------------------
            % Analyze, update worst case stats.
            %--------------------------------------------------------------

%             [Qchk,~,azThChk,elThChk] = arrayAthley(...
%                 egain,SNRdB,lamChk,rUse,eulerUse,M,[azChk elChk],...
%                 truePol,azGrid,elGrid,'Temp');
            
            [Qchk,~,azThChk,elThChk] = arrayAthley(...
                egain,SNRdB,lamChk,rUse,eulerUse,M,[azChk elChk],...
                truePol,azGrid,elGrid);
            
            if (isempty(azThChk) || azThChk(1) <= defdB)
                azThChk = [defdB sqrt(Qchk(1,1,iDefSNR))];
            end
            
            if (isempty(elThChk) || elThChk(1) <= defdB)
                elThChk = [defdB sqrt(Qchk(2,2,iDefSNR))];
            end
            
            if (azPerf(iLam,1) > defdB)   % haven't made threshold yet
                cont = azThChk(1) <= azPerf(iLam,1);
                
                if (cont && (azThChk(1) > azWorst(iLam,1) ...
                        || (azThChk(1) == azWorst(iLam,1) ...
                        && azThChk(2) > azWorst(iLam,2))))
                    azWorst(iLam,:) = azThChk;
                end
            else   % threshold ok, now look for performance
                cont = azThChk(1) == defdB && azThChk(2) <= azPerf(iLam,2);
                
                if (cont && azThChk(2) > azWorst(iLam,2))
                    azWorst(iLam,:) = azThChk;
                end
            end
            
            if (~cont), break, end;   % This config not optimum
            
            if (elPerf(iLam,1) > defdB)   % haven't made threshold yet
                cont = elThChk(1) <= elPerf(iLam,1);
                
                if (cont && (elThChk(1) > elWorst(iLam,1) ...
                        || (elThChk(1) == elWorst(iLam,1) ...
                        && elThChk(2) > elWorst(iLam,2))))
                    elWorst(iLam,:) = elThChk;
                end
            else   % threshold ok, now look for performance
                cont = elThChk(1) == defdB && elThChk(2) <= elPerf(iLam,2);
                
                if (cont && elThChk(2) > elWorst(iLam,2))
                    elWorst(iLam,:) = elThChk;
                end
            end
            
            if (~cont), break, end;   % This config not optimum
            
        end
        
        if (~cont), break, end;   % This config isn't optimum
        
    end
    
    newWinner = cont && iLam == nLam;
    
    if newWinner
        iBest = iTry;
        rBest = r;
        azBest = az;
        azPerf = azWorst;
        elPerf = elWorst;
        
        iHist = iHist + 1;
        hist{iHist} = sprintf('New best design:  #%d',iBest);
        
        if verbose
            fprintf('\n%s...',hist{iHist});
        end
        
        save('opRandTmp','iBest','rBest','azBest','azPerf','elPerf');
    end
    
end

iHist = iHist + 1;
hist{iHist} = sprintf('All done at %s, winner is %d!',datestr(now),iBest);

if verbose
    fprintf('\n%s\n',hist{iHist});
end

end
