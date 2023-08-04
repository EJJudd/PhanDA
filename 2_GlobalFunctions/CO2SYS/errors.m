% errors()
% This subroutine propagates uncertainties for the marine carbonate chemistry calculations
% from errors (or uncertainties) on inputs:
%  - pair of carbonate system variables
%  - nutrient (silicate and phosphate) concentrations
%  - concentrations of ammonium and hydrogen sulfide
%  - temperature and salinity
%  - calcium concentration (optional)
% plus errors in dissociation constants pK0, pK1, pK2, pKb, pKw, pKspa, and pKspc as well as total boron
%
% It calls derivnum, which computes numerical derivatives, and then
% it applies error propagation using the method of moments.
% The latter is a general technique to estimate the 2nd moment of a variable z
% (variance or standard deviation) based on a 1st-order approximation to z.
%
% This subroutine has been modified from its original version (Orr et al.
% 2018) to allow for input variables of CO2, HCO3, and CO3, the inclusion
% of ammonium and hydrogen sulfide, compatibility with CO2SYS.m(v3), and am
% optional input of calcium concentration uncertainty.
%
%**************************************************************************
%
%  **** SYNTAX:
%  [err, headers, units] = errors(PAR1,PAR2,PAR1TYPE,PAR2TYPE,...
%                                 SAL,TEMPIN,TEMPOUT,PRESIN,PRESOUT,SI,PO4,...
%                                 NH4,H2S,ePAR1,ePAR2,eSAL,eTEMP,eSI,ePO4,...
%                                 eNH4,eH2S,epK,eBt,r,pHSCALEIN,K1K2CONSTANTS,...
%                                 KSO4CONSTANT,KFCONSTANT,BORON,eCAL(optional))
% 
%  **** SYNTAX EXAMPLES:
%  [Result]                = errors(2400,2200,1,2,35,10,10,0,0,15,1,0,0,2,2,0.01,0.01,0,0,0,0,0,0,0,1,4,1,1,1)
%  [Result,Headers]        = errors(2400,   8,1,3,35,25,5,0,3000,15,1,0,0,2,0.001,0,0,0,0,0,0,0,0,0,1,4,1,1,1)
%  [Result,Headers,Units]  = errors(500,    8,5,3,35,25,5,0,4000,15,1,0,0,2,0.001,0,0,0,0,0,0,'','',0,1,4,1,1,1)
%  [A]                     = errors(2400,2000:10:2400,1,2,35,10,10,0,0,15,2,0,0,2,0,0,0,0,0,0,0,'','',0,1,4,1,1,1)
%  [A]                     = errors(2400,2200,1,2,0:1:35,0,25,4200,0,15,1,0,0,2,2,0,0,0,0,0,0,'','',0,1,4,1,1,1)
%  epK = [0.002, 0.0075, 0.015, 0.01, 0.01, 0.02, 0.02];
%  eBt = 0.02;
%  [A, hdr, units]         = errors(2400,2200,1,2,35,0,25,0:100:4200,0,15,1,0,0,2,2,0,0,0,0,0,0,epK,eBt,0,1,4,1,1,1)
%  
%**************************************************************************
%
% INPUT:
%
%   - ePAR1, ePAR2   :  uncertainty of PAR1 and PAR2 of input pair of CO2 system variables
%                       * Same units as PAR1 & PAR2
%   - eS, eT         :  uncertainty of Salinity and Temperature (same units as S and T)
%   - ePO4, eSI      :  uncertainty of Phosphate and Silicate total concentrations (same units as PO4 and SI [umol/kg])
%   - eNH4, eH2S     :  uncertainty of Ammonia and Hydrogen Sulfide total concentrations (same units as NH4 and H2S [umol/kg])
%   - epK            :  uncertainty of all seven dissociation constants (a vector) [pK units]
%   - eBt            :  uncertainty of total boron, given as fractional relative error (eBt=0.02 is a 2% error)
%   - r              :  correlation coefficient between PAR1 AND PAR2 (typicaly 0)
%   - others         :  same as input for subroutine  CO2SYS() (version 3, scalar or vectors)
%   - eCAL           :  uncertainty of calcium [mmol/kg] determined by ratio with salinity
%
% All parameters may be scalars or vectors except epK and eBt.
%   * epK must be vector of 7 values : errors of [pK0, pK1, pK2, pKb, pKw, pKspa, pKspc]. 
%     These errors are assumed to be the same for all rows of data.
%     These 7 values are in pK units
%
%     if epK is empty (= ''), this routine specifies default values.
%     These default standard errors are :
%        pK0   :  0.002 
%        pK1   :  0.0075
%        pK2   :  0.015
%        pKb   :  0.01    boric acid
%        pKw   :  0.01    water dissociation
%        pKspa :  0.02    solubility product of Aragonite 
%        pKspc :  0.02    solubility product of Calcite
%
%   * eBt is a scalar real number, fractional relative error (between 0.00 and 1.00)
%     for TB, where the default is eBt=0.02. It is assumed to be the same
%     for all rows of data.
%
% In constrast, ePAR1, ePAR2, eS, eT, ePO4, eSI, eNH4, and eH2S
%   - if vectors, are errors associated with each data point
%   - if scalars, are one error value associated to all data points
% The same for parameter 'r'.
%
% If no value is input for eCAL, it will not be evaluated.
%
% If 'r' is nonzero with a value between -1.0 and 1.0, it indicates the correlation 
% between uncertainties of the input pair of carbonate system variables.
% By default, 'r' is zero. However, for some pairs the user may want to specify a
% different value. For example, measurements of pCO2 and pH are often anti-correlated.
% The same goes for two other pairs: 'CO2 and CO3' and 'pCO2 and
% CO3'. But even for these cases, care is needed when using non-zero values of 'r'.
% 
% When the user propagates errors for an individual
% measurement, 'r' should ALWAYS be zero if each member of the input pair is
% measured independently. In this case, we are interested in the
% correlation between the uncertainties in those measurements, not in
% the correlation between the measurments themselves. Uncertainties from
% those measurements are probably not correlated if they come from
% different instruments. Conversely, if users are interested in the
% error in the mean of a distribution of measurements (i.e., if they are
% propagating standard errors instead of standard deviations), one
% should then also account for the correlation between the measurements of
% the two variables of the input pair.
% 
% For input pairs where one member is pH, this 'errors' routine automatically
% inverses the sign of 'r'.
% That inversion is done because the associated derivatives are computed in terms of 
% the hydrogen ion concentration H+, not pH. Therefore for each of these 6
% flags, if the user wants to compute 'r' that should be done (1) using
% the H+ concentration instead of pH, and (2) the sign of that computed 'r'
% should be inversed when passing it as an argument to this routine.
% To express perfect anticorrelation with pH, the user should 
% use 'r=-1.0'. 
% 
%**************************************************************************
%
% OUTPUT: * an array containing uncertainty for the following variables
%           (one row per sample):
%         *  a cell-array containing crudely formatted headers
%
%    POS  PARAMETER         UNIT
%
%    01 - TAlk              (umol/kgSW)
%    02 - TCO2              (umol/kgSW)
%    03 - [H+] in           (nmol/kgSW)
%    04 - pCO2 in           (uatm)
%    05 - fCO2 in           (uatm)
%    06 - HCO3 in           (umol/kgSW)
%    07 - CO3 in            (umol/kgSW)
%    08 - CO2 in            (umol/kgSW)
%    09 - RF in             ()
%    10 - OmegaCa in        ()
%    11 - OmegaAr in        ()
%    12 - xCO2 in           (ppm)
%    13 - [H+] out          (nmol/kgSW)
%    14 - pCO2 out          (uatm)
%    15 - fCO2 out          (uatm)
%    16 - HCO3 out          (umol/kgSW)
%    17 - CO3 out           (umol/kgSW)
%    18 - CO2 out           (umol/kgSW)
%    19 - RF out            ()
%    20 - OmegaCa out       ()
%    21 - OmegaAr out       ()
%    22 - xCO2 out          (ppm)
%
% NOTE: User-specified uncertainties for the input arguments are provided
%       in addition to uncertainties in output arguments at both input and
%       output conditions.
%

function [total_error, headers, units] = ...
        errors (PAR1, PAR2, PAR1TYPE, PAR2TYPE, SAL, TEMPIN, TEMPOUT, PRESIN, PRESOUT, SI, PO4,...
                NH4, H2S, ePAR1, ePAR2, eSAL, eTEMP, eSI, ePO4, eNH4, eH2S, epK, eBt, r, ...
                pHSCALEIN,K1K2CONSTANTS,KSO4CONSTANT,KFCONSTANT,BORON,eCAL)

    global K0 K1 K2 KW KB KF KS KP1 KP2 KP3 KSi KNH4 KH2S E;
    
    % Input conditioning
    % ------------------
    
    % Handle optional calcium error input
    if (nargin < 30), eCAL=0; end
    
    % Determine lengths of input vectors
    veclengths=[length(PAR1) length(PAR2) length(PAR1TYPE)...
                length(PAR2TYPE) length(SAL) length(TEMPIN)...
                length(TEMPOUT) length(PRESIN) length(PRESOUT)...
                length(SI) length(PO4) length(NH4) length(H2S)...
                length(ePAR1) length(ePAR2) length(eSAL) length(eTEMP)...
                length(eSI) length(ePO4) length(eNH4) length(eH2S)...
                length(r) length(pHSCALEIN) length(K1K2CONSTANTS)...
                length(KSO4CONSTANT) length(KFCONSTANT) length(BORON)...
                length(eCAL)];

    if length(unique(veclengths))>2
            disp(' '); disp('*** INPUT ERROR: Input vectors must all be of same length, or of length 1. ***'); disp(' '); return
    end

    % Make column vectors of all input vectors
    PAR1         =PAR1         (:);
    PAR2         =PAR2         (:);
    PAR1TYPE     =PAR1TYPE     (:);
    PAR2TYPE     =PAR2TYPE     (:);
    SAL          =SAL          (:);
    TEMPIN       =TEMPIN       (:);
    TEMPOUT      =TEMPOUT      (:);
    PRESIN       =PRESIN       (:);
    PRESOUT      =PRESOUT      (:);
    SI           =SI           (:);
    PO4          =PO4          (:);
    NH4          =NH4          (:);
    H2S          =H2S          (:);
    ePAR1        =ePAR1        (:);
    ePAR2        =ePAR2        (:);
    eSAL         =eSAL         (:);
    eTEMP        =eTEMP        (:);
    eSI          =eSI          (:);
    ePO4         =ePO4         (:);
    eNH4         =eNH4         (:);
    eH2S         =eH2S         (:);
    r            =r            (:);
    pHSCALEIN    =pHSCALEIN    (:);
    K1K2CONSTANTS=K1K2CONSTANTS(:);
    KSO4CONSTANT =KSO4CONSTANT (:);
    KFCONSTANT   =KFCONSTANT   (:);
    BORON        =BORON        (:);
    eCAL         =eCAL         (:);
    
    % Find the longest column vector:
    ntps = max(veclengths);

    % Populate column vectors
    PAR1(1:ntps,1)          = PAR1(:)          ;
    PAR2(1:ntps,1)          = PAR2(:)          ;
    PAR1TYPE(1:ntps,1)      = PAR1TYPE(:)      ;
    PAR2TYPE(1:ntps,1)      = PAR2TYPE(:)      ;
    SAL(1:ntps,1)           = SAL(:)           ;
    TEMPIN(1:ntps,1)        = TEMPIN(:)        ;
    TEMPOUT(1:ntps,1)       = TEMPOUT(:)       ;
    PRESIN(1:ntps,1)        = PRESIN(:)        ;
    PRESOUT(1:ntps,1)       = PRESOUT(:)       ;
    SI(1:ntps,1)            = SI(:)            ;
    PO4(1:ntps,1)           = PO4(:)           ;
    NH4(1:ntps,1)           = NH4(:)           ;
    H2S(1:ntps,1)           = H2S(:)           ;
    ePAR1(1:ntps,1)         = ePAR1(:)         ;
    ePAR2(1:ntps,1)         = ePAR2(:)         ;
    eSAL(1:ntps,1)          = eSAL(:)          ;
    eTEMP(1:ntps,1)         = eTEMP(:)         ;
    eSI(1:ntps,1)           = eSI(:)           ;
    ePO4(1:ntps,1)          = ePO4(:)          ;
    eNH4(1:ntps,1)          = eNH4(:)          ;
    eH2S(1:ntps,1)          = eH2S(:)          ;
    r(1:ntps,1)             = r(:)             ;
    pHSCALEIN(1:ntps,1)     = pHSCALEIN(:)     ;
    K1K2CONSTANTS(1:ntps,1) = K1K2CONSTANTS(:) ;
    KSO4CONSTANT(1:ntps,1)  = KSO4CONSTANT(:)  ;
    KFCONSTANT(1:ntps,1)    = KFCONSTANT(:)    ;
    BORON(1:ntps,1)         = BORON(:)         ;
    eCAL(1:ntps,1)          = eCAL(:)          ;

    % Exclude input variables equal to -999
    E = (PAR1~=-999 & PAR2~=-999);
    
    % Default values for epK
    if (isempty(epK))
        epK = [0.002, 0.0075, 0.015, 0.01, 0.01, 0.02, 0.02];
    else
        % Check validity of epK
        if (length(epK) == 1 && epK == 0)
            % this means that the caller does not want to account for errors on dissoc. constants
            epK = [0.0 0.0 0.0 0.0 0.0 0.0 0.0];
        elseif (length(epK) ~= 7)
            error ('invalid parameter epK: ', epK)
        end
    end
    
    % Default value for eBt (also check for incorrectly specified values
    if (isempty(eBt))
        eBt = 0.02;
    elseif ( ~(isscalar(eBt)) )
        error ('invalid parameter eBt (must be scalar): ')
    elseif ( isscalar(eBt))
        if (eBt < 0 || eBt > 1)
           error ('The "eBt" input argument is the fractional error. It must be between 0 and 1. Default is 0.02 (a 2% error).')
	end
    end

    % names of dissociation constants
    Knames = {'K0','K1','K2','Kb','Kw','Kspa', 'Kspc'};

    % Convert error on pH to error on [H+] concentration
    % in case where first input variable is pH
    isH = (E & PAR1TYPE == 3);
    
    pH = PAR1(isH);
    epH = ePAR1(isH);       % Error on pH
    H  = 10.^(-pH);         % H+ concentration
    r(isH) = -r(isH);       % Inverse sign of 'r' if PAR1 is pH

    % dpH = d(-log10[H])
    %     = d(- ln[H] / ln[10] )
    %     = -(1/ln[10]) * d (ln[H])
    %     = -(1/ln[10]) * (dH / H)
    % Thus dH = - ln[1O] * [H] dpH
    eH =  log(10) * (H .* epH);     % Removed the minus sign because all errors (sigmas) are positive by definition
    eH =  eH * 1e9            ;     % Convert from mol/kg to nmol/kg (to have same units as partial derivative)
    ePAR1(isH) = eH;

    % Same conversion for second variable
    isH = (E & PAR2TYPE == 3);
    pH = PAR2(isH);
    epH = ePAR2(isH);       % Error on pH
    H  = 10.^(-pH);         % H+ concentration
    r(isH) = -r(isH);       % Inverse sign of 'r' if PAR2 is pH

    eH =   log(10) * (H .* epH);
    eH =  eH * 1e9;
    ePAR2(isH) = eH;

    % Convert calcium error from mmol/kg to mol/kg
    eCAL = eCAL.*1e-3;

    % initialise total square error
    sq_err = zeros(ntps,1);
    
    % Contribution of PAR1 to squared standard error
    if (any(ePAR1 ~= 0.0))
        % Compute sensitivities (partial derivatives)
        [deriv1(E,:),~,~,headers_err,units_err] = derivnum ('PAR1',...
            PAR1(E),PAR2(E),PAR1TYPE(E),PAR2TYPE(E),SAL(E),TEMPIN(E),...
            TEMPOUT(E),PRESIN(E),PRESOUT(E),SI(E),PO4(E),NH4(E),H2S(E),...
            pHSCALEIN(E),K1K2CONSTANTS(E),KSO4CONSTANT(E),KFCONSTANT(E),BORON(E));
        err = bsxfun(@times,deriv1,ePAR1);
	 sq_err = bsxfun(@plus,err*0., sq_err);
     sq_err = sq_err + err .* err;
    end

    % Contribution of PAR2 to squared standard error
    if (any (ePAR2 ~= 0.0))
        % Compute sensitivities (partial derivatives)
        [deriv2(E,:),~,~,headers_err,units_err] = derivnum ('PAR2',...
            PAR1(E),PAR2(E),PAR1TYPE(E),PAR2TYPE(E),SAL(E),TEMPIN(E),...
            TEMPOUT(E),PRESIN(E),PRESOUT(E),SI(E),PO4(E),NH4(E),H2S(E),...
            pHSCALEIN(E),K1K2CONSTANTS(E),KSO4CONSTANT(E),KFCONSTANT(E),BORON(E));
        err = bsxfun(@times,deriv2,ePAR2);
	 sq_err = bsxfun(@plus,err*0., sq_err);
     sq_err = sq_err + err .* err;
    end

    % Contribution of covariance of PAR1 and PAR2 to squared standard error
    covariance  = nan(ntps,1);
    if (any (r ~= 0.0) && any (ePAR1 ~= 0.0) && any (ePAR2 ~= 0.0))
        % Compute covariance from correlation coeff. & std deviations
        covariance(E) = r(E) .* ePAR1(E) .* ePAR2(E);
        % Contribution to squared error
        err2 = bsxfun(@times,2 * deriv1 .* deriv2, covariance);
        sq_err = sq_err + err2;
    end
    
    % Contribution of Silion (total dissolved inorganic concentration) to squared standard error
    %
    % Remark : does not compute error where SI = 0 
    %          because computation of sensitivity to SI fails in that case
    %
    SI_valid = (E & SI ~= 0 & eSI ~= 0);
    if (any (SI_valid))
        % Compute sensitivities (partial derivatives)
        [deriv, ~, ~, headers_err, units_err] = derivnum ('sil',PAR1(SI_valid),PAR2(SI_valid),PAR1TYPE(SI_valid),PAR2TYPE(SI_valid),...
                   SAL(SI_valid),TEMPIN(SI_valid),TEMPOUT(SI_valid),PRESIN(SI_valid),PRESOUT(SI_valid),...
                   SI(SI_valid),PO4(SI_valid),NH4(SI_valid),H2S(SI_valid),pHSCALEIN(SI_valid),K1K2CONSTANTS(SI_valid),KSO4CONSTANT(SI_valid),...
                   KFCONSTANT(SI_valid),BORON(SI_valid));
        err = bsxfun(@times,deriv,eSI(SI_valid));
        new_size = [ntps size(err,2)];
	sq_err = zeros(new_size) + sq_err;
        sq_err(SI_valid,:) = sq_err(SI_valid,:) + err .* err;
    end

    % Contribution of Phosphorus (total dissoloved inorganic concentration) to squared standard error
    %
    % Remark : does not compute error where PO4 = 0 
    %          because computation of sensitivity to PO4 fails in that case
    %
    PO4_valid = (E & PO4 ~= 0 & ePO4 ~= 0);
    if (any (PO4_valid))
        % Compute sensitivities (partial derivatives)
        [deriv, ~, ~, headers_err, units_err] = derivnum ('phos',PAR1(PO4_valid),PAR2(PO4_valid),PAR1TYPE(PO4_valid),PAR2TYPE(PO4_valid),...
                   SAL(PO4_valid),TEMPIN(PO4_valid),TEMPOUT(PO4_valid),PRESIN(PO4_valid),PRESOUT(PO4_valid),...
                   SI(PO4_valid),PO4(PO4_valid),NH4(PO4_valid),H2S(PO4_valid),pHSCALEIN(PO4_valid),K1K2CONSTANTS(PO4_valid),KSO4CONSTANT(PO4_valid),...
                   KFCONSTANT(PO4_valid),BORON(PO4_valid));
        err = bsxfun(@times,deriv,ePO4(PO4_valid));
        new_size = [ntps size(err,2)];
	sq_err = zeros(new_size) + sq_err;
        sq_err(PO4_valid,:) = sq_err(PO4_valid,:) + err .* err;
    end
    
    % Contribution of Ammonia (total dissolved inorganic concentration) to squared standard error
    %
    % Remark : does not compute error where NH4 = 0 
    %          because computation of sensitivity to SI fails in that case
    %
    NH4_valid = (E & NH4 ~= 0 & eNH4 ~= 0);
    if (any (NH4_valid))
        % Compute sensitivities (partial derivatives)
        [deriv, ~, ~, headers_err, units_err] = derivnum ('amm',PAR1(NH4_valid),PAR2(NH4_valid),PAR1TYPE(NH4_valid),PAR2TYPE(NH4_valid),...
                   SAL(NH4_valid),TEMPIN(NH4_valid),TEMPOUT(NH4_valid),PRESIN(NH4_valid),PRESOUT(NH4_valid),...
                   SI(NH4_valid),PO4(NH4_valid),NH4(NH4_valid),H2S(NH4_valid),pHSCALEIN(NH4_valid),K1K2CONSTANTS(NH4_valid),KSO4CONSTANT(NH4_valid),...
                   KFCONSTANT(NH4_valid),BORON(NH4_valid));
        err = bsxfun(@times,deriv,eNH4(NH4_valid));
        new_size = [ntps size(err,2)];
	sq_err = zeros(new_size) + sq_err;
        sq_err(NH4_valid,:) = sq_err(NH4_valid,:) + err .* err;
    end

    % Contribution of Hydrogen Sulfide (total dissoloved inorganic concentration) to squared standard error
    %
    % Remark : does not compute error where H2S = 0 
    %          because computation of sensitivity to PO4 fails in that case
    %
    H2S_valid = (E & H2S ~= 0 & eH2S ~= 0);
    if (any (H2S_valid))
        % Compute sensitivities (partial derivatives)
        [deriv, ~, ~, headers_err, units_err] = derivnum ('hyd',PAR1(H2S_valid),PAR2(H2S_valid),PAR1TYPE(H2S_valid),PAR2TYPE(H2S_valid),...
                   SAL(H2S_valid),TEMPIN(H2S_valid),TEMPOUT(H2S_valid),PRESIN(H2S_valid),PRESOUT(H2S_valid),...
                   SI(H2S_valid),PO4(H2S_valid),NH4(H2S_valid),H2S(H2S_valid),pHSCALEIN(H2S_valid),K1K2CONSTANTS(H2S_valid),KSO4CONSTANT(H2S_valid),...
                   KFCONSTANT(H2S_valid),BORON(H2S_valid));
        err = bsxfun(@times,deriv,eH2S(H2S_valid));
        new_size = [ntps size(err,2)];
	sq_err = zeros(new_size) + sq_err;
        sq_err(H2S_valid,:) = sq_err(H2S_valid,:) + err .* err;
    end

    % Contribution of T (temperature) to squared standard error
    if (any (eTEMP ~= 0.0))
        % Compute sensitivities (partial derivatives)
        [deriv(E,:),~,~,headers_err,units_err] = derivnum ('T',...
            PAR1(E),PAR2(E),PAR1TYPE(E),PAR2TYPE(E),SAL(E),TEMPIN(E),...
            TEMPOUT(E),PRESIN(E),PRESOUT(E),SI(E),PO4(E),NH4(E),H2S(E),...
            pHSCALEIN(E),K1K2CONSTANTS(E),KSO4CONSTANT(E),KFCONSTANT(E),BORON(E));
        err = bsxfun(@times,deriv,eTEMP);
	 sq_err = err*0. + sq_err;
     sq_err = sq_err + err .* err;
    end

    % Contribution of S (salinity) to squared standard error
    if (any (eSAL ~= 0.0))
        % Compute sensitivities (partial derivatives)
        [deriv(E,:),~,~,headers_err,units_err] = derivnum ('S',...
            PAR1(E),PAR2(E),PAR1TYPE(E),PAR2TYPE(E),SAL(E),TEMPIN(E),...
            TEMPOUT(E),PRESIN(E),PRESOUT(E),SI(E),PO4(E),NH4(E),H2S(E),...
            pHSCALEIN(E),K1K2CONSTANTS(E),KSO4CONSTANT(E),KFCONSTANT(E),BORON(E));
        err = bsxfun(@times,deriv,eSAL);
	 sq_err = err*0. + sq_err;
     sq_err = sq_err + err .* err;
    end

    % Calculate dissociation constants
    data = CO2SYS(PAR1,PAR2,PAR1TYPE,PAR2TYPE,SAL,TEMPIN,TEMPOUT,PRESIN,PRESOUT,SI,PO4,NH4,H2S,pHSCALEIN,K1K2CONSTANTS,KSO4CONSTANT,KFCONSTANT,BORON);
        
    % Calculate [Ca++]
    % '       Riley, J. P. and Tongudai, M., Chemical Geology 2:263-269, 1967:
    % '       this is .010285.*Sali./35
    Ca = 0.02128./40.087.*(SAL./1.80655);% ' in mol/kg-SW

    % Contribution of all pKi to squared standard error
    for i = 1:length(epK)

        % if error on Ki is given
        if (epK(i) ~= 0.0)
            % Select Ki
            switch i
                case 1
                  Ki_in = data(:,63);    % K0in
                  Ki_out = data(:,78);   % K0out
                case 2
                  Ki_in = data(:,64);    % K1in
                  Ki_out = data(:,79);   % K1out
                case 3
                  Ki_in = data(:,65);    % K2in
                  Ki_out = data(:,80);   % K2out
                case 4
                  Ki_in = data(:,69);    % KBin
                  Ki_out = data(:,84);   % KBout
                case 5
                  Ki_in = data(:,68);    % KWin
                  Ki_out = data(:,83);   % KWout
                case 6
                  % Recompute KAr from OmegaAr and ions [Ca++] and [CO3--] concentrations
                  OmegaAr = data(:,18);
                  CO3 = data(:,7) * 1e-6;
                  Ki_in = CO3.*Ca./OmegaAr;  % KspAin
		          Ki_in(isnan(Ki_in)) = 0;
                  OmegaAr = data(:,36);
                  CO3 = data(:,25) * 1e-6;
                  Ki_out = CO3.*Ca./OmegaAr; % KspAout
		          Ki_out(isnan(Ki_out)) = 0;
                case 7
                  % Recompute KCa from OmegaCa and ions [Ca++] and [CO3--] concentrations
                  OmegaCa = data(:,17);
                  CO3 = data(:,7) * 1e-6;
                  Ki_in = CO3.*Ca./OmegaCa;  % KspCin
                  Ki_in(isnan(Ki_in)) = 0;
                  OmegaCa = data(:,35);
                  CO3 = data(:,25) * 1e-6;
                  Ki_out = CO3.*Ca./OmegaCa; % KspCout
                  Ki_out(isnan(Ki_out)) = 0;
            end

            % compute error on Ki from that on pKi
            eKi_in  = - epK(i) * Ki_in * log(10);  % error on K at input conditions
            eKi_out = - epK(i) * Ki_out * log(10); % error on K at output conditions
            % Compute sensitivities (partial derivatives)
            [deriv(E,:),~,~,headers_err,units_err] = derivnum (cell2mat(Knames(1,i)),...
                       PAR1(E),PAR2(E),PAR1TYPE(E),PAR2TYPE(E),...
                       SAL(E),TEMPIN(E),TEMPOUT(E),PRESIN(E),PRESOUT(E),...
                       SI(E),PO4(E),NH4(E),H2S(E),pHSCALEIN(E),K1K2CONSTANTS(E),...
                       KSO4CONSTANT(E),KFCONSTANT(E),BORON(E));
            % Derivatives at input conditions are calculated with respect
            % to constants at input conditions, and derivatives at output
            % conditions are calculated with respect to constants at output
            % conditions. So, error terms are calculated using absolute
            % errors in K values at both input and output conditions.
            err(:, 1:12) = bsxfun(@times, deriv(:, 1:12), eKi_in);
            err(:, 13:22) = bsxfun(@times, deriv(:, 13:22), eKi_out);
	     sq_err = err*0. + sq_err;
         sq_err = sq_err + err .* err;
        end
    end

    % Contribution of Boron (total dissoloved boron concentration) to squared standard error
    if (eBt ~= 0)
        % Compute sensitivities (partial derivatives)
        [deriv(E,:),~,~,headers_err,units_err] = derivnum ('bor',...
            PAR1(E),PAR2(E),PAR1TYPE(E),PAR2TYPE(E),SAL(E),TEMPIN(E),...
            TEMPOUT(E),PRESIN(E),PRESOUT(E),SI(E),PO4(E),NH4(E),H2S(E),...
            pHSCALEIN(E),K1K2CONSTANTS(E),KSO4CONSTANT(E),KFCONSTANT(E),BORON(E));
     err = bsxfun(@times, deriv, (eBt*data(:,93)*1e-6)); % where TB = data(:,93) in umol B/kg
     new_size = [ntps size(err,2)];
	 sq_err = zeros(new_size) + sq_err;
     sq_err = sq_err + err .* err;
    end
    
    % Contribution of Calcium (total dissoloved calcium concentration) to squared standard error
    CAL_valid = (E & eCAL ~= 0);
    if (any (CAL_valid))
        % Compute sensitivities (partial derivatives)
        [deriv,~,~,headers_err,units_err] = derivnum ('cal',...
            PAR1(CAL_valid),PAR2(CAL_valid),PAR1TYPE(CAL_valid),PAR2TYPE(CAL_valid),SAL(CAL_valid),TEMPIN(CAL_valid),...
            TEMPOUT(CAL_valid),PRESIN(CAL_valid),PRESOUT(CAL_valid),SI(CAL_valid),PO4(CAL_valid),NH4(CAL_valid),H2S(CAL_valid),...
            pHSCALEIN(CAL_valid),K1K2CONSTANTS(CAL_valid),KSO4CONSTANT(CAL_valid),KFCONSTANT(CAL_valid),BORON(CAL_valid));
    err = bsxfun(@times,deriv,eCAL(CAL_valid));
    new_size = [ntps size(err,2)];
	sq_err = zeros(new_size) + sq_err;
    sq_err(CAL_valid,:) = sq_err(CAL_valid,:) + err .* err;
    end

    % Compute and return resulting total error (or uncertainty)
    total_error = sqrt (sq_err);
    total_error(isnan(total_error))=-999;
    headers = strcat('u(',headers_err,')');
    units = strcat('(',units_err,')');
end
