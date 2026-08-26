function dualfutile_ma_nohopf()
% To return detailed statistics, replace the first line with
% function stats = dualfutile_ma_nohopf()
% and uncomment the corresponding statistics block below.



%DUALFUTILE_MA_NOHOPF Exact mass-action positivity certificate.
%
% The code constructs
%
%   TtildeStar = delta3^12*xi1^5*eta1^5*
%                (1+xi2)^9*(1+eta2)^9*Ttilde
%
% and verifies coefficientwise that
%
%   TtildeStar - Gamma1*E - Gamma2*Esigma >=_coeff 0.
%
% Every integer operation is certified to remain below flintmax=2^53.
%
% The order of this .m file, after the main function dualfutile_ma_nohopf, is:
% 1) paper-specific auxiliary functions, approximately in their order of use;
% 2) general polynomial, exact-arithmetic, and packed-exponent utilities.

    totalTimer=tic;
    clc;
    global EXACT_STATS VERIFY_EXACT_ARITHMETIC
    VERIFY_EXACT_ARITHMETIC = true;
    EXACT_STATS = init_exact_stats();

    fprintf('============================================================\n');
    fprintf('DUAL FUTILE CYCLE: MASS-ACTION POSITIVITY CERTIFICATE\n');
    fprintf('============================================================\n');
    fprintf('Ttilde = P5*Delta4 - P6*B1plus - N5*B2plus\n');
    fprintf('Packed exponents: 16 fields x 4 bits; overflow checks enabled.\n');
%   fprintf('Exact denominator-clearance checks enabled.\n\n');

    %% 1. Construct Ttilde(a) from the paper's quadratic eigenvalue problem.
    C = build_target_numerator();

    assert(numel(C.k) == 5318338, ...
        'Unexpected support size for the corrected source polynomial.');
    assert(nnz(C.c > 0) == 4737696 && nnz(C.c < 0) == 580642, ...
        'Unexpected sign counts for the corrected source polynomial.');

    fprintf('Ttilde(a) source expansion:\n');
    fprintf('  monomials             = %s\n',comma_integer(numel(C.k)));
    fprintf('  positive coefficients = %s\n',comma_integer(nnz(C.c>0)));
    fprintf('  negative coefficients = %s\n\n',comma_integer(nnz(C.c<0)));

%% 2. Isolate the five-variable substitution from mu' to mu, 
%     collecting the remaining monomial factors into a sparse source matrix A 
%     and updating the exponent of delta3 (in baseKey) induced by the substitution.    

[A,baseKey,NR,NT] = paper_source_matrix(C);
  
assert(NR==5 && NT==5,'Unexpected delta1 or delta4 degree in the mass-action source.');

 
 %  Comment-in if these implementation details are wanted in the output
 %  fprintf('Collected source table:\n');
 %  fprintf('  source monomial groups = %s\n',comma_integer(size(A,1)));
 %  fprintf('  nonzero source entries = %s\n\n',comma_integer(nnz(A)));

 %% 3. Construct the mu'->mu substitution and the 2555544 clearance transform.
 %  Tmass maps source exponents patterns to exp. of xi1 xi2 eta1 eta2
    Tmass = mass_action_clearance_transform();
 %  Sanity checks, a priori:
    assert(isequal(size(Tmass),[900 3600]));
    check_integer_matrix(Tmass,'mass-action clearance transform');
 %  Comment-in if these implementation details are wanted in the output
 %  fprintf('Mass-action transform: 900 x 3600, nnz = %s\n\n', ...
 %       comma_integer(nnz(Tmass)));

  %% 4. Construct E and verify all quadratic identities proved in the paper.
    fprintf('Constructing E and E^sigma from the formulas in the paper ...\n');
    E = paper_E();
    Esigma = swap_arms(E);
    assert(nnz(E.c<0)==39 && nnz(Esigma.c<0)==39, ...
    'Unexpected adverse support in E or E^sigma.');

    fprintf('Construction of E and E^sigma complete.\n\n');

    %% 5. Construct the xi,eta factors in Gamma1,Gamma2.
    %   As for E above, Gamma1Gamm2 are constructed by the definition in the paper 
    % not from Ttildestar
   
    fprintf('Constructing Gamma1 and Gamma2 from the formulas in the paper ...\n');
    [Gamma1Ratios,Gamma2Ratios] = gamma_ratio_factors();
    assert(all(Gamma1Ratios.c>0) && all(Gamma2Ratios.c>0), ...
        'Gamma1 and Gamma2 ratio factors must have positive coefficients.');
  
   % Comment in, if of interest
   % fprintf('Structural factors:\n');
   % fprintf('  E support       = %d (%d adverse monomials)\n',numel(E.k),nnz(E.c<0));
   % fprintf('  E^sigma support = %d (%d adverse monomials)\n',numel(Esigma.k),nnz(Esigma.c<0));
   % fprintf('  Gamma1 ratio factor = %d positive monomials\n',numel(Gamma1Ratios.c));
   % fprintf('  Gamma2 ratio factor = %d positive monomials\n\n',numel(Gamma2Ratios.c));

    %% 6. Remaining monomial factors in Gamma1,Gamma2.
    % Packed source order:
    % 1 delta1, 2 delta3, 3 delta4, 4 alpha, 5 tau1,
    % 6 theta1, 7 tau2, 8 gamma, 9 theta2, 10 j,
    % 11 d, 12 epsilon, 13 z/tau1=w/tau2, 14 r/theta1=u/theta2.
    e1=zeros(1,16);
    e1([2 4 5 8 9 10 12]) = [4 1 1 1 1 4 3];
    Gamma1Base = struct('k',packrow(e1),'c',1);

    e2=zeros(1,16);
    e2([2 4 5 8 9 10 12]) = [5 1 1 1 1 3 4];
    Gamma2Base = struct('k',packrow(e2),'c',1);

    Gamma1BaseE = pmul(Gamma1Base,E);
    Gamma2BaseE = pmul(Gamma2Base,Esigma);

    %% 7. Build Gamma1*E + Gamma2*E^sigma.
    certificateTerms = certificate_matrix(baseKey, ...
    Gamma1BaseE,Gamma1Ratios,Gamma2BaseE,Gamma2Ratios);

fprintf('Construction of Gamma1*E + Gamma2*E^sigma complete.\n\n');

    % The negative source groups must be exactly those predicted by E,E^sigma.
    expectedNegKeys = unique([Gamma1BaseE.k(Gamma1BaseE.c<0); ...
                              Gamma2BaseE.k(Gamma2BaseE.c<0)]);
    assert(numel(expectedNegKeys)==78, ...
        'The expected 39+39 adverse supports overlap unexpectedly.');

    %% 8. Prepare exact clearance checks after full coefficient collection.
    clearance = init_clearance_minimality();
    [EvalXi1,EvalEta1,EvalXi2,EvalEta2,xi2Exp,eta2Exp] = ...
        clearance_test_operators();

    %% 9. Check TtildeStar-Gamma1*E-Gamma2*E^sigma coefficientwise.
    % this is the actual computer-assisted positivity certificate
    % it certifies that TtildeStar>Sigma^+>=0

    nRows=size(A,1);
    block=250;

    rawNegOccurrences=0;
    rawNegRows=false(nRows,1);
    rawPositiveOccurrences=0;
    residualNegativeOccurrences=0;
    residualPositiveOccurrences=0;
    maxAbsTstar=0;
    maxAbsRemainder=0;

    fprintf('Checking that the coefficientwise certificate\n');
    fprintf('  SigmaPlus = TtildeStar - Gamma1*E - Gamma2*Esigma\n');
    fprintf('has only nonnegative coefficients ...\n');

    nextProgress=10;

    for first=1:block:nRows
        rows=first:min(first+block-1,nRows);

        TstarBlk = exact_matrix_product(A(rows,:),Tmass, ...
            'TtildeStar block product');

        clearance = update_clearance_minimality( ...
            clearance,TstarBlk,baseKey(rows),rows, ...
            EvalXi1,EvalEta1,EvalXi2,EvalEta2,xi2Exp,eta2Exp);

        certificateBlk = certificateTerms(rows,:);
        certify_subtraction(TstarBlk,certificateBlk, ...
            'TtildeStar-Gamma1*E-Gamma2*Esigma');
        remainderBlk = TstarBlk-certificateBlk;

        gv=nonzeros(TstarBlk);
        rv=nonzeros(remainderBlk);
        if ~isempty(gv), maxAbsTstar=max(maxAbsTstar,max(abs(gv))); end
        if ~isempty(rv), maxAbsRemainder=max(maxAbsRemainder,max(abs(rv))); end

        negG=(TstarBlk<0);
        negR=(remainderBlk<0);

        rawNegOccurrences=rawNegOccurrences+nnz(negG);
        rawNegRows(rows)=full(any(negG,2));
        rawPositiveOccurrences=rawPositiveOccurrences+nnz(TstarBlk>0);

        residualNegativeOccurrences=residualNegativeOccurrences+nnz(negR);
        residualPositiveOccurrences=residualPositiveOccurrences+nnz(remainderBlk>0);

        assert(~any(negR(:)), ...
            sprintf('The certificate remainder has a negative coefficient in rows %d:%d.', ...
                    rows(1),rows(end)));

        progress=100*rows(end)/nRows;
        if progress>=nextProgress
            fprintf('  %d%% complete\n',nextProgress);
            nextProgress=nextProgress+10;
        end
        
    end

    %% 10. Exact adverse-support checks.
    % These regression checks (not strictly needed) verify that the adverse
    % base rows of TtildeStar are exactly as expected by E and Esigma. For 
    % simplicity - we do not compare the complete expanded negative supports.
    % Step 9 is already itself sufficient for the positivity conclusion

    actualNegKeys=baseKey(rawNegRows);
    assert(numel(actualNegKeys)==78, ...
        sprintf('Expected 78 adverse source groups; found %d.', ...
                numel(actualNegKeys)));
    assert(isequal(sort(actualNegKeys),sort(expectedNegKeys)), ...
        'Raw adverse rows are not exactly the E/E^sigma supports.');

    assert(rawNegOccurrences==1596, ...
        sprintf('Expected 1596 raw adverse occurrences; found %d.', ...
                rawNegOccurrences));

    [tf1,id1]=ismember(Gamma1BaseE.k(Gamma1BaseE.c<0),baseKey);
    [tf2,id2]=ismember(Gamma2BaseE.k(Gamma2BaseE.c<0),baseKey);
    assert(all(tf1) && all(tf2));

    G1=exact_matrix_product(A(id1,:),Tmass,'Gamma1 adverse count');
    G2=exact_matrix_product(A(id2,:),Tmass,'Gamma2 adverse count');
    n1=nnz(G1<0); n2=nnz(G2<0);
    assert(n1==798 && n2==798, ...
        sprintf('Expected a 798+798 adverse split; found %d+%d.',n1,n2));

    assert(residualNegativeOccurrences==0, ...
        'Coefficientwise nonnegativity of the certificate remainder failed.');
    assert(residualPositiveOccurrences>0, ...
        'The certificate remainder is identically zero.');

    %% 11. Complete the exact clearance and minimality checks.
clearance=finalize_clearance_minimality(clearance);

% Uncomment to print the detailed clearance-minimality report.
% The positivity proof requires denominator clearance, but not minimality.
% report_clearance_minimality(clearance);

%% 12. Complete the exact-arithmetic checks and report the conclusion.
report_exact_stats();

fprintf('\n============================================================\n');
fprintf('CERTIFICATE VERIFIED\n');
fprintf('============================================================\n');
fprintf('Exact-integer arithmetic below flintmax          : PASS\n');
fprintf('Packed-exponent overflow checks                  : PASS\n');
fprintf('TtildeStar-Gamma1*E-Gamma2*Esigma >=_coeff 0     : PASS\n');
fprintf('Therefore TtildeStar is positive on the positive orthant.\n\n');

% Uncomment to print detailed coefficient statistics.
%{
totalMonomials=rawPositiveOccurrences+rawNegOccurrences;

fprintf('Monomials of TtildeStar checked                  : %s\n', ...
    comma_integer(totalMonomials));
fprintf('Negative coefficients before subtraction        : %d\n', ...
    rawNegOccurrences);
fprintf('Adverse support split                            : %d + %d\n',n1,n2);
fprintf('Negative coefficients after subtraction         : %d\n', ...
    residualNegativeOccurrences);
fprintf('Positive coefficients after subtraction         : %s\n', ...
    comma_integer(residualPositiveOccurrences));
fprintf('max |coefficient| in TtildeStar                  : %.0f\n', ...
    maxAbsTstar);
fprintf('max |coefficient| in certificate remainder       : %.0f\n', ...
    maxAbsRemainder);
%}

stats=struct();


% Uncomment to return detailed coefficient statistics.
%{
stats.sourceMonomials=numel(C.k);
stats.sourcePositiveCoefficients=nnz(C.c>0);
stats.sourceNegativeCoefficients=nnz(C.c<0);
stats.sourceMonomialGroups=nRows;
stats.massActionTransformNNZ=nnz(Tmass);
stats.TtildeStarMonomials=rawPositiveOccurrences+rawNegOccurrences;
stats.TtildeStarPositiveCoefficients=rawPositiveOccurrences;
stats.TtildeStarNegativeCoefficients=rawNegOccurrences;
stats.Gamma1NegativeOccurrences=n1;
stats.Gamma2NegativeOccurrences=n2;
stats.remainderPositiveCoefficients=residualPositiveOccurrences;
stats.remainderNegativeCoefficients=residualNegativeOccurrences;
stats.maxAbsTtildeStar=maxAbsTstar;
stats.maxAbsRemainder=maxAbsRemainder;
%}

% Uncomment to return the detailed clearance-minimality statistics.
%{
stats.clearanceFactorNames=clearance.names;
stats.clearanceChosenExponents=clearance.chosen;
stats.clearanceNumeratorMultiplicities=clearance.multiplicity;
stats.clearanceReducedExponents=clearance.reduced;
stats.clearanceMinimalLaurentFreeExponents=clearance.minimalLaurentFree;
stats.clearanceChosenIsMinimalLaurentFree= ...
    clearance.chosenIsMinimalLaurentFree;
stats.clearanceWitnessRows=clearance.row;
stats.clearanceWitnessValues=clearance.value;
%}

elapsedTime=toc(totalTimer);
fprintf('Total elapsed time                              : %.1f seconds\n', ...
    elapsedTime);
end
% =======================================================================
% Direct construction of Ttilde(a)
% =======================================================================

function C = build_target_numerator()
    buildTimer = tic;
    fprintf('Constructing det(lambda^2*D+lambda*DN-DLR) ...\n');

    id.delta1=1; id.delta3=2; id.delta4=3; id.alpha=4; id.tau1=5;
    id.theta1=6; id.tau2=7; id.gamma=8; id.theta2=9; id.j=10;
    id.d=11; id.epsilon=12; id.zOverTau1=13; id.rOverTheta1=14;

    O = pconst(1);
    delta1=pvar(id.delta1); delta3=pvar(id.delta3); delta4=pvar(id.delta4);
    alpha=pvar(id.alpha); tau1=pvar(id.tau1); theta1=pvar(id.theta1);
    tau2=pvar(id.tau2); gamma=pvar(id.gamma); theta2=pvar(id.theta2);
    j=pvar(id.j); d=pvar(id.d); epsilon=pvar(id.epsilon);
    zOverTau1=pvar(id.zOverTau1); rOverTheta1=pvar(id.rOverTheta1);

    D = zero_matrix(4,4);
    D{1,1}=delta1; D{2,2}=O; D{3,3}=delta3; D{4,4}=delta4;

    DN = zero_matrix(4,4);
    DN{1,1}=psum(alpha,j,tau1);       DN{1,2}=j;
    DN{2,1}=j;                        DN{2,2}=psum(d,j,theta1);
    DN{2,3}=d;                        DN{3,2}=d;
    DN{3,3}=psum(d,epsilon,tau2);     DN{3,4}=epsilon;
    DN{4,3}=epsilon;                  DN{4,4}=psum(epsilon,gamma,theta2);

    DL = zero_matrix(4,2);
    DL{1,1}=alpha;
    DL{2,1}=pscale(d,-1); DL{2,2}=pscale(d,-1);
    DL{3,1}=pscale(d,-1); DL{3,2}=pscale(d,-1);
    DL{4,2}=gamma;

    R = zero_matrix(2,4);
    R{1,1}=pscale(pprod(zOverTau1,tau1),-1);
    R{1,3}=pprod(zOverTau1,tau2);
    R{2,2}=pprod(rOverTheta1,theta1);
    R{2,4}=pscale(pprod(rOverTheta1,theta2),-1);

    minusDLR=pmatmul(DL,R);
    for r=1:4
        for s=1:4
            minusDLR{r,s}=pscale(minusDLR{r,s},-1);
        end
    end

    Qs = cell(4,4);
    for r=1:4
        for s=1:4
            Qs{r,s} = {minusDLR{r,s}, DN{r,s}, D{r,s}};
        end
    end

    dets = pdet4_s(Qs);                 % coefficients of s^0,...,s^8
    assert(isempty(dets{1}.k) && isempty(dets{2}.k), ...
        'The constant and linear QEP coefficients must vanish.');
    a = cell(1,7);
    for r=0:6
        a{r+1}=dets{9-r};               % a0,...,a6
    end

    expected = [1 12 55 126 154 90 25];
    actual = cellfun(@(x)numel(x.k),a);
    assert(isequal(actual,expected),'Characteristic-coefficient check failed.');
    assert_poly_equal(a{1},pprod(delta1,delta3,delta4), ...
        'a0 must equal det(D)=delta1*delta3*delta4.');

    fprintf('Constructing P5, P6 and Ttilde(a) ...\n');

    
% Compute the positive and negative parts directly from a5 and a6.
P5h = struct( ...
    'k',a{6}.k(a{6}.c>0), ...
    'c',a{6}.c(a{6}.c>0));

N5h = struct( ...
    'k',a{6}.k(a{6}.c<0), ...
    'c',-a{6}.c(a{6}.c<0));

P6h = struct( ...
    'k',a{7}.k(a{7}.c>0), ...
    'c',a{7}.c(a{7}.c>0));

N6h = struct( ...
    'k',a{7}.k(a{7}.c<0), ...
    'c',-a{7}.c(a{7}.c<0));

% Explicit formulas from pen/paper change of variables stated in the paper for sanity check.
N5hFormula = psum( ...
    pprod(delta4,alpha,d,j,zOverTau1,rOverTheta1,theta1,tau2), ...
    pprod(delta1,gamma,d,epsilon,zOverTau1,rOverTheta1,theta1,tau2));

N6hFormula = psum( ...
    pprod(alpha,d,j,zOverTau1,rOverTheta1,theta1,tau2,gamma), ...
    pprod(alpha,d,j,zOverTau1,rOverTheta1,theta1,tau2,epsilon), ...
    pprod(alpha,d,j,zOverTau1,rOverTheta1,theta1,tau2,theta2), ...
    pprod(gamma,d,epsilon,zOverTau1,rOverTheta1,theta1,tau2,alpha), ...
    pprod(gamma,d,epsilon,zOverTau1,rOverTheta1,theta1,tau2,j), ...
    pprod(gamma,d,epsilon,zOverTau1,rOverTheta1,theta1,tau2,tau1), ...
    pprod(alpha,gamma,epsilon,j,zOverTau1,rOverTheta1,theta1,tau2));

% Verify the paper formulas against the parts computed from the QEP.
assert_poly_equal(N5h,N5hFormula, ...
    'Computed N5h disagrees with the paper formula.');

assert_poly_equal(N6h,N6hFormula, ...
    'Computed N6h disagrees with the paper formula.');

assert(numel(P5h.k)==88 && all(P5h.c>0), ...
    'P5 check failed.');

assert(numel(P6h.k)==18 && all(P6h.c>0), ...
    'P6 check failed.');

    a0=a{1}; a1=a{2}; a2=a{3}; a3=a{4};
    a4=a{5}; a5=a{6}; a6=a{7};

    % Numerators after putting every expression over its power of a0:
    %   Delta2 = D2h/a0^2,
    %   Delta3 = D3h/a0^3,
    %   Delta4 = D4h/a0^4,
    %   B1plus = B1h/a0^4,
    %   B2plus = B2h/a0^4.
    D2h = padd(pmul(a1,a2),pmul(a0,a3),-1);

    D3h = padd(pmul(a3,D2h),pmul(pprod(a1,a1),a4),-1);
    D3h = padd(D3h,pprod(a0,a1,a5));

    % Numerator of Delta3-c1*c5.
    D30h = padd(D3h,pprod(a0,a1,a5),-1);

    % Numerators of c2*Delta2-2*c1*c4 and c1*c2-2*c3.
    Lh = padd(pmul(a2,D2h),pprod(a0,a1,a4),-2);
    Ah = padd(pmul(a1,a2),pmul(a0,a3),-2);

    fprintf('  building the numerator of Delta4 ...\n');
    D4h = padd(pmul(a4,D30h),pmul(a5,Lh),-1);
    D4h = padd(D4h,pprod(a0,a0,a5,a5),-1);
    D4h = padd(D4h,pprod(a1,D2h,a6));

    fprintf('  building the numerator of B1plus ...\n');
    B1h = padd(pmul(a3,D3h),pprod(a1,a5,D2h),-1);
    B1h = padd(B1h,pprod(a1,a1,a1,a6));

    fprintf('  building the numerator of B2plus ...\n');
    B2inner = padd(Lh,pprod(a0,a0,padd(P5h,a5)));
    B2h = padd(pmul(P5h,B2inner),pprod(a1,P6h,Ah),-1);

    fprintf('  assembling Ttilde(a) ...\n');
    C = padd(pmul(P5h,D4h),pmul(P6h,B1h),-1);
    C = padd(C,pmul(N5h,B2h),-1);

    check_integer_poly(C);
    fprintf('Construction complete in %.1f seconds.\n',toc(buildTimer));
end

% =======================================================================
% Mass-action source collection and denominator-clearance transform
% =======================================================================

function [A,baseKey,NR,NT] = paper_source_matrix(C)
    delta1=1; delta3=2; delta4=3; zOverTau1=13; rOverTheta1=14;

%keeping track of the exponents of the 5 variables involved in the internal mass action change from mu' to mu
    e1=double(getexp(C.k,delta1));
    e3=double(getexp(C.k,delta3));
    e4=double(getexp(C.k,delta4));
    iXi=double(getexp(C.k,zOverTau1));
    iEta=double(getexp(C.k,rOverTheta1));

    % z/tau1=w/tau2=xi2/(delta3*(1+xi2));
    % r/theta1=u/theta2=eta2/(1+eta2).
    er=e1;
    es=e4;
    assert(all(iXi<=4) && all(iEta<=4), ...
        'Source degree exceeds the mass-action transform table.');
    keptExp=e1+e3+2-iXi;
    assert(all(er>=0) && all(es>=0) && all(keptExp>=0), ...
        'Invalid exponents in the mass-action substitution.');
    assert(all(er<=15) && all(es<=15) && all(keptExp<=15), ...
        'Packed-exponent overflow in the mass-action substitution.');

    NR=max(er); NT=max(es);
    assert(NR<=5 && NT<=5,'Delta1 or delta4 degree exceeds the transform table.');

    base=C.k;
    for v=[delta1 delta3 delta4 zOverTau1 rOverTheta1]
        base=clearexp(base,v);
    end
    base=setexp(base,delta3,keptExp);

    % Column order: z/tau1, r/theta1, delta1, delta4.
    src=1+iXi+5*iEta+25*er+25*(NR+1)*es;
    nsrc=25*(NR+1)*(NT+1);

    [base,ord]=sort(base);
    src=src(ord);
    val=C.c(ord);
    first=[true;base(2:end)~=base(1:end-1)];
    row=cumsum(first);
    baseKey=base(first);
    nrows=double(row(end));

    certify_sparse_accumulation(row,src,val,nrows,nsrc, ...
        'mass-action source assembly');
    A=sparse(row,src,val,nrows,nsrc);
end


function T = mass_action_clearance_transform()
    % Output powers: xi1=0:5, xi2=0:9, eta1=0:5, eta2=0:9.
    I=[]; J=[]; V=[];

    for d=0:5
    for c=0:5
    for b=0:4
    for a=0:4
        row=1+a+5*b+25*c+150*d;

        n1=5-c;
        n2=4-a+c;
        n3=5-d;
        n4=4-b+d;
        assert(min([n1 n2 n3 n4])>=0);

        for r1=0:n1
        for r2=0:n2
        for r3=0:n3
        for r4=0:n4
            eXi1=c+r1;
            eXi2=5+a-c+r2;
            eEta1=d+r3;
            eEta2=5+b-d+r4;
            col=ratio_index(eXi1,eXi2,eEta1,eEta2);

            val=nchoosek(n1,r1)*nchoosek(n2,r2)* ...
                nchoosek(n3,r3)*nchoosek(n4,r4);

            I(end+1,1)=row; %#ok<AGROW>
            J(end+1,1)=col; %#ok<AGROW>
            V(end+1,1)=val; %#ok<AGROW>
        end
        end
        end
        end
    end
    end
    end
    end

    certify_sparse_accumulation(I,J,V,900,3600, ...
        'mass-action transform assembly');
    T=sparse(I,J,V,900,3600);
end

function idx = ratio_index(xi1,xi2,eta1,eta2)
    assert(xi1>=0 && xi1<=5 && xi2>=0 && xi2<=9 && ...
           eta1>=0 && eta1<=5 && eta2>=0 && eta2<=9);
    idx=1+xi1+6*xi2+60*eta1+360*eta2;
end

% =======================================================================
% Structural polynomial E and arm symmetry
% =======================================================================

function E = paper_E()
    persistent saved
    if ~isempty(saved), E=saved; return; end

    alpha=pvar(4); tau1=pvar(5); theta1=pvar(6); tau2=pvar(7);
    gamma=pvar(8); theta2=pvar(9); j=pvar(10); epsilon=pvar(12);
    X=padd(alpha,tau1);
    U=padd(X,theta1);
    Y=padd(gamma,theta2);
    Z=padd(Y,tau2);

    F=psum( ...
        pscale(pprod(X,X,X),4), ...
        pscale(pprod(X,X,theta1),14), ...
        pscale(pprod(X,theta1,theta1),12), ...
        pscale(pprod(theta1,theta1,theta1),2));

    Q11=psum(pscale(pprod(U,U),2), ...
             pscale(pprod(Z,Z),2), ...
             pscale(pprod(U,Z),-2));

    Q10=psum( ...
        pscale(pprod(Z,Z,Z),2), ...
        pscale(pprod(tau2,Y,Z),4), ...
        pscale(pprod(U,Z,Z),-1), ...
        pscale(pprod(U,tau2,tau2),2), ...
        pscale(pprod(U,U,padd(Z,tau2)),6));

    leading=padd(pscale(X,9),pscale(theta1,16));
    middle=psum(pprod(theta1,theta1), ...
                pscale(pprod(theta1,X),-4), ...
                pscale(pprod(X,X),-3));
    % The quadratic variable is Z=gamma+theta2+tau2 (not Y=gamma+theta2).
    Q01=psum(pprod(leading,Z,Z),pprod(middle,Z),F);

    E=psum(pmul(j,Q10),pmul(epsilon,Q01),pprod(j,epsilon,Q11));
    assert(numel(E.k)==113 && nnz(E.c<0)==39, ...
        'Main-envelope support check failed.');

    % Q11 = U^2+Z^2+(U-Z)^2 > 0.
    UmZ=padd(U,Z,-1);
    rhs11=psum(pprod(U,U),pprod(Z,Z),pprod(UmZ,UmZ));
    assert_poly_equal(Q11,rhs11,'Q11 identity failed.');

    % Q10 is a sum of positive terms.  The only quadratic needed is
    % 2*Z^2-U*Z+6*U^2, and
    % 8*(2*Z^2-U*Z+6*U^2)=(4*Z-U)^2+47*U^2.
    quad=psum(pscale(pprod(Z,Z),2),pscale(pprod(U,Z),-1), ...
              pscale(pprod(U,U),6));
    rhs10=psum(pprod(Z,quad),pscale(pprod(tau2,Y,Z),4), ...
               pscale(pprod(U,tau2,tau2),2),pscale(pprod(U,U,tau2),6));
    assert_poly_equal(Q10,rhs10,'Q10 identity failed.');
    fourZmU=padd(pscale(Z,4),U,-1);
    quadSquare=psum(pprod(fourZmU,fourZmU),pscale(pprod(U,U),47));
    assert_poly_equal(pscale(quad,8),quadSquare,'Q10 square identity failed.');

    % Q01=leading*Z^2+middle*Z+F.  Its discriminant is negative because
    % 4*leading*F-middle^2 equals the following positive polynomial.
    discGap=padd(pscale(pprod(leading,F),4),pprod(middle,middle),-1);
    discPositive=psum( ...
        pscale(pprod(X,X,X,X),135), ...
        pscale(pprod(theta1,X,X,X),736), ...
        pscale(pprod(theta1,theta1,X,X),1318), ...
        pscale(pprod(theta1,theta1,theta1,X),848), ...
        pscale(pprod(theta1,theta1,theta1,theta1),127));
    assert_poly_equal(discGap,discPositive,'Q01 discriminant identity failed.');
    assert(all(F.c>0) && all(leading.c>0) && all(discPositive.c>0), ...
        'Q01 positivity check failed.');

    saved=E;
end

function P = swap_arms(P)
    % alpha<->gamma, tau1<->theta2, theta1<->tau2, j<->epsilon.
    exponents=decodekeys(P.k);
    pairs=[4 8;5 9;6 7;10 12];
    for r=1:size(pairs,1)
        i=pairs(r,1);
        k=pairs(r,2);
        tmp=exponents(:,i);
        exponents(:,i)=exponents(:,k);
        exponents(:,k)=tmp;
    end
    P=pnormalize(packrows(exponents),P.c);
end



% =======================================================================
% Gamma1, Gamma2 and the certificate matrix
% =======================================================================

function [Gamma1,Gamma2] = gamma_ratio_factors()
    % Exponent columns are [xi1 xi2 eta1 eta2].

    e=[]; c=[];
    % xi2^6*eta1*eta2^5*(1+xi2)^3*(1+eta2)^4.
    for a=0:3
        for b=0:4
            e(end+1,:)=[0,6+a,1,5+b]; %#ok<AGROW>
            c(end+1,1)=nchoosek(3,a)*nchoosek(4,b); %#ok<AGROW>
        end
    end
    % xi2^9*eta1*eta2^8*(eta1*(1+eta2)+xi1*eta2).
    e=[e; 0 9 2 8; 0 9 2 9; 1 9 1 9];
    c=[c;1;1;1];
    Gamma1=ratio_canon(e,c);

    e=[]; c=[];
    % xi1*xi2^5*eta2^6*(1+xi2)^4*(1+eta2)^3.
    for a=0:4
        for b=0:3
            e(end+1,:)=[1,5+a,0,6+b]; %#ok<AGROW>
            c(end+1,1)=nchoosek(4,a)*nchoosek(3,b); %#ok<AGROW>
        end
    end
    % xi1*xi2^8*eta2^9*(xi1*(1+xi2)+xi2*eta1).
    e=[e; 2 8 0 9; 2 9 0 9; 1 9 1 9];
    c=[c;1;1;1];
    Gamma2=ratio_canon(e,c);
end

function W = ratio_canon(e,c)
    [eu,~,g]=unique(e,'rows','sorted');
    cu=accumarray(g,c,[],@sum);
    use=cu~=0;
    W=struct('e',eu(use,:),'c',cu(use));
    assert(all(W.c==fix(W.c)) && all(abs(W.c)<flintmax));
end

function S = certificate_matrix(baseKey,Gamma1E,Gamma1Ratios, ...
        Gamma2E,Gamma2Ratios)
    [I1,J1,V1]=separable_entries(baseKey,Gamma1E,Gamma1Ratios);
    [I2,J2,V2]=separable_entries(baseKey,Gamma2E,Gamma2Ratios);

    I=[I1;I2]; J=[J1;J2]; V=[V1;V2];
    certify_sparse_accumulation(I,J,V,numel(baseKey),3600, ...
        'Gamma1*E+Gamma2*Esigma assembly');
    S=sparse(I,J,V,numel(baseKey),3600);
    check_integer_matrix(S,'certificate matrix');
end

function [I,J,V] = separable_entries(baseKey,KE,W)
    [tf,row]=ismember(KE.k,baseKey);
    assert(all(tf),'A Gamma*E source monomial is absent from TtildeStar.');

    nw=numel(W.c); ne=numel(KE.c);
    I=zeros(ne*nw,1); J=zeros(ne*nw,1); V=zeros(ne*nw,1);
    z=0;
    for r=1:ne
        rr=z+(1:nw);
        I(rr)=row(r);
        J(rr)=1+W.e(:,1)+6*W.e(:,2)+60*W.e(:,3)+360*W.e(:,4);
        certify_outer_products(KE.c(r),W.c,'Gamma*E coefficient products');
        V(rr)=KE.c(r)*W.c;
        z=z+nw;
    end
end




% =======================================================================
% Exact denominator-clearance multiplicities
% =======================================================================

function C = init_clearance_minimality()
    C.names = {'delta3','xi2','eta2','1+xi1','1+eta1','1+xi2','1+eta2'};
    C.chosen = [2 5 5 5 5 4 4];
    C.multiplicity = inf(1,7);
    C.binomialWitness = false(1,7);
    C.row = zeros(1,7);
    C.value = zeros(1,7);
end

function [EvalXi1,EvalEta1,EvalXi2,EvalEta2,xi2Exp,eta2Exp] = ...
        clearance_test_operators()
    n=3600;
    cols=(1:n)';
    q=cols-1;
    xi1=mod(q,6); q=floor(q/6);
    xi2Exp=mod(q,10); q=floor(q/10);
    eta1=mod(q,6); q=floor(q/6);
    eta2Exp=mod(q,10);

    EvalXi1=sparse(cols,1+xi2Exp+10*eta1+60*eta2Exp, ...
        (-1).^xi1,n,600);
    EvalEta1=sparse(cols,1+xi1+6*xi2Exp+60*eta2Exp, ...
        (-1).^eta1,n,600);
    EvalXi2=sparse(cols,1+xi1+6*eta1+36*eta2Exp, ...
        (-1).^xi2Exp,n,360);
    EvalEta2=sparse(cols,1+xi1+6*xi2Exp+60*eta1, ...
        (-1).^eta2Exp,n,360);

    check_integer_matrix(EvalXi1,'xi1=-1 evaluation map');
    check_integer_matrix(EvalEta1,'eta1=-1 evaluation map');
    check_integer_matrix(EvalXi2,'xi2=-1 evaluation map');
    check_integer_matrix(EvalEta2,'eta2=-1 evaluation map');
end

function C = update_clearance_minimality(C,TstarBlk,blockKeys,globalRows, ...
        EvalXi1,EvalEta1,EvalXi2,EvalEta2,xi2Exp,eta2Exp)

    rowNonzero=full(any(TstarBlk~=0,2));

    rows=find(rowNonzero);
    if ~isempty(rows)
        delta3Exp=double(getexp(blockKeys,2));
        [m,z]=min(delta3Exp(rows));
        if m<C.multiplicity(1)
            r=rows(z); [~,~,v]=find(TstarBlk(r,:),1);
            C.multiplicity(1)=m; C.row(1)=globalRows(r); C.value(1)=v;
        end
    end

    cols=find(full(any(TstarBlk~=0,1)));
    if ~isempty(cols)
        [m,z]=min(xi2Exp(cols));
        if m<C.multiplicity(2)
            c=cols(z); [r,~,v]=find(TstarBlk(:,c),1);
            C.multiplicity(2)=m; C.row(2)=globalRows(r); C.value(2)=v;
        end
        [m,z]=min(eta2Exp(cols));
        if m<C.multiplicity(3)
            c=cols(z); [r,~,v]=find(TstarBlk(:,c),1);
            C.multiplicity(3)=m; C.row(3)=globalRows(r); C.value(3)=v;
        end
    end

    if ~C.binomialWitness(4)
        Q=exact_matrix_product(TstarBlk,EvalXi1,'clearance: xi1=-1');
        [r,~,v]=find(Q,1);
        if ~isempty(r)
            C.binomialWitness(4)=true; C.multiplicity(4)=0;
            C.row(4)=globalRows(r); C.value(4)=v;
        end
    end
    if ~C.binomialWitness(5)
        Q=exact_matrix_product(TstarBlk,EvalEta1,'clearance: eta1=-1');
        [r,~,v]=find(Q,1);
        if ~isempty(r)
            C.binomialWitness(5)=true; C.multiplicity(5)=0;
            C.row(5)=globalRows(r); C.value(5)=v;
        end
    end
    if ~C.binomialWitness(6)
        Q=exact_matrix_product(TstarBlk,EvalXi2,'clearance: xi2=-1');
        [r,~,v]=find(Q,1);
        if ~isempty(r)
            C.binomialWitness(6)=true; C.multiplicity(6)=0;
            C.row(6)=globalRows(r); C.value(6)=v;
        end
    end
    if ~C.binomialWitness(7)
        Q=exact_matrix_product(TstarBlk,EvalEta2,'clearance: eta2=-1');
        [r,~,v]=find(Q,1);
        if ~isempty(r)
            C.binomialWitness(7)=true; C.multiplicity(7)=0;
            C.row(7)=globalRows(r); C.value(7)=v;
        end
    end
end

function C = finalize_clearance_minimality(C)
    assert(all(isfinite(C.multiplicity(1:3))), ...
        'A monomial clearance multiplicity was not determined.');
    assert(all(C.binomialWitness(4:7)), ...
        'A binomial nondivisibility witness was not found.');
    assert(all(C.multiplicity==fix(C.multiplicity)) && ...
           all(C.multiplicity>=0), ...
        'Invalid numerator multiplicity.');

    C.reduced=C.chosen-C.multiplicity;
    assert(all(C.reduced>=0),'The chosen factor does not clear the denominator.');
    assert(isequal(C.reduced,[2 4 4 5 5 4 4]), ...
        'Unexpected reduced denominator clearance.');

    % The fifth power of det(D) shifts these seven factor exponents by
    % [10,-5,-5,-5,-5,5,5].  A displayed multiplier of Ttilde must have
    % no negative exponent after this shift.
    detD5Shift=[10 -5 -5 -5 -5 5 5];
    laurentFreeLower=max(zeros(1,7),-detD5Shift);
    C.minimalLaurentFree=max(C.reduced,laurentFreeLower);
    C.chosenIsMinimalLaurentFree=isequal(C.chosen,C.minimalLaurentFree);
    assert(C.chosenIsMinimalLaurentFree, ...
        'The paper clearance is not the minimal Laurent-free clearance.');
end

function report_clearance_minimality(C)
    fprintf('\n============================================================\n');
    fprintf('DENOMINATOR CLEARANCE VERIFIED\n');
    fprintf('============================================================\n');
    fprintf('Factor order:\n');
    fprintf('  delta3, xi2, eta2, 1+xi1, 1+eta1, 1+xi2, 1+eta2\n');
    fprintf('Chosen exponents                 : [%s]\n',num2str(C.chosen));
    fprintf('Numerator multiplicities         : [%s]\n',num2str(C.multiplicity));
    fprintf('Reduced denominator exponents    : [%s]\n',num2str(C.reduced));
    fprintf('Minimal Laurent-free exponents   : [%s]\n',num2str(C.minimalLaurentFree));
    fprintf('Exact witnesses:\n');
    for i=1:numel(C.names)
        fprintf('  %-7s multiplicity %d  (group %d, witness %.0f)\n', ...
            C.names{i},C.multiplicity(i),C.row(i),C.value(i));
    end
    fprintf('The chosen 2555544 clearance is minimal among Laurent-free choices.\n');
end




% =======================================================================
% Sparse polynomial arithmetic
% =======================================================================

function assert_poly_equal(A,B,message)
    D=padd(A,B,-1);
    assert(isempty(D.k),message);
end

function P = pzero()
    P=struct('k',zeros(0,1,'uint64'),'c',zeros(0,1));
end

function P = pconst(c)
    if c==0, P=pzero(); else, P=struct('k',uint64(0),'c',double(c)); end
end

function P = pvar(i)
    P=struct('k',bitshift(uint64(1),4*(i-1)),'c',1);
end

function P = pscale(P,s)
    if s==0, P=pzero(); return; end
    P.c=exact_scale_vector(P.c,s,'pscale');
end

function C = padd(A,B,scaleB)
    if nargin<3, scaleB=1; end
    if isempty(A.k), C=pscale(B,scaleB); return; end
    if isempty(B.k), C=A; return; end
    scaledB=exact_scale_vector(B.c,scaleB,'padd scaling');
    C=pnormalize([A.k;B.k],[A.c;scaledB]);
end

function C = psum(varargin)
    C=pzero();
    for z=1:nargin, C=padd(C,varargin{z}); end
end

function C = pprod(varargin)
    C=pconst(1);
    sizes=cellfun(@(x)numel(x.k),varargin);
    [~,ord]=sort(sizes);
    for z=ord, C=pmul(C,varargin{z}); end
end

function C = pmul(A,B)
    if isempty(A.k)||isempty(B.k), C=pzero(); return; end

    % Packed keys use one 4-bit field per exponent.  Direct uint64 key
    % addition is valid if and only if no individual exponent sum exceeds
    % 15; otherwise a carry would silently corrupt the next field.
    assert_pack_add_safe(A.k,B.k,'pmul');

    if numel(A.k)>numel(B.k), T=A; A=B; B=T; end

    if numel(A.k)==1
        C.k=B.k+A.k;
        C.c=exact_scale_vector(B.c,A.c,'pmul scalar product');
        return
    end

    pairBudget=1.5e6;
    nb=numel(B.k);
    block=max(1,floor(pairBudget/nb));
    nchunk=ceil(numel(A.k)/block);
    chunks=cell(nchunk,1);
    z=0;
    for first=1:block:numel(A.k)
        z=z+1;
        ii=first:min(first+block-1,numel(A.k));
        kk=A.k(ii)+B.k.';
        certify_outer_products(A.c(ii),B.c,'pmul coefficient products');
        cc=A.c(ii)*B.c.';
        chunks{z}=pnormalize(kk(:),cc(:));
    end

    while numel(chunks)>1
        next=cell(ceil(numel(chunks)/2),1);
        w=0;
        for z=1:2:numel(chunks)
            w=w+1;
            if z==numel(chunks), next{w}=chunks{z};
            else, next{w}=padd(chunks{z},chunks{z+1});
            end
        end
        chunks=next;
    end
    C=chunks{1};
end

function P = pnormalize(k,c)
    use=c~=0;
    k=k(use); c=c(use);
    if isempty(k), P=pzero(); return; end
    [k,ord]=sort(k); c=c(ord);
    first=[true;k(2:end)~=k(1:end-1)];
    group=cumsum(first);
    certify_group_accumulation(group,c,'pnormalize coefficient collection');
    c=accumarray(group,c,[],@sum);
    k=k(first);
    use=c~=0;
    P=struct('k',k(use),'c',c(use));
    check_integer_poly(P);
end

function check_integer_poly(P)
    assert(all(P.c==fix(P.c)) && all(abs(P.c)<flintmax), ...
        'A coefficient is no longer an exactly represented integer.');
end

function C = pmatmul(A,B)
    [m,k]=size(A); [k2,n]=size(B); assert(k==k2);
    C=zero_matrix(m,n);
    for i=1:m
        for j=1:n
            for z=1:k
                C{i,j}=padd(C{i,j},pmul(A{i,z},B{z,j}));
            end
        end
    end
end

function M = zero_matrix(m,n)
    M=cell(m,n);
    for i=1:m, for j=1:n, M{i,j}=pzero(); end, end
end

function C = sconv(A,B)
    C=cell(1,numel(A)+numel(B)-1);
    for z=1:numel(C), C{z}=pzero(); end
    for i=1:numel(A)
        for j=1:numel(B)
            C{i+j-1}=padd(C{i+j-1},pmul(A{i},B{j}));
        end
    end
end

function dets = pdet4_s(Qs)
    dets=cell(1,9);
    for z=1:9, dets{z}=pzero(); end
    P=perms(1:4);
    for r=1:size(P,1)
        term={pconst(permutation_sign(P(r,:)))};
        for i=1:4, term=sconv(term,Qs{i,P(r,i)}); end
        for z=1:numel(term), dets{z}=padd(dets{z},term{z}); end
    end
end

function s = permutation_sign(p)
    inv=0;
    for i=1:numel(p), for j=i+1:numel(p), inv=inv+(p(i)>p(j)); end, end
    s=1-2*mod(inv,2);
end

% =======================================================================
% Exact-integer verification helpers
% =======================================================================

function S = init_exact_stats()
    S.maxScalarProduct = 0;
    S.maxPairProduct = 0;
    S.maxAbsAccumulation = 0;
    S.maxSparseAccumulation = 0;
    S.maxMatrixAccumulation = 0;
end

 function report_exact_stats()
    global EXACT_STATS VERIFY_EXACT_ARITHMETIC
    if ~VERIFY_EXACT_ARITHMETIC, return; end

    M=max([EXACT_STATS.maxScalarProduct,EXACT_STATS.maxPairProduct, ...
           EXACT_STATS.maxAbsAccumulation,EXACT_STATS.maxSparseAccumulation, ...
           EXACT_STATS.maxMatrixAccumulation]);

    assert(M<flintmax, ...
        'Exact-integer safety bound exceeded flintmax.');
end

function y = exact_scale_vector(x,s,label)
    global EXACT_STATS VERIFY_EXACT_ARITHMETIC
    if isempty(x), y=x; return; end
    if VERIFY_EXACT_ARITHMETIC
        assert(s==fix(s) && abs(s)<flintmax,[label ': noninteger/oversized scalar.']);
        ax=max(abs(x)); as=abs(s);
        if as~=0
            assert(ax < flintmax/as,[label ': scalar product may exceed flintmax.']);
        end
    end
    y=x*s;
    if VERIFY_EXACT_ARITHMETIC
        assert(all(y==fix(y)) && all(abs(y)<flintmax), ...
            [label ': loss of exact integer arithmetic.']);
        EXACT_STATS.maxScalarProduct=max(EXACT_STATS.maxScalarProduct,max(abs(y)));
    end
end

function certify_outer_products(a,b,label)
    global EXACT_STATS VERIFY_EXACT_ARITHMETIC
    if ~VERIFY_EXACT_ARITHMETIC || isempty(a) || isempty(b), return; end
    ma=max(abs(a)); mb=max(abs(b));
    if mb~=0
        assert(ma < flintmax/mb,[label ': coefficient product may exceed flintmax.']);
    end
    EXACT_STATS.maxPairProduct=max(EXACT_STATS.maxPairProduct,ma*mb);
end

function certify_group_accumulation(group,c,label)
    global EXACT_STATS VERIFY_EXACT_ARITHMETIC
    if ~VERIFY_EXACT_ARITHMETIC || isempty(c), return; end
    % Since every summand is an exactly represented integer, the positive
    % absolute sum bounds every possible signed partial sum, independent of
    % the order used internally by accumarray.
    absSum=accumarray(group,abs(c),[],@sum);
    m=max(absSum);
    assert(m<flintmax,[label ': an absolute coefficient accumulation reaches flintmax.']);
    EXACT_STATS.maxAbsAccumulation=max(EXACT_STATS.maxAbsAccumulation,m);
end

function certify_sparse_accumulation(i,j,v,m,n,label)
    global EXACT_STATS VERIFY_EXACT_ARITHMETIC
    if ~VERIFY_EXACT_ARITHMETIC || isempty(v), return; end
    % sparse(i,j,v,...) may combine duplicate entries.  Summing absolute
    % values gives an order-independent bound for every signed partial sum.
    Aabs=sparse(i,j,abs(v),m,n);
    nz=nonzeros(Aabs);
    if isempty(nz), mx=0; else, mx=max(nz); end
    assert(mx<flintmax,[label ': sparse accumulation reaches flintmax.']);
    EXACT_STATS.maxSparseAccumulation=max(EXACT_STATS.maxSparseAccumulation,mx);
end

function C = exact_matrix_product(A,B,label)
    global EXACT_STATS VERIFY_EXACT_ARITHMETIC
    if VERIFY_EXACT_ARITHMETIC
        av=nonzeros(A); bv=nonzeros(B);
        if ~isempty(av) && ~isempty(bv)
            ma=max(abs(av)); mb=max(abs(bv));
            assert(ma < flintmax/mb,[label ': matrix entry product may exceed flintmax.']);
        end
        % Positive absolute matrix product bounds every partial sum in the
        % signed product.  Because all contributions are nonnegative, if the
        % final bound is below flintmax then every partial sum is too.
        Bound=abs(A)*abs(B);
        bv2=nonzeros(Bound);
        if isempty(bv2), mx=0; else, mx=max(bv2); end
        assert(mx<flintmax,[label ': matrix accumulation reaches flintmax.']);
        EXACT_STATS.maxMatrixAccumulation=max(EXACT_STATS.maxMatrixAccumulation,mx);
    end
    C=A*B;
    if VERIFY_EXACT_ARITHMETIC
        cv=nonzeros(C);
        assert(all(cv==fix(cv)) && all(abs(cv)<flintmax), ...
            [label ': loss of exact integer arithmetic.']);
    end
end

function check_integer_matrix(A,label)
    global VERIFY_EXACT_ARITHMETIC
    if ~VERIFY_EXACT_ARITHMETIC, return; end
    v=nonzeros(A);
    assert(all(v==fix(v)) && all(abs(v)<flintmax), ...
        [label ': matrix is not exactly integer-valued below flintmax.']);
end

function certify_subtraction(A,B,label)
    global EXACT_STATS VERIFY_EXACT_ARITHMETIC
    if ~VERIFY_EXACT_ARITHMETIC, return; end
    Bound=abs(A)+abs(B);
    v=nonzeros(Bound);
    if isempty(v), mx=0; else, mx=max(v); end
    assert(mx<flintmax,[label ': subtraction may exceed flintmax.']);
    EXACT_STATS.maxMatrixAccumulation=max(EXACT_STATS.maxMatrixAccumulation,mx);
end

% =======================================================================
% Packed exponent utilities
% =======================================================================

function e = getexp(k,i)
    assert(isscalar(i) && i==fix(i) && i>=1 && i<=16, ...
        'getexp: packed-field index must be an integer in 1:16.');
    e=uint8(bitand(bitshift(k,-4*(i-1)),uint64(15)));
end

function k = clearexp(k,i)
    assert(isscalar(i) && i==fix(i) && i>=1 && i<=16, ...
        'clearexp: packed-field index must be an integer in 1:16.');
    mask=bitshift(uint64(15),4*(i-1));
    k=bitand(k,bitcmp(mask));
end

function k = setexp(k,i,e)
    assert(isscalar(i) && i==fix(i) && i>=1 && i<=16, ...
        'setexp: packed-field index must be an integer in 1:16.');
    ed=double(e);
    assert(all(isfinite(ed(:))) && all(ed(:)==fix(ed(:))) && ...
           all(ed(:)>=0) && all(ed(:)<=15), ...
        'setexp: every exponent must be an integer in 0:15.');
    k=clearexp(k,i);
    k=bitor(k,bitshift(uint64(ed),4*(i-1)));
end

function E = decodekeys(k)
    E=zeros(numel(k),16);
    for i=1:16, E(:,i)=double(getexp(k,i)); end
end

function k = packrow(e)
    e=double(e(:).');
    assert(numel(e)==16,'packrow: expected exactly 16 exponents.');
    assert(all(isfinite(e)) && all(e==fix(e)) && all(e>=0) && all(e<=15), ...
        'packrow: every exponent must be an integer in 0:15.');
    k=uint64(0);
    for i=1:16
        % bitor is deliberate: the validated 4-bit fields are disjoint,
        % so packing itself cannot create an arithmetic carry.
        k=bitor(k,bitshift(uint64(e(i)),4*(i-1)));
    end
end

function k = packrows(E)
    E=double(E);
    assert(size(E,2)==16,'packrows: expected exactly 16 exponent columns.');
    assert(all(isfinite(E(:))) && all(E(:)==fix(E(:))) && ...
           all(E(:)>=0) && all(E(:)<=15), ...
        'packrows: every exponent must be an integer in 0:15.');
    k=zeros(size(E,1),1,'uint64');
    for i=1:16
        k=bitor(k,bitshift(uint64(E(:,i)),4*(i-1)));
    end
end

function assert_pack_add_safe(kA,kB,label)
    % A packed uint64 key contains 16 independent 4-bit exponent fields.
    % For nonnegative exponents, the largest possible exponent produced by
    % multiplying any monomial from A with any monomial from B is
    % max_A e_i + max_B e_i.  Requiring this to be <=15 for every field
    % is therefore necessary and sufficient to rule out nibble carries in
    % the direct packed-key additions used by pmul.
    if nargin<3, label='packed-key addition'; end
    for i=1:16
        maxA=max(double(getexp(kA,i)));
        maxB=max(double(getexp(kB,i)));
        assert(maxA+maxB<=15, ...
            sprintf('%s: exponent overflow in packed field %d (%g+%g>15).', ...
                    label,i,maxA,maxB));
    end
end

function s = comma_integer(n)
    s=sprintf('%.0f',n);
    first=mod(numel(s),3); if first==0, first=3; end
    parts={s(1:first)};
    for k=first+1:3:numel(s), parts{end+1}=s(k:k+2); end %#ok<AGROW>
    s=strjoin(parts,',');
end
