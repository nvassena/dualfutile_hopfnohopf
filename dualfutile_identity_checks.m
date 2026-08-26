function dualfutile_identity_checks()
%IDENTITY_CHECKS Exact symbolic checks of three identities used in the paper.

    clc;
    timer=tic;
    syms a b c d f g h o r u v w x y z real
    j=sym('j','real');

    %% 1. Full and reduced characteristic polynomials
    G=[ ...
        -a, 0, x, 0, -h, 0, 0, w, 0;
        0, -b-d, z, y, -j, -f, 0, v, u;
        a, 0, -x-z, 0, h, 0, 0, 0, 0;
        0, d, 0, -y-r, j, 0, 0, 0, 0;
        -a, -d, x+z, y+r, -h-j, 0, 0, 0, 0;
        0, -b, 0, 0, 0, -f-g, -c, w+v, u+o;
        0, 0, 0, r, 0, -g, -c, 0, o;
        0, b, 0, 0, 0, f, 0, -v-w, 0;
        0, 0, 0, 0, 0, g, c, 0, -u-o];

    Gred=[ ...
        0, 0, -z, 0, w, 0;
        0, 0, 0, r, 0, -u;
        a, 0, -a-h-x-z, -h, 0, 0;
        -d, -d, -j, -d-j-y-r, -d, 0;
        -b, -b, 0, -b, -v-b-f-w, -f;
        0, c, 0, 0, -g, -o-c-g-u];

    fullCoefficients=charpoly(G);
    reducedCoefficients=charpoly(Gred);
    assert(all(isAlways(fullCoefficients(1:7)==reducedCoefficients)) && ...
           all(isAlways(fullCoefficients(8:10)==0)), ...
        'The characteristic polynomials of G and Gred do not agree.');
    fprintf('1. det(lambda*I-G)=lambda^3 det(lambda*I-Gred): PASS\n');

    %% 2. Hurwitz identity
    syms c1 c2 c3 c4 c5 c6 real
    H=[ ...
        c1 c3 c5 0  0  0;
        1  c2 c4 c6 0  0;
        0  c1 c3 c5 0  0;
        0  1  c2 c4 c6 0;
        0  0  c1 c3 c5 0;
        0  0  1  c2 c4 c6];

    Delta2=expand(det(H(1:2,1:2)));
    Delta3=expand(det(H(1:3,1:3)));
    Delta4=expand(det(H(1:4,1:4)));
    Delta5=expand(det(H(1:5,1:5)));

    hurwitzIdentity=expand(Delta2*Delta5 ...
        -(c5*Delta2-c1^2*c6)*Delta4+c6*Delta3^2);
    assert(isAlways(hurwitzIdentity==0),'The Hurwitz identity failed.');
fprintf('2. Delta2*Delta5=(c5*Delta2-c1^2*c6)*Delta4-c6*Delta3^2: PASS\n');
    %% 3. P5,N5,P6,N6 rearrangement
    B1plus=expand(c3*Delta3-c1*c5*Delta2+c1^3*c6);
    assert(isAlways(expand(Delta5-c5*Delta4+c6*B1plus)==0), ...
        'The formula Delta5=c5*Delta4-c6*B1plus failed.');

    syms P5 N5 P6 N6 real
    decomposition=P6*Delta5-(N6*P5-N5*P6)*Delta4 ...
        -c6*(P5*Delta4-P6*B1plus);
    decomposition=expand(subs(decomposition,[c5 c6], ...
        [P5-N5 P6-N6]));
    assert(isAlways(decomposition==0),'The P5,N5,P6,N6 identity failed.');
fprintf('3. P6*Delta5=(N6*P5-N5*P6)*Delta4+c6*(P5*Delta4-P6*B1plus): PASS\n');
    fprintf('\nAll exact identity checks passed in %.1f seconds.\n',toc(timer));
end
