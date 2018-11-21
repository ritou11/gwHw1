% calcY0
%% import case
mpc = case9();
N = size(mpc.bus,1);
b = size(mpc.branch,1);
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
%% A0
nodeA = [mpc.branch(:,1);mpc.branch(:,2);(N+1)*ones(N,1);(1:N)'];
branchA = [1:b,1:b,b+(1:N),b+(1:N)]';
valueA = [ones(b,1);-ones(b,1);-ones(N,1);ones(N,1)];
A0 = sparse(nodeA, branchA, valueA);
%% yb
A = sparse(mpc.branch(:,1:2),[1:b,1:b]',ones(b,2));
y0 = A*1j*mpc.branch(:,5)/2 + mpc.bus(mpc.bus(:,1),5) + 1j*mpc.bus(mpc.bus(:,1),6);
yb = spdiags([1./(mpc.branch(:,3)+1j*mpc.branch(:,4));y0],0,b+N,b+N);
%% Y0 Y
Y0 = A0*yb*A0';
Y = Y0(1:N, 1:N);
%% makeYbus
disp('makeYbus');
[Ybus, ~, ~] = makeYbus(mpc);
disp(Ybus - Y);
%% feature of Y
disp('sparse');
fprintf('sparse degree = %.3f\n', nnz(Y)/numel(Y));
disp('diagonal priority');
f=figure();
yyaxis left;
histDiagPri(mpc);
ylabel('num');
yyaxis right;
histDiagPri(case118());
xlabel('diagElem / max(abs(otherColumnElem))');
ylabel('num');
legend('case9','','case118','');
saveas(f, [pwd '\meta\diagpri.png']);
close(f);
disp('non-singularity');
fprintf('|det(Y)| = %s\n', abs(det(Y)));
%% LDU
disp('LDU');
[L, U] = calcLU(Y);
fprintf('LU error = %s\n', norm(full(L*U - Y)));
[L, D, U] = calcLDU(Y);
fprintf('LDU error = %s\n', norm(full(L*D*U - Y)));
%% Modify Y
mpc_m = case9_modified();
mpYm = makeYbus(mpc_m);

fbn = mpc_m.branch(end, F_BUS);
tbn = mpc_m.branch(end, T_BUS);
ybr = 1./(mpc_m.branch(end, BR_R) + 1j*mpc_m.branch(end, BR_X));
bbn = mpc_m.branch(end, BR_B) / 2;
Ml = sparse([fbn tbn], [1 1], [1 -1], N, 1);
Ym = Y + Ml * ybr * Ml' + sparse([fbn tbn], [fbn tbn], [1j 1j] * bbn, N, N);
disp('modify Y');
fprintf('Ym error = %s\n', norm(full(Ym - mpYm)));
%% Modify Z
Z = inv(Y);
zaa = mpc_m.branch(end, BR_R) + 1j*mpc_m.branch(end, BR_X) + Ml' * Z * Ml;
Zm = Z - Z * Ml / zaa * Ml' * Z;
zaa = 1 / (1j * bbn) + Zm(fbn, fbn);
Zm = Zm - Zm(:, fbn) / zaa * Zm(fbn, :);
zaa = 1 / (1j * bbn) + Zm(tbn, tbn);
Zm = Zm - Zm(:, tbn) / zaa * Zm(tbn, :);
disp('modify Z');
fprintf('Zm error = %s\n', norm(full(Zm - inv(mpYm))));
%% Modify LDU - Rank1
% symmetry matrix
[L, D, U] = calcLDU(mpYm);
tic;
[Lm1, Dm1, Um1] = modifyLDUr1(D, U, Ml, -ybr);
[Lm1, Dm1, Um1] = modifyLDUr1(Dm1, Um1, sparse(fbn, 1, 1, N, 1), -1j * bbn);
[Lm1, Dm1, Um1] = modifyLDUr1(Dm1, Um1, sparse(tbn, 1, 1, N, 1), -1j * bbn);
time_r1 = toc;
disp('modify LDU - Rank 1');
fprintf('ldu-r1 error = %s\n', norm(full(Y - Lm1*Dm1*Um1)));
fprintf('time_r1 = %s\n', time_r1);
%% Modify LDU - Local Re
dY = sparse([fbn, tbn], [fbn, tbn], [-1j * bbn, -1j * bbn], N, N);
dY = dY + Ml * -ybr * Ml.';
tic;
[Lm2, Dm2, Um2] = modifyLDUlr(D, U, dY);
time_lr = toc;
disp('modify LDU - Local Re');
fprintf('ldu-lr error = %s\n', norm(full(Y - Lm2*Dm2*Um2)));
fprintf('time_lr = %s\n', time_lr);
%% Compensation method
% Compensation current
V = sparse(1:N, 1, 1, N, 1);
I = Ym * V;
deltay = sparse(1:3, 1:3, -[ybr, 1j * bbn, 1j * bbn], 3, 3);
M = sparse([fbn, tbn, fbn, tbn], [1,1,2,3], [1,-1,1,1], N, 3);
c = inv(inv(deltay) + M.' * inv(Ym) * M);
deltaI = - M * c * M.' * inv(Ym) * I;
Is = I + deltaI;
Vs = Ym \ Is;
Vsorg = Y \ I;
disp('Compensation current');
fprintf('V error = %s\n', norm(full(Vsorg - Vs)));
% steps
Eij = M.' * V;
full(Eij)
Zt = M.' * inv(Ym) * M;
full(Zt)
invc = inv(deltay) + Zt;
full(invc)
Iij = invc \ Eij;
full(Iij)
DeltaI = - M * Iij;
fprintf('physics error = %s\n',norm(full(Ym\(I+DeltaI) - Vsorg)));
