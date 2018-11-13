% calcY0
%% import case
mpc = case9();
N = size(mpc.bus,1);
b = size(mpc.branch,1);
%% A0
nodeA = [mpc.branch(:,1);mpc.branch(:,2);(N+1)*ones(N,1);(1:N)'];
branchA = [1:b,1:b,b+(1:N),b+(1:N)]';
valueA = [ones(b,1);-ones(b,1);ones(N,1);-ones(N,1)];
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
disp(nnz(Y)/numel(Y));
disp('diagonal priority');
figure();
yyaxis left;
histDiagPri(mpc);
ylabel('num');
yyaxis right;
histDiagPri(case118());
xlabel('diagElem / max(abs(otherColumnElem))');
ylabel('num');
legend('case9','','case118','');
saveas(gcf, [pwd '\meta\diagpri.png']);
disp('non-singularity');
fprintf('|det(Y)| = %s\n', abs(det(Y)));
