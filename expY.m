%% Import case
mpc = case9();

N = size(mpc.bus,1);
b = size(mpc.branch,1);
%% Build Y – Step 1
A = sparse(mpc.branch(:,1:2),[1:b,1:b]',[ones(b,1),-ones(b,1)]);
yb = spdiags(1./(mpc.branch(:,3)+1j*mpc.branch(:,4)),0,b,b);
Y = A*yb*A';
display(det(Y));
%% Build Y – Step 2
A = sparse(mpc.branch(:,1:2),[1:b,1:b]',ones(b,2));
y0 = A*1j*mpc.branch(:,5)/2 + mpc.bus(mpc.bus(:,1),5)+1j*mpc.bus(mpc.bus(:,1),6);
Y = Y + spdiags(y0,0,N,N);
%% Test
Y0 = [Y,-y0;-y0.',sum(y0)];
display(det(Y0));
