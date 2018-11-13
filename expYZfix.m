%% Import case
% ???case9_modified??
mpc_m = case9_modified();
inn = mpc_m.branch(end,1);
jnn = mpc_m.branch(end,2);
%% Modify Y
yl = 1./(mpc_m.branch(end,3)+1j*mpc_m.branch(end,4));
bl = mpc_m.branch(end,5)/2;
Ml = sparse([inn,jnn],1,[1,-1],N,1);
Ym = Y + Ml*yl*Ml' + sparse([inn,jnn],[inn,jnn],[1j;1j]*bl,N,N);
%% Modify Z
% Add line reactance
zl = 1/yl + Ml'*Zi*Ml;
Zm = Zi - Zi*Ml/zl*Ml'*Zi;
% Add shunt reactance at initial
zl = -1j/bl + Zm(inn,inn);
Zm = Zm - Zm(:,inn)/zl*Zm(inn,:);
% Add shunt reactance at end
zl = -1j/bl + Zm(jnn,jnn);
Zm = Zm - Zm(:,jnn)/zl*Zm(jnn,:);
%% Test
display(max(max(abs(Ym-makeYbus(case9_modified())))));
display(max(max(abs(inv(Ym)-Zm))));