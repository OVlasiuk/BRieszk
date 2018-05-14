function cnf = f_cnfinit(N, surfF)
% F_CNFINIT
% cnf = f_cnfinit(N, surfF)
% Example of how to produce a random configuration on a complicated surface:
warning('off','optim:fsolve:NonSquareSystem')
cnf = zeros(3,N);
x0 = 5*rand(3,N)-2.5;
for i=1:N
cnf(:,i) = fsolve(surfF, x0(:,i),optimoptions('fsolve','Display','off'));
end
