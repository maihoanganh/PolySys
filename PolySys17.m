function PolySys17(x,g,h,k)

while 1>0

format long



% Find a few feasible real points x* in semialgebraic set S(g,h), where
%  g = {g0,g1,...,gm} are inequalities with g0=1 and h = {h1,...,hm} are equalities.
% 
% Adding sphere inequalities method:
% 1) Let a0,...,an in R^n such that a1-a0,...,an-a0 are linear independent
% 2) Solve L = min {|x|^2: x in S(g,h)} by SDP relaxation
%          Lk = min Ly(theta^k(|x|^2 + eps theta^2))
%               s.t. y in R^(2(k+2))
%                    M(k+2-u_j)(gj y)>=0
%                    M(k+2-w_j)(hj y)=0
%                    Ly(theta^k)=1
% where u_j=ceil(deg(g_j)/2) and w_j=ceil(deg(h_j)/2). Set g_m = Lk - |x|^2.
% 3) Solve xit = min{|x-at|^2: x in S(g cup {xij-|x-aj|^2: j=0,...,t-1},h)},
% t=0,...,n by numerical scheme of Lasserre's hierarchies
%          omegak_t = min Ly(|x-at|^2)
%                   s.t. y in R^(2(k+2))
%                        M(k+2-u_j)(gj y)>=0
%                        M(k+1)((omegaj-|x-aj|^2) y)>=0, j=0,...,t-1
%                        M(k+2-w_j)(hj y)=0              
%                        Ly(1)=1
% for t=0,...,n.
% 4) Check rank condition of each problem of optimal value omegat t in
% {0,...,n} to extract supp(mu) if solution y has representing atomic
% measure mu. In this case, x* in supp(mu).
% 5) Solve non-singular linear system Ax*=b where A=(a1-a0 ... an-a0) and
% b= (-(|aj|^2 - |a0|^2 + omegaj - omega0)/2)_{j=1,...,n}.


tic


n=length(x);

m=length(g);

l=length(h);




%degree
w=[];
for j=1:m
    w=[w;ceil(.5*degree(g(j)))];
end


%degree
u=[];
for j=1:l
    u=[u;ceil(.5*degree(h(j)))];
end



%small parameter for omega0_upper
epsilon=1e-2;


%to check zero eigenvalues
tol=1e-2;

% for pivot
TOL=1e-3;



lw=[];
for j=1:m
    lw=[lw;nchoosek(n+2+k-w(j),n)];
end

lu=[];
for j=1:l
    lu=[lu;nchoosek(n+4+2*k-2*u(j),n)];
end

%vector of monomial
v=monolist(x,2+k);

v2=monolist(x,4+2*k);






%finding upper bound of radii r0


%Gram matrix of sos
G0=sdpvar(length(v),length(v));%vaiable of sdp (*)

%building the element in quadratic module
quadratic_module=v'*G0*v;
for j=1:m
  Q{j}=sdpvar(lw(j),lw(j));
  quadratic_module=quadratic_module+v(1:lw(j))'*Q{j}*v(1:lw(j))*g(j);
end

for j=1:l
  c{j}=sdpvar(lu(j),1);
  quadratic_module=quadratic_module+v2(1:lu(j))'*c{j}*h(j);
end

theta=1+x'*x;
lambda=sdpvar(1);

F=[coefficients(theta^k*(x'*x-lambda+epsilon*theta^2)-quadratic_module,x)==0,G0>=0];

for j=1:m
    F=[F,Q{j}>=0];
end

ops = sdpsettings('solver','mosek','verbose',0);

diagnostics=optimize(F,-lambda,ops);

lambda=double(lambda);

omega0_upper=lambda;

  
fprintf('omega0_upper          problem       \n');
fprintf('%0f                 %0d\n',omega0_upper,diagnostics.problem);
fprintf('\n');

if diagnostics.problem~=0 && diagnostics.problem~=4
    fprintf('Increase order k!\n');
    break
end
    
    
%centers of sequence of balls 
a0=zeros(n,1);


k_min=max([w;u;1]);

lw=[];
for j=1:m
    lw=[lw;nchoosek(n+k_min+k-w(j),n)];
end

lu=[];
for j=1:l
    lu=[lu;nchoosek(n+2*k_min+2*k-2*u(j),n)];
end
lu_sphere=nchoosek(n+2*k_min+2*k-2,n);
%vector of monomial
v=monolist(x,k_min+k);

v2=monolist(x,2*k_min+2*k);

u0=monolist(x,k+k_min-1);

%finding radii r0


%Gram matrix of sos
G0=sdpvar(length(v),length(v));%vaiable of sdp (*)

Q0=sdpvar(length(u0));
 
lambda=sdpvar(1);

%building the element in quadratic module
quadratic_module=v'*G0*v+(omega0_upper-(x-a0)'*(x-a0))*u0'*Q0*u0;
for j=1:m
    Q{j}=sdpvar(lw(j),lw(j));
    quadratic_module=quadratic_module+v(1:lw(j))'*Q{j}*v(1:lw(j))*g(j);
end

for j=1:l
    c{j}=sdpvar(lu(j),1);
    quadratic_module=quadratic_module+v2(1:lu(j))'*c{j}*h(j);
end



F=[coefficients((x-a0)'*(x-a0)-lambda-quadratic_module,x)==0,G0>=0,Q0>=0];

for j=1:m
    F=[F,Q{j}>=0];
end

  
  
diagnostics=optimize(F,-lambda,ops);

lambda=double(lambda);

omega0=lambda;

M = dual(F(2));

dim_subM=nchoosek(n+k,k);

subM=M(1:dim_subM,1:dim_subM);


rank_M=rank(M,tol);

rank_subM=rank(subM,tol);

fprintf('omega0          problem       rank of moment matrix      rank of moment submatrix\n');
fprintf('%0f         %0d              %0d                        %0d\n',omega0,diagnostics.problem,rank_M,rank_subM);
fprintf('\n');

if diagnostics.problem~=0 && diagnostics.problem~=4
    fprintf('Increase order k!\n');
    break
end

if rank_subM==rank_M
    [V,~] = eig(double(G0));
    r=rank_M;
    V=V(:,1:r);
    V=V';
    [U,pivot] = rref(V,TOL);
    U=U';
    % Figure out multiplying matrices using YALMIP code
    w = v(pivot);
    for i = 1:n
        xw = x(i)*w;
        k = [];
        for j = 1:length(xw)
            k = [k;find(ismember(xw(j),v))];           
        end
        N{i} = U(k,:);
    end
    
    

    % Create random convex combination
    rands = rand(n,1);rands = rands/sum(rands);
    M = 0;
    for i = 1:n
        M = M + rands(i)*N{i};
    end

    [L,T] = schur(M);
    % Extract solution
    for i = 1:r
        solution=zeros(n,1);
        for j = 1:n
            solution(j) =  L(:,i)'*N{j}*L(:,i);
        end
        
        solution=solution
        if m~=0 
            check_ineq=replace(g,x, solution)
        end

        if l~=0
            check_eq=replace(h,x, solution)    
        end
        toc
    end
    break
end
   













%centers of sequence of balls 
a=3*sqrt(omega0)*eye(n);

%find radius
omega=zeros(n,1);
for t=1:n
  
    %Gram matrix of sos
    G0=sdpvar(length(v),length(v));%vaiable of sdp (*)

    Q0=sdpvar(lu_sphere,1);



    lambda=sdpvar(1);

    quadratic_module=v'*G0*v+(omega0-(x-a0)'*(x-a0))*v2(1:lu_sphere)'*Q0;
    for j=1:m
        Q{j}=sdpvar(lw(j),lw(j));
        quadratic_module=quadratic_module+v(1:lw(j))'*Q{j}*v(1:lw(j))*g(j);
    end

    for j=1:l
        c{j}=sdpvar(lu(j),1);
        quadratic_module=quadratic_module+v2(1:lu(j))'*c{j}*h(j);
    end



    if t>1
        for j=1:t-1
           H{j}=sdpvar(lu_sphere,1);
           quadratic_module=quadratic_module+v2(1:lu_sphere)'*H{j}*(omega(j)-(x-a(:,j))'*(x-a(:,j)));
        end
    end

    F=[coefficients((x-a(:,t))'*(x-a(:,t))-lambda-quadratic_module,x)==0,G0>=0];
    for j=1:m
        F=[F,Q{j}>=0];
    end


    diagnostics=optimize(F,-lambda,ops);

    lambda=double(lambda);

    omega(t)=lambda;

    M = dual(F(2));

    dim_subM=nchoosek(n+k,k);

    subM=M(1:dim_subM,1:dim_subM);


    rank_M=rank(M,tol);

    rank_subM=rank(subM,tol);

    fprintf('omega%0d         problem       rank of moment matrix      rank of moment submatrix\n',t);
    fprintf('%0f         %0d              %0d                        %0d\n',omega(t),diagnostics.problem,rank_M,rank_subM);
    fprintf('\n');
    
    if diagnostics.problem~=0 && diagnostics.problem~=4
        fprintf('Increase order k!\n');
        break
    end

    if rank_subM==rank_M
        [V,~] = eig(double(G0));
        r=rank_M;
        V=V(:,1:r);
        V=V';
        [U,pivot] = rref(V,TOL);
        U=U';
        % Figure out multiplying matrices using YALMIP code
        w = v(pivot);
        for i = 1:n
            xw = x(i)*w;
            k = [];
            for j = 1:length(xw)
                k = [k;find(ismember(xw(j),v))];     
            end
            N{i} = U(k,:);
        end



        % Create random convex combination
        rands = rand(n,1);rands = rands/sum(rands);
        M = 0;
        for i = 1:n
            M = M + rands(i)*N{i};
        end

        [L,T] = schur(M);
        % Extract solution
        solution=zeros(n,1);
        for i = 1:r
            for j = 1:n
                solution(j) =  L(:,i)'*N{j}*L(:,i);
            end
            
            solution=solution
            if m~=0 
                check_ineq=replace(g,x, solution)
            end

            if l~=0
                check_eq=replace(h,x, solution)   
            end
            toc
        end
        break
    end
   

end





 
    
fprintf('Computing approximate root from sphere equations...\n');

%solve linear programming
A=[];
b=[];
for j=1:n
    A=[A;-(a(:,j)-a0)'];
    b=[b;omega(j)-omega0-norm(a(:,j))^2+norm(a0)^2];
end
b=.5*b;

solution=inv(A)*b
   
if m~=0 
    check_ineq=replace(g,x, solution)
end

if l~=0
   check_eq=replace(h,x, solution)    
end
toc
break

end