% Stagewise Orthogonal Matching Pursuit
% This function returns a column vector recovered using StOMP algorithm.
% length(x)=size(A,2);
% A--> measurement matrix
% y--> test vector: y=Ax
% N--> Number of iterations/stages
function x = StOMP(A,y,N)
    r=y; % initial Residue
    O=[]; % initialisation of Support Vectors
    n=size(A,2);
    x=zeros(n,1); % initialisation of x
    t=3; % threshold parameter
    sd=norm(r);
    for i=1:N
        if A*x==y
            break
        end
        c=A'*r; % Correlation
        sd=norm(r)/sqrt(n); % noise level(standard deviation)
        ind=find(abs(c)>=t*sd); % find the desired indices greater than threshold
        O=union(O,ind); % Update Support
        Ao=A(:,O); % Updated measurement matrix
        x1=Ao\y; % min ||y-Ao*x||
        r=y-Ao*x1; % ith step residual
        x(O)=x1;
    end
end
