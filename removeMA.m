function [w] = removeMA(X,pmlA)
    Q = zeros(3,1000);
    R = zeros(3,3);
    for j = 1:3
        v=pmlA(j,:);
        for i=1:j-1
            R(j,i)=Q(i,:)*pmlA(j,:)';
            v = v-R(j,i)*Q(i,:);
        end
        R(j,j)=norm(v);
        Q(j,:)=v/R(j,j);
    end;
    XX = X-(X*Q')*Q;
    w = fastica (XX, 'numOfIC', 2);
    while size(w) == [0,0]
        w = fastica (XX, 'numOfIC', 2);
    end
end