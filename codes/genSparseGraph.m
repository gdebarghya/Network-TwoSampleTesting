function A = genSparseGraph(m,model)

A = cell(m,1);

for i = 1:m
    switch model.name
        
        case 'ER'
            A1 = triu(ceil(sprand(model.n,model.n,model.p)),1);
            A{i} = A1|A1';
            
        case '2SBM'
            n1 = floor(model.n/2);
            n2 = model.n-n1;
            A11 = triu(ceil(sprand(n1,n1,model.p)),1);
            A22 = triu(ceil(sprand(n2,n2,model.p)),1);
            A12 = ceil(sprand(n1,n2,model.q));
            A1 = [A11 A12; sparse(n2,n1) A22];
            A{i} = A1|A1';
            
        case 'SBM'
            inc = floor(model.n/model.k);
            A1 = sparse(model.n,model.n);
            ni = 0;
            for ci = 1:model.k
                nj = ni;
                A1((ni+1):(ni+inc),(nj+1):(nj+inc)) = triu(ceil(sprand(inc,inc,model.p)),1);
                for cj = (ci+1):model.k
                    nj = nj+inc;
                    A1((ni+1):(ni+inc),(nj+1):(nj+inc)) = ceil(sprand(inc,inc,model.q));
                end
                ni = ni+inc;
            end
            A{i} = A1|A1';
            
        case 'IER'
            EA = triu(model.P,1);
            A1 = sparse(rand(model.n)<EA);
            A{i} = A1|A1';
            
            
        otherwise
            error('model name unavailable.')
    end
end
