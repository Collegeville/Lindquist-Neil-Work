function csrVScscTest

csrSum = 0;
cscSum = 0;

for ii = 1:5
    A = sprand(10000, 10000, 25/10000);
    x = rand(10000, 1);
    [csrR, csrC, csrV, cscR, cscC, cscV] = toFormats(A);

    csrSum = csrSum + timeit(@() csrMatvec(csrR, csrC, csrV, x));
    cscSum = cscSum + timeit(@() cscMatvec(cscR, cscC, cscV, x));
end

disp(['csr: ' num2str(csrSum/ii)])
disp(['csc: ' num2str(cscSum/ii)])

function [csrR, csrC, csrV, cscR, cscC, cscV] = toFormats(A)
    [m, n] = size(A);
    nz = nnz(A);
    [Ar, Ac, Av] = find(A);
    Ar = single(Ar);
    Ac = single(Ac);
    Av = single(Av);
    
    cscC = single(zeros(n, 1));
    c = 1;
    cscC(1) = 1;
    i = 1;
    while i <= nz
        if Ac(i) > c
            c = c+1;
            cscC(c) = i;
        else
            i = i+1;
        end
    end
    cscC(c+1:n+1) = nz+1;

    cscR = single(Ar);
    cscV = single(Av);

    rowOrdered = sortrows([Ar, Ac, Av]);
    
    csrR = single(zeros(m, 1));
    r = 1;
    csrR(1) = 1;
    i = 1;
    while i <= nz
        if rowOrdered(i, 1) > r
            r = r+1;
            csrR(r) = i;
        else
            i = i+1;
        end
    end
    csrR(r+1:m+1) = nz+1;

    csrC = rowOrdered(:, 2);
    csrV = rowOrdered(:, 3);
end

function b = csrMatvec(csrR, csrC, csrV, x)
    b = single(zeros(size(x)));

    nextEnd = csrR(2);
    row = 1;
    temp = 0;
    i = 1;
    while i <= nnz(A)
        if i == nextEnd
            b(row) = single(temp);
            temp = 0;
            row = row + 1;
            nextEnd = csrR(row+1);
        else
            temp = temp + double(csrV(i))*double(x(csrC(i)));
            i = i+1;
        end
    end
    b(row) = single(temp);
end

function b = cscMatvec(cscR, cscC, cscV, x)
    b = zeros(size(x));

    nextEnd = cscC(2);
    col = 1;
    xVal = double(x(col));
    i = 1;
    while i <= nnz(A)
        if i == nextEnd
            col = col + 1;
            xVal = double(x(col));
            nextEnd = cscC(col+1);
        else
            r = cscR(i);
            b(r) = b(r) + double(cscV(i))*xVal;
            i = i + 1;
        end
    end
    b = single(b);
end

end