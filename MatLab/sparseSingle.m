classdef sparseSingle
    %Contains a sparse matrix using CSR format
    properties
        rows    %the start index of each row
        cols    %the column indices that match the respective values
        vals    %the entries of the matrix
        m       %the number of rows in the matrix
        n       %the number of columns in the matrix
    end
    
    methods
        function self = sparseSingle(rows, cols, vals, m, n)
            if nargin == 1
                %matrix copy
                org = rows;
                [self.m, self.n] = size(org);
                
                self.rows = uint32(zeros(self.m+1, 1));
                nz = nnz(org);
                self.cols = uint32(zeros(nz, 1));
                self.vals = single(zeros(nz, 1));
                
                
                [ogRows, ogCols, ogVals] = find(org);
                
                i = 1;
                for r = 1:self.m
                    self.rows(r) = i;
                    for j = 1:nz
                        if ogRows(j) == r
                            self.cols(i) = ogCols(j);
                            self.vals(i) = ogVals(j);
                            i = i+1;
                        end
                    end
                end
                self.rows(self.m+1) = nz+1;
                
            else
                %creates an mxn sparse matrix with single precision
                if m <= 0 || n <= 0
                    error('Matrix dimentions must be positive')
                end
                if ~isequal(size(rows), [1, m])
                    if isequal(size(rows), [m, 1])
                        rows = rows';
                    else
                        error('Incorrect number of row indices')
                    end
                end
                if size(cols) ~= size(vals)
                    error('Number of row indices does not match number of values')
                end

                self.rows = uint32([rows, length(cols)+1]);
                self.cols = uint32(cols);
                self.vals = single(vals);
                self.m = m;
                self.n = n;
            end
        end
        
        function count = nnz(self)
            count = size(self.vals);
            count = count(1);
        end
        
        function varargout = size(self)
            if nargout <= 1
                varargout{1} = [self.m self.n];
            else
                varargout{1} = self.m;
                varargout{2} = self.n;
            end
        end
        
        function b = mtimes(A, x)
            
            %Computes the matrix-vector product.
            %The matrix should be a single precision sparse matrix
            %The vector should be a single precision dense vector

            if isscalar(x)
                b = A.vals*x;
            elseif isscalar(A);
                b = x.vals*A;
            else
            
                [aRows, aCols] = size(A);
                [xRows, xCols] = size(x);

                if aCols ~= xRows
                    error('Inner matrix dimensions must agree.');
                end
                if xCols ~= 1
                    error('Second operand must be a column vector');
                end

                r = A.rows;
                c = A.cols;
                v = A.vals;

                b = single(zeros(aRows, xCols));

                nextEnd = r(2);
                row = 1;
                temp = 0;
                i = 1;
                while i <= nnz(A)
                    if i == nextEnd
                        b(row) = single(temp);
                        temp = 0;
                        row = row + 1;
                        nextEnd = r(row+1);
                    end
                    temp = temp + double(v(i))*double(x(c(i)));
                    i = i+1;
                end
                b(row) = single(temp);
            end
        end
        
        function C = plus(A, B)
            if isscalar(B)
                if B ~= 0
                    C = full(A)+B;
                else
                    C = sparseSingle(A);
                end
            elseif isscalar(A)
                if A ~= 0
                    C = full(B)+A;
                else
                    C = sparseSingle(B);
                end
            else
                if issparse(A) && issparse(B)
                    
                    if ~isequal(size(A), size(B))
                        error('Matrices must be the same size for addition')
                    end
                    
                    if isa(A, 'sparseSingle')
                        self = A;
                        other = B;
                    else
                        self = B;
                        other = A;
                    end
                    
                    [r1, c1, v1] = find(self);
                    [r2, c2, v2] = find(other);
                    
                    entries = [r1, c1, v1; r2, c2, v2];
                    
                    entries = sortrows(entries, [1, 2]);
                    duplicates = 0;
                    len = size(entries, 1);
                    for i = 2:len
                        if entries(i-1, 1) == entries(i, 1) && entries(i-1, 2) == entries(i, 2)
                            entries(i-1, 3) = entries(i-1, 3)+entries(i, 3);
                            entries(i, 1) = -1;
                            entries(i, 2) = -1;
                            duplicates = duplicates + 1;
                        end
                    end
                    entries = sortrows(entries, [1, 2]);
                    entries = entries(duplicates+1:len, :);
                    r = entries(:, 1);
                    c = entries(:, 2);
                    v = entries(:, 3);
                    row = zeros(self.m, 1);
                    cRow = 1;
                    i = 1;
                    while i <= size(entries, 1);
                        if r(i) > cRow
                            cRow = cRow+1;
                            row(cRow) = i;
                        else
                            i = i + 1;
                        end
                    end
                    
                    C = sparseSingle(row, c, v, self.m, self.n);
                else
                    C = full(A)+full(B);
                end
            end
        end
        
        function varargout = subsref(obj, s)
            switch s(1).type
                case '()'
                    subs = s(1).subs;
                    if length(subs) == 1
                        %linear indexed, just get as list
                        [r, c] = ind2sub(size(obj), subs{1});
                        v = zeros(length(r), 1);
                        for i = 1:length(r)
                            
                            v(i) = subsref(obj, struct('type', {'()'}, 'subs', {{r(i), c(i)}, 0}));
                        end
                        varargout{1} = v;
                    else
                        r = subs{1};
                        c = subs{2};

                        if r == ':'
                            r = 1:obj.m;
                        end
                        if c == ':'
                            c = 1:obj.n;
                        end

                        if length(r) * length(c) == 1
                            i = obj.rows(r);
                            varargout{1} = 0;
                            while i < obj.rows(r+1)
                                if obj.cols(i) == c
                                    varargout{1} = obj.vals(i);
                                    break
                                end
                                i = i+1;
                            end
                        else
                            i = obj.rows(r(1));
                            j = 1;
                            newRows = zeros(1, length(r));
                            nnz = obj.rows(r(end)+1)-obj.rows(r(1));
                            newCols = zeros(1, nnz);
                            newVals = zeros(1, nnz);
                            currentRow = r(1)-1;
                            while i < obj.rows(r(end)+1)
                                 if i >= obj.rows(currentRow+1)
                                     currentRow = currentRow + 1;
                                     newRows(currentRow-r(1)+1) = j;
                                 else
                                     newCols(j) = obj.cols(i);
                                     newVals(j) = obj.vals(i);
                                     i = i+1;
                                     j = j+1;
                                 end
                            end
                            varargout{1} = sparseSingle(newRows, newCols, newVals, length(r), length(c));
                        end
                    end
                    
                case '.'
                    
                    varargout{1} = builtin('subsref', obj, s);
                case '{}'
                    varargout{1} = builtin('subsref', obj, s);
                otherwise
                    error('Not a valid indexing expression')
            end
        end
        
        function s = sparse(self)
            %converts the single precision sparse matrix to the builtin,
            %double precision sparse matrix
            r = zeros(size(self.cols));
            currentRow = 1;
            i = 1;
            while i <= length(self.vals)
                if i >= self.rows(currentRow+1)
                    currentRow = currentRow + 1;
                else
                    if currentRow > 479
                        disp(i);
                    end
                    r(i) = currentRow;
                    i = i + 1;
                end
            end
            s = sparse(r, double(self.cols), double(self.vals), self.m, self.n);
        end
        
        function s = full(self)
            s = single(zeros(self.m, self.n));
            
            r = 1;
            i = 1;
            while i<=nnz(self)
                if i >= self.rows(r+1)
                    r = r + 1;
                else
                    s(r, self.cols(i)) = self.vals(i);
                    i = i + 1;
                end
            end
        end
        
        function [row, col, v] = find(self, n, direction)
            
            
            nz = nnz(self);
            
            if nargin < 3
                direction = 'first';
                
                if nargin < 2
                    n = nz;
                end
            end
            
            if n > nz
                error('n must be less than the number of non-zero elements in the matrix');
            end
            
            
            row = zeros(nz, 1);
            col = zeros(nz, 1);
            v = zeros(nz, 1);
                
            r = 1;
            i = 1;
            while i <= nz
                if i >= self.rows(r+1)
                    r = r + 1;
                else
                    row(i) = r;
                    col(i) = self.cols(i);
                    v(i)   = self.vals(i);
                    i = i+1;
                end
            end
            
            sorted = sortrows([col, row, v]);
            if strcmp('first', direction)
                col = sorted(1:n, 1);
                row = sorted(1:n, 2);
                v   = sorted(1:n, 3);
            else
                col = sorted(nz-n+1:nz, 1);
                row = sorted(nz-n+1:nz, 2);
                v   = sorted(nz-n+1:nz, 3);
            end

            if nargout <= 1
                row = sub2ind(size(self), row, col);
            end
            
            if isrow(self)
                row = row';
                col = col';
                v   = v';
            end
        end
    
        function disp(self)
            if self.m ~= 0 && self.n ~= 0
                [r, c, v] = find(self);
                for i = 1:length(r)
                    fprintf('   (%g,%g)     %g\n', r(i), c(i), v(i));
                end
            end
            %if has a dimension 0, print nothing
        end
    end
    
    methods(Static)
        function TF = issparse()
            TF = true;
        end
    end
end

    
    