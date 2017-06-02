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
                disp('copying')
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
                if size(rows) ~= [1, m]
                    error('Incorrect number of row indices')
                end
                if size(cols) ~= size(vals)
                    error('Number of row indices does not match number of values')
                end

                self.rows = [rows m+1];
                self.cols = cols;
                self.vals = single(vals);
                self.m = m;
                self.n = n;
            end
        end
        
        function count = nnz(self)
            count = size(self.vals);
            count = count(1);
        end
        
        function [m n] = size(self)
            if nargout == 1
                m = [self.m self.n];
            else
                m = self.m;
                n = self.n;
            end
        end
        
        function b = mtimes(A, x)
            
            %Computes the matrix-vector product.
            %The matrix should be a single precision sparse matrix
            %The vector should be a single precision dense vector


            r = A.rows;
            c = A.cols;
            v = A.vals;

            b = single(zeros(size(x)));

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
end