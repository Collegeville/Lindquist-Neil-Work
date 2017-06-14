classdef residualTest < matlab.unittest.TestCase
    
    properties
    end
    
    methods(TestMethodSetup)
    end
    
    methods(Test)
        function testNoResidual(testCase)
            A = sparseSingle(eye(10));
            x = single([1; 5; 6.002; 7; 3; 123; 98; 3.222; 8; 1]);
            b = single(x);

            testCase.verifyEqual(residual(b, A.rows, A.cols, A.vals, x), ...
                zeros(10, 1));
        end
        
        function testResidual(testCase)
            sA = sparseSingle([.75, 1, 0, 0, 0; 0, .5, 1, 0, 0; 0, 0, 1, 2, 3; 0, 0, 0, 1, 0; .11, 0, 0, 0, .25]);
            sx = single([.654; 12; 98; 12; 9]);
            sb = single([4; 5; 69; .25; 4]);
            A = sparse(sA);
            x = double(sx);
            b = double(sb);
            testCase.verifyEqual(residual(sb, sA.rows, sA.cols, sA.vals, sx), ...
                b-A*x);
        end
    end
end