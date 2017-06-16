classdef doubleMatvecTest < matlab.unittest.TestCase
    
    properties
        m1      % a test matrix
        m2      % a test matrix
        one     % a 1x1 matrix
        empty   % a 6x8 zero filled matrix
    end
    
    methods(TestMethodSetup)
        function createTestObjects(testCase)
            testCase.m1 = sparse([1, 1, 0; 0, 0, 2; 1, 0, 0]);
            testCase.m2 = sparse([.75, 1, 0, 0, 0; 0, .5, 1, 0, 0; 0, 0, 1, 2, 3; 0, 0, 0, 1, 0]);
            testCase.one = sparse(5);
            testCase.empty = sparse(zeros(6, 8));
        end
    end
    
    methods(Test)
        
        function testMTimesOnes(testCase)
            b = doubleMatvec(testCase.m1', ones(3, 1));
            testCase.verifyEqual(b, testCase.m1*ones(3, 1));

            b = doubleMatvec(testCase.m2', ones(5, 1));
            testCase.verifyEqual(b, testCase.m2*ones(5, 1));
            
            b = doubleMatvec(testCase.one', 1);
            testCase.verifyEqual(b, 5);

            b = doubleMatvec(testCase.empty', ones(8, 1));
            testCase.verifyEqual(b, testCase.empty*ones(8, 1));
        end
        
        function testMTimes(testCase)
            testCase.verifyEqual(doubleMatvec(testCase.m1', [4; 9; 6]), ...
                    testCase.m1*[4; 9; 6]);
                    
            testCase.verifyEqual(doubleMatvec(testCase.m2', [1.02; -3; 4.05; -.95; 10]), ...
                    testCase.m2*[1.02; -3; 4.05; -.95; 10])
                
            testCase.verifyEqual(doubleMatvec(testCase.m2', [-3; -4; -.002; -12.103; -3.4]), ...
                    testCase.m2*[-3; -4; -.002; -12.103; -3.4]);
        end
    end
end
