classdef normSTest < matlab.unittest.TestCase
    
    methods(Test)
        function testLength3_4_5(testCase)
            testCase.verifyEqual(normS(single([3; 4])), 5);
        end
        
        function testLong(testCase)
            a = single([2; 3; 1.256; 1.3; 8.1; 100; -6; -8; 9; 10]);
            
            testCase.verifyEqual(normS(a), norm(double(a)));
        end
        
        function testSingle(testCase)
            testCase.verifyEqual(normS(single(1.23568)), double(single(1.23568)));
        end
    end
end