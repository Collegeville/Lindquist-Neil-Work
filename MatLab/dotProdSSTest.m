classdef dotProdSSTest < matlab.unittest.TestCase
    
    properties
    end
    
    methods(TestMethodSetup)
    end
    
    methods(Test)
        function testLength4(testCase)
            testCase.verifyEqual(dotProdSS(single([5; 3.5; 1; 2]), single([7; -2; 8; -6])), 24);
        end
        
        function testLong(testCase)
            a = single([2; 3; 1.256; 1.3; 8.1; 100; -6; -8; 9; 10]);
            b = single([3.486; 3.6; 5.69; -3; 7; .75; 6; 2; 1.25; 2.25]);
            
            testCase.verifyEqual(dotProdSS(a, b), double(a')*double(b));
        end
        
        function testSingle(testCase)
            testCase.verifyEqual(dotProdSS(single(1.23568), single(6.547)), ...
                8.089996771007577);
        end
    end
end