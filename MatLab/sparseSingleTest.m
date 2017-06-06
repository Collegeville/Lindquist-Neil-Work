classdef sparseSingleTest < matlab.unittest.TestCase
    
    properties
        f1      % the full representation of ssp1
        ssp1    % a 3x3 sparseSingle matrix
        f2      % the full representation of ssp2
        ssp2    % a 4x5 sparseSingle matrix
        one     % a 1x1 sparseSingle matrix
        empty   % a 6x8 zero filled sparseSingle matrix
    end
    
    methods(TestMethodSetup)
        function createTestObjects(testCase)
            testCase.f1 = [1, 1, 0; 0, 0, 2; 1, 0, 0];
            testCase.ssp1 = sparseSingle(testCase.f1);
            testCase.f2 = [.75, 1, 0, 0, 0; 0, .5, 1, 0, 0; 0, 0, 1, 2, 3; 0, 0, 0, 1, 0];
            testCase.ssp2 = sparseSingle(single(testCase.f2));
            testCase.one = sparseSingle(5);
            testCase.empty = sparseSingle(zeros(6, 8));
        end
    end
    
    methods(Test)
        
        function testCopyConstuctorDenseDouble(testCase)
            testCase.verifyEqual(testCase.ssp1.m, 3);
            testCase.verifyEqual(testCase.ssp1.n, 3);
            
            testCase.verifyEqual(testCase.ssp1.rows, uint32([1; 3; 4; 5]));
            testCase.verifyEqual(testCase.ssp1.cols, uint32([1; 2; 3; 1]));
            testCase.verifyEqual(testCase.ssp1.vals, single([1; 1; 2; 1]));
        end
        
        function testCopyConstuctorDenseSingle(testCase)
            testCase.verifyEqual(testCase.ssp2.m, 4);
            testCase.verifyEqual(testCase.ssp2.n, 5);
            
            testCase.verifyEqual(testCase.ssp2.rows, uint32([1; 3; 5; 8; 9]));
            testCase.verifyEqual(testCase.ssp2.cols, uint32([1; 2; 2; 3; 3; 4; 5; 4]));
            testCase.verifyEqual(testCase.ssp2.vals, single([.75; 1; .5; 1; 1; 2; 3; 1]));
        end
        
        function testCopyConstructorOneElement(testCase)
            testCase.verifyEqual(testCase.one.m, 1);
            testCase.verifyEqual(testCase.one.n, 1);
            
            testCase.verifyEqual(testCase.one.rows, uint32([1; 2]));
            testCase.verifyEqual(testCase.one.cols, uint32(1));
            testCase.verifyEqual(testCase.one.vals, single(5));
        end
        
        function testCopyConstructorEmpty(testCase)
            testCase.verifyEqual(testCase.empty.m, 6);
            testCase.verifyEqual(testCase.empty.n, 8);
            
            testCase.verifyEqual(testCase.empty.rows, uint32([1;1;1;1;1;1;1;]));
            testCase.verifyEqual(testCase.empty.cols, uint32(zeros(0, 1)));
            testCase.verifyEqual(testCase.empty.vals, single(zeros(0, 1)));
        end
        
        function testCopyConstructorSparseDouble(testCase)
            ssp = sparseSingle(sparse([3, 0; 0, 1; 0, 0; 0, 2; 1, 1]));
            
            testCase.verifyEqual(ssp.m, 5);
            testCase.verifyEqual(ssp.n, 2);
            
            testCase.verifyEqual(ssp.rows, uint32([1; 2; 3; 3; 4; 6]));
            testCase.verifyEqual(ssp.cols, uint32([1; 2; 2; 1; 2]));
            testCase.verifyEqual(ssp.vals, single([3; 1; 2; 1; 1]));
        end
        
        function testCopyConstructorSparseSingle(testCase)
            ssp = sparseSingle(testCase.ssp1);
            
            testCase.verifyEqual(ssp.m, 3);
            testCase.verifyEqual(ssp.n, 3);
            
            testCase.verifyEqual(ssp.rows, uint32([1; 3; 4; 5]));
            testCase.verifyEqual(ssp.cols, uint32([1; 2; 3; 1]));
            testCase.verifyEqual(ssp.vals, single([1; 1; 2; 1]));
        end
        
        function testCSRConstructorDouble(testCase)
            r = [1; 3; 6; 10; 12;];
            c = [1; 4; 1; 2; 4; 1; 3; 4; 5; 3; 4; 5;];
            v = [1; 2; 3; 4; 5; 6; 7; 8; 9; 10; 11; 12;];
            
            ssp = sparseSingle(r, c, v, 5, 5);
            
            testCase.verifyEqual(ssp.m, 5);
            testCase.verifyEqual(ssp.n, 5);
            
            testCase.verifyEqual(ssp.rows, uint32([r; 13]));
            testCase.verifyEqual(ssp.cols, uint32(c));
            testCase.verifyEqual(ssp.vals, single(v));
        end
        
        function testCSRConstructorSingle(testCase)
            r = single([1; 3; 6; 9]);
            c = single([1; 3; 1; 4; 5; 2; 3; 6; 2; 5]);
            v = single([1; 2; .25; 1; 1; 1; 1; 3; .3; 1]);
            
            ssp = sparseSingle(r, c, v, 4, 6);
            
            testCase.verifyEqual(ssp.m, 4);
            testCase.verifyEqual(ssp.n, 6);
            
            testCase.verifyEqual(ssp.rows, uint32([r; 11]));
            testCase.verifyEqual(ssp.cols, uint32(c));
            testCase.verifyEqual(ssp.vals, v);
        end
    
        function testNNZ(testCase)
            testCase.verifyEqual(nnz(testCase.ssp1), 4);

            testCase.verifyEqual(nnz(testCase.ssp2), 8);

            testCase.verifyEqual(nnz(testCase.one), 1);

            testCase.verifyEqual(nnz(testCase.empty), 0);
        end
        
        function testSizeSingleOut(testCase)
            sz = size(testCase.ssp1);
            testCase.verifyEqual(sz, [3 3]);

            sz = size(testCase.ssp2);
            testCase.verifyEqual(sz, [4 5]);

            sz = size(testCase.one);
            testCase.verifyEqual(sz, [1 1]);
            
            sz = size(testCase.empty);
            testCase.verifyEqual(sz, [6 8]);
        end
        
        function testSizeDoubleOut(testCase)
            [m, n] = size(testCase.ssp1);
            testCase.verifyEqual(m, 3);
            testCase.verifyEqual(n, 3);
            
            [m, n] = size(testCase.ssp2);
            testCase.verifyEqual(m, 4);
            testCase.verifyEqual(n, 5);
            
            [m, n] = size(testCase.one);
            testCase.verifyEqual(m, 1);
            testCase.verifyEqual(n, 1);
            
            [m, n] = size(testCase.empty);
            testCase.verifyEqual(m, 6);
            testCase.verifyEqual(n, 8);
        end
        
        function testSizeDim1(testCase)
            szdim = size(testCase.ssp1, 1);
            testCase.verifyEqual(szdim, 3);
            
            szdim = size(testCase.ssp2, 1);
            testCase.verifyEqual(szdim, 4);
            
            szdim = size(testCase.one, 1);
            testCase.verifyEqual(szdim, 1);
            
            szdim = size(testCase.empty, 1);
            testCase.verifyEqual(szdim, 6);
        end
        
        function testSizeDim2(testCase)
            m = size(testCase.ssp1, 2);
            testCase.verifyEqual(m, 3);

            szdim = size(testCase.ssp2, 2);
            testCase.verifyEqual(szdim, 5);
            
            szdim = size(testCase.one, 2);
            testCase.verifyEqual(szdim, 1);

            szdim = size(testCase.empty, 2);
            testCase.verifyEqual(szdim, 8);
        end
        
        function testMTimesOnes(testCase)
            b = testCase.ssp1*ones(3, 1);
            testCase.verifyEqual(b, single([2; 2; 1]));

            b = testCase.ssp2*ones(5, 1);
            testCase.verifyEqual(b, single([1.75; 1.5; 6; 1]));
            
            b = testCase.one*1;
            testCase.verifyEqual(b, sparseSingle(5));

            b = testCase.empty*ones(8, 1);
            testCase.verifyEqual(b, single(zeros(6, 1)));
        end
        
        function testMTimesScalard(testCase)
            C = testCase.ssp1*5;
            testCase.verifyEqual(C, sparseSingle(testCase.f1*5));

            C = testCase.ssp1*12;
            testCase.verifyEqual(C, sparseSingle(testCase.f1*12));

            C = testCase.one*9;
            testCase.verifyEqual(C, sparseSingle(45));

            C = testCase.empty*12.5;
            testCase.verifyEqual(C, sparseSingle(zeros(6, 8)));

            C = 8*testCase.ssp1;
            testCase.verifyEqual(C, sparseSingle(testCase.f1*8));
        end
        
        function testPlusdoubled(testCase)
            testCase.verifyEqual(testCase.ssp1+testCase.ssp1, testCase.ssp1*2);

            testCase.verifyEqual(testCase.ssp2+testCase.ssp2, testCase.ssp2*2);
        end
        
        function testPlus2Ones(testCase)
            testCase.verifyEqual(testCase.ssp2+ones(4, 5), single(testCase.f2+1));
        end
        
        function testPlus2LeftOnes(testCase)
            testCase.verifyEqual(ones(4, 5)+testCase.ssp2, single(testCase.f2+1));
        end
        
        function testPlusScalar(testCase)
            testCase.verifyEqual(testCase.ssp2+1, single(testCase.f2+1));

            testCase.verifyEqual(1+testCase.ssp2, single(testCase.f2+1));
        end
        
        function testSubsrefDot(testCase)
            testCase.verifyEqual(testCase.ssp1.rows, uint32([1; 3; 4; 5]));
        end
        
        function testSubsrefSingleElementSubs(testCase)
            testCase.verifyEqual(testCase.ssp1(1, 1), single(1));
            testCase.verifyEqual(testCase.ssp1(2, 1), single(0));
            testCase.verifyEqual(testCase.ssp1(3, 3), single(0));
        end
        
        function testSubsrefSingleElementLinearIndex(testCase)
            testCase.verifyEqual(testCase.ssp1(1), single(1));
            testCase.verifyEqual(testCase.ssp1(2), single(0));
            testCase.verifyEqual(testCase.ssp1(9), single(0));
        end
        
        function testSubsrefRangeSubs(testCase)
            testCase.verifyEqual(testCase.ssp1(1:2, 1:2), sparseSingle(testCase.f1(1:2, 1:2)));
            testCase.verifyEqual(testCase.ssp2(3:4, :), sparseSingle(testCase.f2(3:4, :)));
            testCase.verifyEqual(testCase.ssp2(:, 2:4), sparseSingle(testCase.f2(:, 2:4)));
        end
        
        function testSubsrefVectLinearIndexes(testCase)
            testCase.verifyEqual(testCase.ssp1([2, 8, 1]), single([0, 2, 1]));
        end
        
        function testFindAllElements(testCase)
            act = find(testCase.ssp2);
            exp = find(testCase.f2);
            testCase.verifyEqual(act, exp);
            
            [actR, actC] = find(testCase.ssp2);
            [expR, expC] = find(testCase.f2);
            testCase.verifyEqual(actR, expR);
            testCase.verifyEqual(actC, expC);
            
            [actR, actC, actV] = find(testCase.ssp2);
            [expR, expC, expV] = find(testCase.f2);
            testCase.verifyEqual(actR, expR);
            testCase.verifyEqual(actC, expC);
            testCase.verifyEqual(actV, expV);
        end
        
        function testFindNElements(testCase)
            n = nnz(testCase.ssp2);
            act = find(testCase.ssp2, n);
            exp = find(testCase.f2, n);
            testCase.verifyEqual(act, exp);
            
            [actR, actC] = find(testCase.ssp2, n);
            [expR, expC] = find(testCase.f2, n);
            testCase.verifyEqual(actR, expR);
            testCase.verifyEqual(actC, expC);
            
            [actR, actC, actV] = find(testCase.ssp2, n);
            [expR, expC, expV] = find(testCase.f2, n);
            testCase.verifyEqual(actR, expR);
            testCase.verifyEqual(actC, expC);
            testCase.verifyEqual(actV, expV);
        end
        
        function testFind4Elements(testCase)
            act = find(testCase.ssp2, 4);
            exp = find(testCase.f2, 4);
            testCase.verifyEqual(act, exp);
            
            [actR, actC] = find(testCase.ssp2, 4);
            [expR, expC] = find(testCase.f2, 4);
            testCase.verifyEqual(actR, expR);
            testCase.verifyEqual(actC, expC);
            
            [actR, actC, actV] = find(testCase.ssp2, 4);
            [expR, expC, expV] = find(testCase.f2, 4);
            testCase.verifyEqual(actR, expR);
            testCase.verifyEqual(actC, expC);
            testCase.verifyEqual(actV, expV);
        end
        
        function testFindFirstNElements(testCase)
            n = nnz(testCase.ssp2);
            act = find(testCase.ssp2, n, 'first');
            exp = find(testCase.f2, n);
            testCase.verifyEqual(act, exp);
            
            [actR, actC] = find(testCase.ssp2, n, 'first');
            [expR, expC] = find(testCase.f2, n);
            testCase.verifyEqual(actR, expR);
            testCase.verifyEqual(actC, expC);
            
            [actR, actC, actV] = find(testCase.ssp2, n, 'first');
            [expR, expC, expV] = find(testCase.f2, n);
            testCase.verifyEqual(actR, expR);
            testCase.verifyEqual(actC, expC);
            testCase.verifyEqual(actV, expV);
        end
        
        function testFindFirst4Elements(testCase)
            act = find(testCase.ssp2, 4, 'first');
            exp = find(testCase.f2, 4);
            testCase.verifyEqual(act, exp);
            
            [actR, actC] = find(testCase.ssp2, 4, 'first');
            [expR, expC] = find(testCase.f2, 4);
            testCase.verifyEqual(actR, expR);
            testCase.verifyEqual(actC, expC);
            
            [actR, actC, actV] = find(testCase.ssp2, 4, 'first');
            [expR, expC, expV] = find(testCase.f2, 4);
            testCase.verifyEqual(actR, expR);
            testCase.verifyEqual(actC, expC);
            testCase.verifyEqual(actV, expV);
        end
        
        function testFindLastNElements(testCase)
            n = nnz(testCase.ssp2);
            act = find(testCase.ssp2, n, 'last');
            exp = find(testCase.f2, n, 'last');
            testCase.verifyEqual(act, exp);
            
            [actR, actC] = find(testCase.ssp2, n, 'last');
            [expR, expC] = find(testCase.f2, n, 'last');
            testCase.verifyEqual(actR, expR);
            testCase.verifyEqual(actC, expC);
            
            [actR, actC, actV] = find(testCase.ssp2, n, 'last');
            [expR, expC, expV] = find(testCase.f2, n, 'last');
            testCase.verifyEqual(actR, expR);
            testCase.verifyEqual(actC, expC);
            testCase.verifyEqual(actV, expV);
        end
        
        function testFindLast4Elements(testCase)
            act = find(testCase.ssp2, 4, 'last');
            exp = find(testCase.f2, 4, 'last');
            testCase.verifyEqual(act, exp);
            
            [actR, actC] = find(testCase.ssp2, 4, 'last');
            [expR, expC] = find(testCase.f2, 4, 'last');
            testCase.verifyEqual(actR, expR);
            testCase.verifyEqual(actC, expC);
            
            [actR, actC, actV] = find(testCase.ssp2, 4, 'last');
            [expR, expC, expV] = find(testCase.f2, 4, 'last');
            testCase.verifyEqual(actR, expR);
            testCase.verifyEqual(actC, expC);
            testCase.verifyEqual(actV, expV);
        end
        
        function testSparse2(testCase)
            dsp = sparse(testCase.ssp2);
            
            testCase.verifyTrue(isa(dsp, 'double'));
            testCase.verifyTrue(issparse(dsp));
            
            [actR, actC, actV] = find(dsp);
            [expR, expC, expV] = find(testCase.ssp2);
            
            testCase.verifyEqual(actR, expR);
            testCase.verifyEqual(actC, expC);
            testCase.verifyEqual(actV, expV);
        end
        
        function testFull2(testCase)
            sf = full(testCase.ssp2);
            
            testCase.verifyTrue(isa(sf, 'single'));
            testCase.verifyFalse(issparse(sf));
            
            testCase.verifyEqual(sf, single(testCase.f2));
        end
        
        function testIssparse(testCase)
            testCase.verifyTrue(issparse(testCase.ssp1));
            testCase.verifyTrue(issparse(testCase.ssp2));
            testCase.verifyTrue(issparse(testCase.one));
            testCase.verifyTrue(issparse(testCase.empty));
        end
        
    end
end