function tests = matlabTest
    tests = functiontests(localfunctions);
end

function testMakeCorrMatrixZeros(testCase)
    a = [0 0 1 -1 0 0; 1 -1 -1 0 1 1; 0 -1 1 -1 1 1];
    assertEqual(testCase,makeCorrMatrix(a),[0 0 0;0 1 .5774;0 .5774 1],'AbsTol',1e-2);
end

function testMakeCorrMatrixSame(testCase)
    a = [0 0 1 -1 0 0; 1 1 1 0 1 1; 0 -1 1 -1 1 1];
    assertEqual(testCase,makeCorrMatrix(a),[0 0 0;0 0 0;0 0 1],'AbsTol',1e-2);
end

function testPriorFileFormat(testCase)
    prior = dlmread('../ChubukovTest/prior.txt');
    assertTrue(testCase, max(prior(:,1))<=44 && max(prior(:,2))<=44);
    assertTrue(testCase, size(prior,2)==3);
    assertTrue(testCase, all(all(mod(prior,1)==0)));
end

function testPertFileValues(testCase)
    pert = dlmread('../ChubukovTest/pert.txt');
    first8MaxVals = max(pert(:,1:8),2);
    first8MaxVals = first8MaxVals/max(first8MaxVals);
    assertTrue(testCase, all(first8MaxVals>=1/3));
    assertTrue(testCase, all(all(pert>=0)));
end