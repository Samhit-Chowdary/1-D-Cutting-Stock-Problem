logLength = 100;
lengthlist = [14; 31; 36; 45];
quantity = [211; 395; 610; 97];
nLengths = length(lengthlist);
patterns = diag(floor(logLength./lengthlist));
nPatterns = size(patterns,2);
lb2 = zeros(nLengths,1);
A2 = lengthlist';
b2 = logLength;
lpopts = optimoptions('linprog','Display','off');
ipopts = optimoptions('intlinprog',lpopts);
reducedCost = -Inf;
reducedCostTolerance = -0.0001;
exitflag = 1;
itr=0;
while reducedCost < reducedCostTolerance && exitflag > 0
    itr=itr+1;
    lb = zeros(nPatterns,1);
    f = lb + 1;
    A = -patterns;
    b = -quantity;
    
    [values,nLogs,exitflag,~,lambda] = linprog(f,A,b,[],[],lb,[],lpopts); 
    fprintf('iteration %g\nvalues of x\n',itr);
    disp(values);
    fprintf('values of dual variables\n')
    disp(lambda.ineqlin);
    if exitflag > 0
        fprintf('Using %g logs\n',nLogs);
        % Now generate a new pattern, if possible
        f2 = -lambda.ineqlin;
        [values,reducedCost,pexitflag] = intlinprog(f2,1:nLengths,A2,b2,[],[],lb2,[],ipopts);
        reducedCost = 1 + reducedCost; % continue if this reducedCost is negative
        newpattern = round(values);
        if pexitflag > 0 && reducedCost < reducedCostTolerance
            patterns = [patterns newpattern];
            fprintf('new pattern\n');
            disp(newpattern);
            nPatterns = nPatterns + 1;
        end
    end
end
if exitflag <= 0 
    disp('Error in column generation phase')
else
    [values,logsUsed,exitflag] = intlinprog(f,1:length(lb),A,b,[],[],lb,[],[],ipopts);
    if exitflag > 0
        values = round(values);
        logsUsed = round(logsUsed);
        fprintf('Optimal solution uses %g logs\n', logsUsed);
        totalwaste = sum((patterns*values - quantity).*lengthlist); % waste due to overproduction
        for j = 1:size(values)
            if values(j) > 0
                fprintf('Cut %g logs with pattern\n',values(j));
                for w = 1:size(patterns,1)
                    if patterns(w,j) > 0
                        fprintf('    %d cut(s) of length %d\n', patterns(w,j),lengthlist(w));
                    end
                end
                wastej = logLength - dot(patterns(:,j),lengthlist); % waste due to pattern inefficiency
                totalwaste = totalwaste + wastej;
            fprintf('    Waste of this pattern is %g\n', wastej);
            end
        end
        fprintf('Total waste in this problem is %g.\n',totalwaste);
    else 
        disp('Error in final optimization')
    end
end