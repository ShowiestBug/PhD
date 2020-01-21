function [x,fval,StopReason,IterationCounter,Time] = PgaV4r22(fun,nvars,lb,ub,Workers,Range,PopulationSize,ToVal,CountVal)
    if(matlabpool('size')<2)  
        msg = 'Error occurred. Matlabpool not open';
        error(msg)
    end
    tstart = tic;
    if nargin < 9, CountVal = 100;
        if nargin < 8, ToVal = -inf;
            if nargin < 7, PopulationSize = 20;
                if nargin <6, Range = [];
                    if nargin <5, Workers = matlabpool('size');
                        if nargin < 4, lb = [];
                            if nargin < 3,ub  = [];
                            end
                        end
                    end
                end
            end 
        end
    end
    mode = 2;   
    EliteKids = 3;
    MutaKids = 6;  
    CrossParents = (2*(PopulationSize - EliteKids - MutaKids))+ MutaKids;
    ExitFlag = [];
    IterationCounter = 0;
    GenBest = [];
    [Population,Score,Pass] = PgaInI(PopulationSize,nvars,Workers,lb,ub,Range);
    while isempty(ExitFlag)
        IterationCounter = IterationCounter + 1;   
        [Population,Score] = PgaSolver(fun,Population,Score,EliteKids,MutaKids,CrossParents,IterationCounter,Workers,mode,ub,lb,CountVal);
        [TempVal1,FirstIndex] = min(Score);
        [TempVal2,SecondIndex] = min(TempVal1);  
        GenBest(IterationCounter) = TempVal2;
        if(IterationCounter < 2)
            BestScore = Score(FirstIndex(SecondIndex),SecondIndex);
            BestGenes = Population(FirstIndex(SecondIndex),:,SecondIndex);
        elseif(IterationCounter > 1)
            if(Score(FirstIndex(SecondIndex),SecondIndex) < BestScore)
                BestScore = Score(FirstIndex(SecondIndex),SecondIndex);
                BestGenes = Population(FirstIndex(SecondIndex),:,SecondIndex);
            end
        end
        [ExitFlag,StopReason] = PgaStop(IterationCounter,GenBest,CountVal,ToVal,1e-6,Workers);
        if(mode <3)
            [Population,Score] = PgaPass(Population,Score,Pass,Workers);
        end  
    end
    Time = toc(tstart);
    x = BestGenes;
    fval = BestScore;
end

function [Population,Score,Pass] = PgaInI(PopulationSize,nvars,Workers,lb,ub,Range)
    if(isempty(Range) && isempty(lb) && isempty(ub))
        lr = zeros(1,nvars);
        ur = ones(1,nvars);
    elseif(~isempty(Range) && isempty(lb) && isempty(ub))
        lr = zeros(1,nvars)+ Range(1);
        ur = zeros(1,nvars)+ Range(2);
    elseif(isempty(Range) && ~isempty(lb) && ~isempty(ub))
        lr = lb;
        ur = ub;
    elseif(~isempty(Range) && ~isempty(lb) && ~isempty(ub))
        lr = zeros(1,nvars)+ Range(1);
        ur = zeros(1,nvars)+ Range(2);
    end 
    Population = zeros(PopulationSize,nvars,Workers);
    for i=1:nvars
        Population(:,i,:) = lr(i) + (ur(i) - lr(i))* rand(PopulationSize,1,Workers);        
    end
    Score = zeros(PopulationSize,Workers);
    Pass = Workers - 1;
    if(Pass > 0.2 * int8(PopulationSize))
        Pass = 0.2 * int8(PopulationSize);
    end 
end

function [Population,Score] = PgaSolver(fun,Population,Score,EliteKids,MutaKids,CrossParents,IterationCounter,Workers,mode,ub,lb,CountVal)

    if(mode == 1)   
        spmd
            % Rozdzielenie populacji i przystosowania
            SubPopulation = Population(:,:,labindex);

            % Przystosowanie populacji startowej            
            if(IterationCounter < 2)
                SubScore = PgaFval(fun,SubPopulation,1);                
            elseif(IterationCounter > 1)
                SubScore = Score(:,labindex);
            end
            
            % Wyznaczenie Elity z danej populacji
            [EliteScore,EliteIndex] = sort(SubScore); 
            Elite = SubPopulation(EliteIndex(1:EliteKids),:);
            EliteScore = EliteScore(1:EliteKids,1);
            
            % Przeskalowanie przystosowania
            Expectation = PgaRank(SubScore,EliteIndex,CrossParents);
            
            % Selekcja osobnikow do krzyzowania i mutacji
            Selection = PgaSelect(Expectation,CrossParents);
            
            % krzyzowanie
            CrossChild = PgaCross(Selection,CrossParents,MutaKids,SubPopulation);
            
            % mutacja do sprawdzenia
            MutaChild = PgaMuta(Selection,SubPopulation,MutaKids,(IterationCounter-1)*Workers,CountVal); 
            
            % Odtworzenie populacji i przystosowania %% bounds
            SubPopulation = [Elite; CrossChild; MutaChild];
            SubPopulation = PgaBounds(SubPopulation,ub,lb);
            SubScore = PgaFval(fun,SubPopulation,2,EliteScore);  
        end
        % zkladanie populacji i scora
        Population = SubPopulation{1};
        for i = 2 : Workers
            Population = cat(3,Population,SubPopulation{i});
        end
        Score = SubScore{1};
        for i = 2 : Workers 
            Score = cat(2,Score,SubScore{i});
        end
    elseif(mode == 2)
        parfor W = 1 : Workers
            
            %Rozdzielenie populacji i przystosowania
            SubPopulation = Population(:,:,W);            
            if(IterationCounter < 2)
                SubScore = PgaFval(fun,SubPopulation,1);                
            elseif(IterationCounter > 1)
                SubScore = Score(:,W);
            end
            
            %Wyznaczenie Elity
            [EliteScore,EliteIndex] = sort(SubScore);
            Elite = SubPopulation(EliteIndex(1:EliteKids),:);
            EliteScore = EliteScore(1:EliteKids,1);
            
            %Przeskalowanie przystosowania
            Expectation = PgaRank(SubScore,EliteIndex,CrossParents);
            
            %Selekcja do mutacji i krzyzowania
            Selection = PgaSelect(Expectation,CrossParents);           
            
            %krzyzowanie
            CrossChild = PgaCross(Selection,CrossParents,MutaKids,SubPopulation);            
            
            %mutacja
            MutaChild = PgaMuta(Selection,SubPopulation,MutaKids,(IterationCounter-1)*Workers,CountVal);             
            
            %Odtworzenie populacji
            SubPopulation = [Elite; CrossChild; MutaChild]; 
            SubPopulation = PgaBounds(SubPopulation,ub,lb);
            SubScore = PgaFval(fun,SubPopulation,2,EliteScore);            
%             SubPopulation = [Elite; SubPopulation];
%             SubScore = [EliteScore; SubScore];
            
            Population(:,:,W) = SubPopulation;
            Score(:,W) = SubScore;    
        end
    elseif(mode == 3)
        SubPopulation = Population(:,:,1);
        if(IterationCounter < 2)
            SubScore = Pgafval2(fun,SubPopulation,1);                
        elseif(IterationCounter > 1)
            SubScore = Score(:,1);
        end
        
        %Wyznaczenie Elity
        [EliteScore,EliteIndex] = sort(SubScore); 
        Elite = SubPopulation(EliteIndex(1:EliteKids),:);
        EliteScore = EliteScore(1:EliteKids,1);

        %Przeskalowanie przystosowania
        Expectation = PgaRank(SubScore,EliteIndex,CrossParents);

        %Selekcja do mutacji i krzyzowania
        Selection = PgaSelect(Expectation,CrossParents);           

        %krzyzowanie
        CrossChild = PgaCross(Selection,CrossParents,MutaKids,SubPopulation);            

        %mutacja
        MutaChild = PgaMuta(Selection,SubPopulation,MutaKids);             

        %Odtworzenie populacji
        SubPopulation = [Elite; CrossChild; MutaChild];
        SubPopulation = PgaBounds(SubPopulation,ub,lb);
        SubScore = Pgafval2(fun,SubPopulation,2,EliteScore);            

        for WW = 1 : Workers
            Population(:,:,WW) = SubPopulation;
            Score(:,WW) = SubScore;
        end
    end 
end

function Expectation = PgaRank(Score,EliteIndex,CrossParents)

    Expectation(EliteIndex) = 1 ./ ((1:length(Score))  .^ 0.5);
    Expectation = CrossParents * Expectation ./ sum(Expectation);
    
end

function Score = PgaFval(fun,Population,mode,EliteScore)

    if nargin < 4
        EliteScore = [1; 1];
    end
    [Rows,~] = size(Population);
    Score = zeros(Rows,1);
    if(mode == 1)
        for i = 1 : Rows
            Score(i,1) = feval(fun,Population(i,:));
        end
    elseif(mode == 2)
        [Rows1,~] = size(EliteScore);
        Score(1:Rows1,1) = EliteScore;
        for i = Rows1 + 1 : Rows
            Score(i,1) = feval(fun,Population(i,:));
        end   
    end   
end

function Score = Pgafval2(fun,Population,mode,EliteScore)
    if nargin < 4
        EliteScore = [1; 1];
    end
    [Rows,~] = size(Population);
    Score = zeros(Rows,1);
    if(mode == 1)
        parfor i = 1 : Rows
            Score(i,1) = feval(fun,Population(i,:));
        end
    elseif(mode == 2)
        [Rows1,~] = size(EliteScore);
        Score(1:Rows1,1) = EliteScore;
        parfor i = Rows1 + 1 : Rows
            Score(i,1) = feval(fun,Population(i,:));
        end   
    end   
end

function CrossChild = PgaCross(Selection,CrossParents,MutaKids,Population)

    [~,Columns] = size(Population);
    CrossChild = zeros((CrossParents - MutaKids) / 2,Columns);
    index = 1;
    for i = 1:2:CrossParents - MutaKids - 1                  
        Parent1 = Population(Selection(i),:);
        Parent2 = Population(Selection(i+1),:);
        for j = 1 : Columns
            if(rand < 0.5)
                CrossChild(index,j) = Parent1(1,j);
            else
                CrossChild(index,j) = Parent2(1,j);
            end
        end
        index = index + 1;
    end 


end

function MutaChild = PgaMuta(Selection,Population,MutaKids,sIterationCounter,CountVal)
    
    scale = 1 - (sIterationCounter+labindex)/CountVal;
    [Rows,Columns] = size(Population);
    MutaChild = zeros(MutaKids,Columns);
    for m = 1 : MutaKids
        MutaChild(m,:) = Population(Selection(Rows+m),:) + scale .* randn(1,Columns);
        
        
    end
end

function [Population,Score] = PgaPass(Population,Score,Pass,Workers)

    [~,FirstIndex] = min(Score);
    for i = 1 : Workers
        TempGuy = Population(FirstIndex(i),:,i);
        TempScore = Score(FirstIndex(i),i);
        if(i + Pass <= Workers)
            for j = 1 : Pass
                [~,ThirdIndex] = max(Score); 
                Population(ThirdIndex(i + j),:,i + j) = TempGuy;
                Score(ThirdIndex(i + j),i + j) = TempScore;
            end
        elseif(i + Pass > Workers)
            for k = 1 : (Workers - i)
                [~,ThirdIndex] = max(Score);
                Population(ThirdIndex(i + k),:,i + k) = TempGuy;
                Score(ThirdIndex(i + k),i+k) = TempScore;
            end
            for l = 1 : (Pass - (Workers - i))
                [~,ThirdIndex] = max(Score);
                Population(ThirdIndex(l),:,l) = TempGuy;
                Score(ThirdIndex(l),l) = TempScore;
            end         
        end 
    end
end

function Selection = PgaSelect(Expectation,CrossParents)

    Selection = zeros(CrossParents,1);
    wheel = cumsum(Expectation) / CrossParents;
    for i = 1:CrossParents
        r = rand;
        for j = 1:length(wheel)
            if(r < wheel(j))
                Selection(i) = j;
                break;
            end
        end
    end 

end

function [ExitFlag,StopReason] = PgaStop(IterationCounter,BestScore,CountVal,ToVal,AveChan,Workers)

    FunChange = Inf;
    Window = int8(50/Workers)-1;
    ExitFlag = [];
    StopReason = '';
    if IterationCounter > Window
        Bestfvals =  BestScore((IterationCounter - Window):end);
        FunChange = 0;
        Weight = 0.5;
        for i = 1:Window
            FunChange = FunChange + Weight^(Window-i)*(abs(Bestfvals(i+1) - Bestfvals(i))/(abs(Bestfvals(i))+1));
        end
        FunChange = FunChange/Window;
    end

    if(IterationCounter >= CountVal)
        ExitFlag = 1;
        StopReason ='Iteration Limit Exceeded';
    elseif(BestScore(IterationCounter) <= ToVal)
        ExitFlag = 2;
        StopReason = 'Target Function Value Achieved';
    elseif(FunChange <= AveChan)
        ExitFlag = 3;
        StopReason = 'Minimal Change in Fitness Value Achieved';
    end 
end

function Population = PgaBounds(Population,ub,lb)%mozliweosc ulepszenia nie wszysktie musza miec ograniczenia
    [Rows,Columns] = size(Population);
    if(~isempty(ub) && ~isempty(lb))
        for rr=1:Rows
            bounds = max([Population(rr,:)' - ub(:);lb(:) - Population(rr,:)'],0);
            Population(rr,:) = Population(rr,:) - bounds(1:Columns)';
            Population(rr,:) = Population(rr,:) + bounds(Columns+1:Columns*2)';
        end
    elseif(~isempty(ub) && isempty(lb))
        for rr=1:Rows
            bounds = max(Population(rr,:)' - ub(:),0);
            Population(rr,:) = Population(rr,:) - bounds(1:Columns)';
        end
    elseif(isempty(ub) && ~isempty(lb))
        for rr=1:Rows
            bounds = max(lb(:) - Population(rr,:)',0);
            Population(rr,:) = Population(rr,:) + bounds(Columns+1:Columns*2)';
        end
    end
end

