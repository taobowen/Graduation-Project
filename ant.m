clear all; %清除所有变量
close all; %清图
clc ;      %清屏

global artifacts;
global process;
global timeData;
global e;
global w;
global fragmentNum;
% 每个工件在每道工序上加工时间的信息数组
timeData =  [5 76 74 99 26;
  74 21 83 52 90;
  67 48 6 66 38;
  97 36 71 68 81;
  87 86 64 11 31;
  1 42 20 90 23;
  69 32 99 26 57;
  69 12 54 80 16;
  11 63 24 16 89;
  87 52 43 10 26;
  25 59 88 87 40;
  50 42 72 77 29;
  58 76 71 82 94;
  79 48 20 63 97;
  35 57 78 99 80;
  70 76 53 2 19;
  79 22 77 74 95;
  34 99 49 3 61;
  37 24 32 35 4;
  50 88 46 63 76;]

artifacts = size(timeData,1);    %% 工件个数
process = size(timeData,2);     %% 工序/机器个数
m = 50;    %% m 蚂蚁个数
Alpha = 2;  %% Alpha 表征信息素重要程度的参数 信息启发因子
Beta = 1;  %% Beta 表征启发式因子重要程度的参数 期望启发因子
Rho = 0.2; %% Rho 信息素蒸发系数 信息挥发因子
Q = 6000;     %%信息素增加强度系数
NC_max = 400; %%最大迭代次数
iterations = 10; %% 全局信息素更新迭代周期
unchangeMaxValue = 20;  %% 在进行unchangeMaxValue次迭代后全局最优解没改变，则进行邻域搜索
unchangeNum = 0; %% 记录全局最优解未更新的次数
isDeepUpdate = 1; %% 是否满足邻域搜索条件

% 正交表使用 L(16，4，5)
L_s = 16; % 测试数
L_f = 4;  % 水平数
L_v = 5;  % 同一水平下的片段数
L_table = [1,1,1,1,1;
           1,2,2,2,2;
           1,3,3,3,3;
           1,4,4,4,4;
           2,1,2,3,4;
           2,2,1,4,3;
           2,3,4,1,2;
           2,4,3,2,1;
           3,1,3,4,2;
           3,2,4,3,1;
           3,3,1,2,4;
           3,4,2,1,3;
           4,1,4,2,3;
           4,2,3,1,4;
           4,3,2,4,1;
           4,4,1,3,2;]; %正交表
fragmentNum = artifacts / L_v;         % 每个片段中的工件数
fragment = zeros(L_f,fragmentNum,L_v); % 存放不同水平下的所有片段

% 构造相邻工件i与j之间的完工时间差和等待时间差的信息矩阵
e=zeros(artifacts,artifacts); % 完工时间差
w=zeros(artifacts,artifacts); % 等待时间差
for i=1:artifacts
    for j=1:artifacts
        if i~=j
            diffTime = timechange(i,j);
            e(i,j) = diffTime(1);  %计算第i个工件与第j个工件之间的完工时间差
            w(i,j) = diffTime(2);  %计算第i个工件与第j个工件之间的等待时间差
        end
    end
end

Tau=ones(artifacts,artifacts);     %Tau为信息素矩阵
Tabu=zeros(m,artifacts);   %存储并记录路径的生成
NC=1;               %迭代计数器，记录迭代次数
R_best=zeros(NC_max,artifacts);       %各代最佳路线
L_best=inf.*ones(NC_max,1);   %各代最佳路线的长度
L_ave=zeros(NC_max,1);        %各代路线的平均长度

%% 主要符号说明
%% NC_max 最大迭代次数
%% m 蚂蚁个数
%% artifacts 工件个数
%% process 工序个数
%% Alpha 表征信息素重要程度的参数 信息启发因子
%% Beta 表征启发式因子重要程度的参数 期望启发因子
%% Rho 信息素蒸发系数 信息挥发因子
%% Q 信息素增加强度系数
%% R_best 各代最佳工件序列
%% L_best 各代最佳工件序列的目标函数值
%% Tabu 存储并记录一轮迭代中所有蚂蚁路径的生成
%% Tau 信息素矩阵

while NC<=NC_max        %停止条件之一：达到最大迭代次数，停止
    %%将m只蚂蚁放到n个工件上
    Randpos=[];   %随即存取
    for i=1:(ceil(m/artifacts))
        Randpos=[Randpos,randperm(artifacts)];
    end
    Tabu(:,1)=(Randpos(1,1:m))';
    
    %%m只蚂蚁按概率函数选择下一个工件，完成各自的周游
    for j=2:artifacts     %所在城市不计算
        for i=1:m
            visited=Tabu(i,1:(j-1)); %记录已访问的工件，避免重复访问
            J=zeros(1,(artifacts-j+1));       %待访问的工件
            P=J;                      %待访问工件的选择概率分布
            Jc=1;                     %待访问的工件序号
            r = getThreshold(NC);                 %表示是否选用轮盘赌的随机变量
            % 生成待访问的工件矩阵J
            for k=1:artifacts
                if isempty(find(visited==k, 1))      %开始时置0 如果已访问的工件里没有序号为k的
                    J(Jc)=k;
                    Jc=Jc+1;                         %待访问的工件序号自加1
                end
            end
            %下面计算待选工件的概率分布
            for k=1:length(J)
                P(k)=(Tau(visited(end),J(k))^Alpha)*(Eta(visited(end),J(k))^Beta);
            end
            [~,maxindex] = max(P);
            P=P/(sum(P));
            %按概率原则选取下一个工件
            Pcum=cumsum(P);     %cumsum，元素累加即求和
            if rand > r         %轮盘赌选择下一个工件
                Select = find(Pcum >= rand);
                to_visit=J(Select(1));
            else %取使函数值最大的工件
                to_visit=J(maxindex);
            end
            Tabu(i,j)=to_visit;
        end
    end

    %%记录本次迭代最佳路线
    L=zeros(m,1);     %记录每只蚂蚁选择路线的总完工时间，m*1的列向量
    for i=1:m
        R=Tabu(i,:);
        L(i) = totaltime(R);
    end
    [minvalue,minindex] = min(L);
    model_R_best = Tabu(minindex,:); %此轮蚂蚁搜索过后得出的最佳路线
    model_L_best = minvalue;         %此轮蚂蚁搜索过后得出的最小完工时间

    
    % 正交测试
    
    % 构建不同水平的片段序列
    for i=1:L_v
        fragment(1,:,i) = model_R_best(1,((i-1)*fragmentNum+1):(i*fragmentNum));
        for j=2:(L_f-1)
            basic_fragment_R = fragment((j-1),:,i);
            basic_fragment_L = totaltime(basic_fragment_R);
            new_fragment = localSearch(basic_fragment_R, basic_fragment_L);
            if basic_fragment_L == new_fragment(1,1) % 如果局部搜索没找到更好的解
                randFragment = randperm(fragmentNum);
                index = randFragment(1:2);
                fragment(j,:,i) = insertSerialnumber(basic_fragment_R,index(1),index(2));
            else
                fragment(j,:,i) = new_fragment(1,2:end);
            end
        end
        initial_fragment_R = fragment(1,:,i);
        initial_fragment_L = totaltime(initial_fragment_R);
        last_fragment = localBestSearch(initial_fragment_R, initial_fragment_L);
        if initial_fragment_L == last_fragment(1,1) % 如果局部搜索没找到更好的解
            randFragment=randperm(fragmentNum);
            index=randFragment(1:2);
            fragment(L_f,:,i) = insertSerialnumber(initial_fragment_R,index(1),index(2));
        else
            fragment(L_f,:,i) = last_fragment(1,2:end);
        end
    end
    
    % 将各个水平的序列片段进行重组试验
    L_s_time = zeros(1,L_s);
    for i=1:L_s
        totalFragment = [];
        for j=1:L_v
            totalFragment = [totalFragment fragment(L_table(i,j),:,j)];
        end
        L_s_time(1,i) = totaltime(totalFragment);
    end
    averageFragmentTime = zeros(L_f,L_v);
    averageFragmentTime(1,1) = (L_s_time(1,1) + L_s_time(1,2) + L_s_time(1,3) + L_s_time(1,4)) / L_f;
    averageFragmentTime(1,2) = (L_s_time(1,1) + L_s_time(1,5) + L_s_time(1,9) + L_s_time(1,13)) / L_f;
    averageFragmentTime(1,3) = (L_s_time(1,1) + L_s_time(1,6) + L_s_time(1,11) + L_s_time(1,16)) / L_f;
    averageFragmentTime(1,4) = (L_s_time(1,1) + L_s_time(1,7) + L_s_time(1,12) + L_s_time(1,14)) / L_f;
    averageFragmentTime(1,5) = (L_s_time(1,1) + L_s_time(1,8) + L_s_time(1,10) + L_s_time(1,15)) / L_f;
    averageFragmentTime(2,1) = (L_s_time(1,5) + L_s_time(1,6) + L_s_time(1,7) + L_s_time(1,8)) / L_f;
    averageFragmentTime(2,2) = (L_s_time(1,2) + L_s_time(1,6) + L_s_time(1,10) + L_s_time(1,14)) / L_f;
    averageFragmentTime(2,3) = (L_s_time(1,2) + L_s_time(1,5) + L_s_time(1,12) + L_s_time(1,15)) / L_f;
    averageFragmentTime(2,4) = (L_s_time(1,2) + L_s_time(1,8) + L_s_time(1,11) + L_s_time(1,13)) / L_f;
    averageFragmentTime(2,5) = (L_s_time(1,2) + L_s_time(1,7) + L_s_time(1,9) + L_s_time(1,16)) / L_f;
    averageFragmentTime(3,1) = (L_s_time(1,9) + L_s_time(1,10) + L_s_time(1,11) + L_s_time(1,12)) / L_f;
    averageFragmentTime(3,2) = (L_s_time(1,3) + L_s_time(1,7) + L_s_time(1,11) + L_s_time(1,15)) / L_f;
    averageFragmentTime(3,3) = (L_s_time(1,3) + L_s_time(1,8) + L_s_time(1,9) + L_s_time(1,14)) / L_f;
    averageFragmentTime(3,4) = (L_s_time(1,3) + L_s_time(1,5) + L_s_time(1,10) + L_s_time(1,16)) / L_f;
    averageFragmentTime(3,5) = (L_s_time(1,3) + L_s_time(1,6) + L_s_time(1,12) + L_s_time(1,13)) / L_f;
    averageFragmentTime(4,1) = (L_s_time(1,13) + L_s_time(1,14) + L_s_time(1,15) + L_s_time(1,16)) / L_f;
    averageFragmentTime(4,2) = (L_s_time(1,4) + L_s_time(1,8) + L_s_time(1,12) + L_s_time(1,16)) / L_f;
    averageFragmentTime(4,3) = (L_s_time(1,4) + L_s_time(1,7) + L_s_time(1,10) + L_s_time(1,13)) / L_f;
    averageFragmentTime(4,4) = (L_s_time(1,4) + L_s_time(1,6) + L_s_time(1,9) + L_s_time(1,15)) / L_f;
    averageFragmentTime(4,5) = (L_s_time(1,4) + L_s_time(1,5) + L_s_time(1,11) + L_s_time(1,14)) / L_f;
    
    bestLevel = zeros(1,L_v); % 最优解水平组合
    for i=1:L_v
         [maxvalue,maxindex] = min(averageFragmentTime(:,i));
         bestLevel(1,i) = maxindex;
    end
    model_R_best = [fragment(bestLevel(1,1),:,1) fragment(bestLevel(1,2),:,2) fragment(bestLevel(1,3),:,3) fragment(bestLevel(1,4),:,4) fragment(bestLevel(1,5),:,5)];
    model_L_best = totaltime(model_R_best);
    
    
    % 当前最佳路线和最小目标函数值
    if NC >= 2 && model_L_best >= L_best(NC-1)
        prepare_R_best = R_best(NC-1,:);
        prepare_L_best = L_best(NC-1);
        unchangeNum = unchangeNum + 1;
        
        if isDeepUpdate == 1 && unchangeNum >= unchangeMaxValue               %领域搜索
            PZ = zeros(artifacts*artifacts,3);    %工件序列所有使最大完成时间减小的交换移动集合
            global betterPairNum;
            betterPairNum = 0;
            new_L = 0;

            % 初始化PZ
            for q=1:artifacts
               for p=q:artifacts
                   if q~=p
                        new_R = swapSerialnumber(prepare_R_best,q,p);
                        new_L = totaltime(new_R);
                        if new_L < prepare_L_best %如果新序列比当前最佳序列完工时间更少
                            betterPairNum = betterPairNum + 1;
                            PZ(betterPairNum,1) = q;
                            PZ(betterPairNum,2) = p;
                        end
                   end
               end
            end

            if betterPairNum  > 0
                while 1
                    for j = 1:betterPairNum
                        PZ(j,3) = dependentPointNum(PZ,j);
                    end
                    [maxvalue,maxindex] = max(PZ(1:betterPairNum,3));
                    if maxvalue > 0
                        PZ(maxindex,:) = [];
                        betterPairNum = betterPairNum - 1;
                    else
                        break;
                    end
                end
                for i = 1:betterPairNum
                    prepare_R_best = swapSerialnumber(prepare_R_best,PZ(i,1),PZ(i,2));
                end
                prepare_L_best = totaltime(prepare_R_best);
            else                         % 若此次领域搜索对解的改进没有任何效果，那么以后将停止对该全局最优解变化情况的监视，不再进行领域搜索
                isDeepUpdate = 0; 
            end
        end
    else
        prepare_R_best = model_R_best;
        prepare_L_best = model_L_best;
        unchangeNum = 0;
        isDeepUpdate = 1;
    end
    
    R_best(NC,:)=prepare_R_best; %此轮迭代后的最佳路线
    L_best(NC)=prepare_L_best;   %此轮迭代后的最小目标函数值  
    L_ave(NC)=mean(L);           %此轮迭代后的平均距离
    NC = NC + 1;                 %迭代继续

    %%全局信息素更新
    if mod(NC,iterations) ~= 0
        for i=1:(artifacts-1)
            preparetau = (1-Rho).*Tau(model_R_best(i),model_R_best(i+1))+Rho*(Q/model_L_best);
            Tau(model_R_best(i),model_R_best(i+1)) = tauchange(preparetau); %考虑信息素挥发，更新后的信息素
        end
    else  %每迭代十次进行一次全局信息素更新
        pathset = zeros(artifacts,artifacts);
        for i = 1:(artifacts-1)
            pathset(prepare_R_best(1,i),prepare_R_best(1,i+1)) = 1;
        end
        for i=1:artifacts
            for j=1:artifacts
                if pathset(i,j) ~= 1
                    preparetau = (1-Rho).*Tau(i,j);
                else
                    preparetau = (1-Rho).*Tau(i,j)+Rho*(Q/prepare_L_best);
                end
                Tau(i,j) = tauchange(preparetau);
            end
        end
    end
    
    %%禁忌表清零
    Tabu=zeros(m,artifacts);             %%直到最大迭代次数
end
%%第七步：输出结果
Pos=find(L_best==min(L_best)); %找到最佳路径（非0为真）
Shortest_Route = R_best(Pos(1),:); % 最佳路径
Shortest_Length = L_best(Pos(1)); % 最短距离

figure(1) 
plot(L_best)
xlabel('迭代次数')
ylabel('目标函数值')
title('适应度进化曲线')

figure(2)
plot(L_best)
hold on                         %保持图形
plot(L_ave,'r')
title('平均完工时间和最小完工时间')     %标题


% 输出启发值函数 输入参数为工件序列对<ij>之间的等待时间
function inspire = Eta(i,j)
    global w;
    wt = w(i,j);
    if wt <= 0
        inspire = 2;
    elseif wt >= 200
        inspire = 0.7;
    else
        inspire = 0.7 + 1/wt;
    end
end

% 输出工件i,j之间的完工时间差 输入参数为前后两个工件的序号

function Timedifference = timechange(i,j)
    global process;
    global timeData;
    s = zeros(process,2); % 工件完成时间矩阵
    c = zeros(process,2); % 工件开始时间矩阵
    s(1,1) = 0;
    c(1,1) = s(1,1) + timeData(i,1);
    for k = 2:process
        s(k,1) = c(k-1,1);
        c(k,1) = s(k,1) + timeData(i,k);
    end
    s(1,2) = c(1,1);
    c(1,2) = s(1,2) + timeData(j,1);
    for m = 2:process
        s(m,2) = max(c(m-1,2),c(m,1));
        c(m,2) = s(m,2) + timeData(j,m);
    end
    completeDiff = c(process,2) - c(process,1);
    waitDiff = c(process,2) - timeData(i,1);
    for n = 1:process
        waitDiff = waitDiff - timeData(j,n);
    end
    Timedifference = [completeDiff, waitDiff];
end

% 输出最大完工时间的目标函数值 输入参数为1*n的工件序列矩阵

function Completion = totaltime(serialnumber)
    global process;
    global timeData;
    global e;
    serialLen = length(serialnumber);
    output = 0;
    for i = 1:process
        output = output + timeData(serialnumber(1,1),i);
    end
    for j = 2:serialLen
        output = output + e(serialnumber(1,j-1),serialnumber(1,j));
    end
    Completion = output;
end

% 更新信息素的转化方程
function Tau = tauchange(preparetau)
     if preparetau > 1
         Tau = 1;
     elseif preparetau < 0.1
         Tau = 0.1;
     else
         Tau = preparetau;
     end
end

%获取经过交换操作后的工件序列 输出参数为1*n的工件序列
function newserialnumber = swapSerialnumber(serialnumber,outindex,swapindex)
    mediation = serialnumber(1,outindex);
    serialnumber(1,outindex) = serialnumber(1,swapindex);
    serialnumber(1,swapindex) = mediation;
    newserialnumber = serialnumber;
end

%获取经过插入操作后的工作序列 输出参数为1*n的工件序列
function newserialnumber = insertSerialnumber(serialnumber,outindex,insertindex)
    outNumber = serialnumber(1,outindex);
    if outindex > insertindex
        serialnumber(outindex) = [];
        serialnumber(insertindex+1:end+1) = serialnumber(insertindex:end);
        serialnumber(insertindex) = outNumber;
    else
        serialnumber(insertindex+1:end+1) = serialnumber(insertindex:end);
        serialnumber(insertindex) = outNumber;
        serialnumber(outindex) = [];
    end
    newserialnumber = serialnumber;
end

% 求PZ中与第i个移动不相互独立的移动的个数
function noindependentpoint = dependentPointNum(PZ,i)
    global betterPairNum;
    output = 0;
    for j = 1:betterPairNum
        if j~=i && ((PZ(i,1)-PZ(j,1))^2 <= 1 || (PZ(i,2)-PZ(j,1))^2 <= 1 || (PZ(i,1)-PZ(j,2))^2  <= 1 || (PZ(i,2)-PZ(j,2))^2 <= 1)
            output = output + 1;
        end
    end
    noindependentpoint = output;
end

% 输出轮盘赌的阈值 输入参数为迭代次数
function randomThreshold = getThreshold(NC)
    if NC < 200
        randomThreshold = 0.1 + NC/250;
    else
        randomThreshold = 0.9;
    end
end

% 搜索局部更优解
function localSolution = localSearch(model_R_best, model_L_best)
    global fragmentNum;
    for q=1:fragmentNum
       for p=q:fragmentNum
           isBreak = 0;
           if q~=p
                new_R = insertSerialnumber(model_R_best,q,p);
                new_L = totaltime(new_R);
                if new_L < model_L_best   %如果新序列比当前最佳序列完工时间更少
                    localSolution = [new_L new_R]; %更新最佳路线
                    isBreak = 1;
                    break;
                end
           end
       end
       if isBreak == 1
           break;
       end
    end
    if isBreak == 0
        localSolution = [model_L_best model_R_best];
    end
end

% 搜索局部最优解
function localBestSolution = localBestSearch(model_R_best, model_L_best)
    global fragmentNum;
    prepare_R_best = model_R_best;
    prepare_L_best = model_L_best;
    for q=1:fragmentNum
       for p=q:fragmentNum
           if q~=p
                new_R = swapSerialnumber(model_R_best,q,p);
                new_L = totaltime(new_R);
                if new_L < prepare_L_best   %如果新序列比当前最佳序列完工时间更少
                    prepare_R_best = new_R;
                    prepare_L_best = new_L;
                end
           end
       end
    end
    localBestSolution = [prepare_L_best prepare_R_best]; %更新最佳路线
end








