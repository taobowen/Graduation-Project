clear all; %������б���
close all; %��ͼ
clc ;      %����

global artifacts;
global process;
global timeData;
global e;
global w;
global fragmentNum;
% ÿ��������ÿ�������ϼӹ�ʱ�����Ϣ����
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

artifacts = size(timeData,1);    %% ��������
process = size(timeData,2);     %% ����/��������
m = 50;    %% m ���ϸ���
Alpha = 2;  %% Alpha ������Ϣ����Ҫ�̶ȵĲ��� ��Ϣ��������
Beta = 1;  %% Beta ��������ʽ������Ҫ�̶ȵĲ��� ������������
Rho = 0.2; %% Rho ��Ϣ������ϵ�� ��Ϣ�ӷ�����
Q = 6000;     %%��Ϣ������ǿ��ϵ��
NC_max = 400; %%����������
iterations = 10; %% ȫ����Ϣ�ظ��µ�������
unchangeMaxValue = 20;  %% �ڽ���unchangeMaxValue�ε�����ȫ�����Ž�û�ı䣬�������������
unchangeNum = 0; %% ��¼ȫ�����Ž�δ���µĴ���
isDeepUpdate = 1; %% �Ƿ�����������������

% ������ʹ�� L(16��4��5)
L_s = 16; % ������
L_f = 4;  % ˮƽ��
L_v = 5;  % ͬһˮƽ�µ�Ƭ����
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
           4,4,1,3,2;]; %������
fragmentNum = artifacts / L_v;         % ÿ��Ƭ���еĹ�����
fragment = zeros(L_f,fragmentNum,L_v); % ��Ų�ͬˮƽ�µ�����Ƭ��

% �������ڹ���i��j֮����깤ʱ���͵ȴ�ʱ������Ϣ����
e=zeros(artifacts,artifacts); % �깤ʱ���
w=zeros(artifacts,artifacts); % �ȴ�ʱ���
for i=1:artifacts
    for j=1:artifacts
        if i~=j
            diffTime = timechange(i,j);
            e(i,j) = diffTime(1);  %�����i���������j������֮����깤ʱ���
            w(i,j) = diffTime(2);  %�����i���������j������֮��ĵȴ�ʱ���
        end
    end
end

Tau=ones(artifacts,artifacts);     %TauΪ��Ϣ�ؾ���
Tabu=zeros(m,artifacts);   %�洢����¼·��������
NC=1;               %��������������¼��������
R_best=zeros(NC_max,artifacts);       %�������·��
L_best=inf.*ones(NC_max,1);   %�������·�ߵĳ���
L_ave=zeros(NC_max,1);        %����·�ߵ�ƽ������

%% ��Ҫ����˵��
%% NC_max ����������
%% m ���ϸ���
%% artifacts ��������
%% process �������
%% Alpha ������Ϣ����Ҫ�̶ȵĲ��� ��Ϣ��������
%% Beta ��������ʽ������Ҫ�̶ȵĲ��� ������������
%% Rho ��Ϣ������ϵ�� ��Ϣ�ӷ�����
%% Q ��Ϣ������ǿ��ϵ��
%% R_best ������ѹ�������
%% L_best ������ѹ������е�Ŀ�꺯��ֵ
%% Tabu �洢����¼һ�ֵ�������������·��������
%% Tau ��Ϣ�ؾ���

while NC<=NC_max        %ֹͣ����֮һ���ﵽ������������ֹͣ
    %%��mֻ���Ϸŵ�n��������
    Randpos=[];   %�漴��ȡ
    for i=1:(ceil(m/artifacts))
        Randpos=[Randpos,randperm(artifacts)];
    end
    Tabu(:,1)=(Randpos(1,1:m))';
    
    %%mֻ���ϰ����ʺ���ѡ����һ����������ɸ��Ե�����
    for j=2:artifacts     %���ڳ��в�����
        for i=1:m
            visited=Tabu(i,1:(j-1)); %��¼�ѷ��ʵĹ����������ظ�����
            J=zeros(1,(artifacts-j+1));       %�����ʵĹ���
            P=J;                      %�����ʹ�����ѡ����ʷֲ�
            Jc=1;                     %�����ʵĹ������
            r = getThreshold(NC);                 %��ʾ�Ƿ�ѡ�����̶ĵ��������
            % ���ɴ����ʵĹ�������J
            for k=1:artifacts
                if isempty(find(visited==k, 1))      %��ʼʱ��0 ����ѷ��ʵĹ�����û�����Ϊk��
                    J(Jc)=k;
                    Jc=Jc+1;                         %�����ʵĹ�������Լ�1
                end
            end
            %��������ѡ�����ĸ��ʷֲ�
            for k=1:length(J)
                P(k)=(Tau(visited(end),J(k))^Alpha)*(Eta(visited(end),J(k))^Beta);
            end
            [~,maxindex] = max(P);
            P=P/(sum(P));
            %������ԭ��ѡȡ��һ������
            Pcum=cumsum(P);     %cumsum��Ԫ���ۼӼ����
            if rand > r         %���̶�ѡ����һ������
                Select = find(Pcum >= rand);
                to_visit=J(Select(1));
            else %ȡʹ����ֵ���Ĺ���
                to_visit=J(maxindex);
            end
            Tabu(i,j)=to_visit;
        end
    end

    %%��¼���ε������·��
    L=zeros(m,1);     %��¼ÿֻ����ѡ��·�ߵ����깤ʱ�䣬m*1��������
    for i=1:m
        R=Tabu(i,:);
        L(i) = totaltime(R);
    end
    [minvalue,minindex] = min(L);
    model_R_best = Tabu(minindex,:); %����������������ó������·��
    model_L_best = minvalue;         %����������������ó�����С�깤ʱ��

    
    % ��������
    
    % ������ͬˮƽ��Ƭ������
    for i=1:L_v
        fragment(1,:,i) = model_R_best(1,((i-1)*fragmentNum+1):(i*fragmentNum));
        for j=2:(L_f-1)
            basic_fragment_R = fragment((j-1),:,i);
            basic_fragment_L = totaltime(basic_fragment_R);
            new_fragment = localSearch(basic_fragment_R, basic_fragment_L);
            if basic_fragment_L == new_fragment(1,1) % ����ֲ�����û�ҵ����õĽ�
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
        if initial_fragment_L == last_fragment(1,1) % ����ֲ�����û�ҵ����õĽ�
            randFragment=randperm(fragmentNum);
            index=randFragment(1:2);
            fragment(L_f,:,i) = insertSerialnumber(initial_fragment_R,index(1),index(2));
        else
            fragment(L_f,:,i) = last_fragment(1,2:end);
        end
    end
    
    % ������ˮƽ������Ƭ�ν�����������
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
    
    bestLevel = zeros(1,L_v); % ���Ž�ˮƽ���
    for i=1:L_v
         [maxvalue,maxindex] = min(averageFragmentTime(:,i));
         bestLevel(1,i) = maxindex;
    end
    model_R_best = [fragment(bestLevel(1,1),:,1) fragment(bestLevel(1,2),:,2) fragment(bestLevel(1,3),:,3) fragment(bestLevel(1,4),:,4) fragment(bestLevel(1,5),:,5)];
    model_L_best = totaltime(model_R_best);
    
    
    % ��ǰ���·�ߺ���СĿ�꺯��ֵ
    if NC >= 2 && model_L_best >= L_best(NC-1)
        prepare_R_best = R_best(NC-1,:);
        prepare_L_best = L_best(NC-1);
        unchangeNum = unchangeNum + 1;
        
        if isDeepUpdate == 1 && unchangeNum >= unchangeMaxValue               %��������
            PZ = zeros(artifacts*artifacts,3);    %������������ʹ������ʱ���С�Ľ����ƶ�����
            global betterPairNum;
            betterPairNum = 0;
            new_L = 0;

            % ��ʼ��PZ
            for q=1:artifacts
               for p=q:artifacts
                   if q~=p
                        new_R = swapSerialnumber(prepare_R_best,q,p);
                        new_L = totaltime(new_R);
                        if new_L < prepare_L_best %��������бȵ�ǰ��������깤ʱ�����
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
            else                         % ���˴����������Խ�ĸĽ�û���κ�Ч������ô�Ժ�ֹͣ�Ը�ȫ�����Ž�仯����ļ��ӣ����ٽ�����������
                isDeepUpdate = 0; 
            end
        end
    else
        prepare_R_best = model_R_best;
        prepare_L_best = model_L_best;
        unchangeNum = 0;
        isDeepUpdate = 1;
    end
    
    R_best(NC,:)=prepare_R_best; %���ֵ���������·��
    L_best(NC)=prepare_L_best;   %���ֵ��������СĿ�꺯��ֵ  
    L_ave(NC)=mean(L);           %���ֵ������ƽ������
    NC = NC + 1;                 %��������

    %%ȫ����Ϣ�ظ���
    if mod(NC,iterations) ~= 0
        for i=1:(artifacts-1)
            preparetau = (1-Rho).*Tau(model_R_best(i),model_R_best(i+1))+Rho*(Q/model_L_best);
            Tau(model_R_best(i),model_R_best(i+1)) = tauchange(preparetau); %������Ϣ�ػӷ������º����Ϣ��
        end
    else  %ÿ����ʮ�ν���һ��ȫ����Ϣ�ظ���
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
    
    %%���ɱ�����
    Tabu=zeros(m,artifacts);             %%ֱ������������
end
%%���߲���������
Pos=find(L_best==min(L_best)); %�ҵ����·������0Ϊ�棩
Shortest_Route = R_best(Pos(1),:); % ���·��
Shortest_Length = L_best(Pos(1)); % ��̾���

figure(1) 
plot(L_best)
xlabel('��������')
ylabel('Ŀ�꺯��ֵ')
title('��Ӧ�Ƚ�������')

figure(2)
plot(L_best)
hold on                         %����ͼ��
plot(L_ave,'r')
title('ƽ���깤ʱ�����С�깤ʱ��')     %����


% �������ֵ���� �������Ϊ�������ж�<ij>֮��ĵȴ�ʱ��
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

% �������i,j֮����깤ʱ��� �������Ϊǰ���������������

function Timedifference = timechange(i,j)
    global process;
    global timeData;
    s = zeros(process,2); % �������ʱ�����
    c = zeros(process,2); % ������ʼʱ�����
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

% �������깤ʱ���Ŀ�꺯��ֵ �������Ϊ1*n�Ĺ������о���

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

% ������Ϣ�ص�ת������
function Tau = tauchange(preparetau)
     if preparetau > 1
         Tau = 1;
     elseif preparetau < 0.1
         Tau = 0.1;
     else
         Tau = preparetau;
     end
end

%��ȡ��������������Ĺ������� �������Ϊ1*n�Ĺ�������
function newserialnumber = swapSerialnumber(serialnumber,outindex,swapindex)
    mediation = serialnumber(1,outindex);
    serialnumber(1,outindex) = serialnumber(1,swapindex);
    serialnumber(1,swapindex) = mediation;
    newserialnumber = serialnumber;
end

%��ȡ�������������Ĺ������� �������Ϊ1*n�Ĺ�������
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

% ��PZ�����i���ƶ����໥�������ƶ��ĸ���
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

% ������̶ĵ���ֵ �������Ϊ��������
function randomThreshold = getThreshold(NC)
    if NC < 200
        randomThreshold = 0.1 + NC/250;
    else
        randomThreshold = 0.9;
    end
end

% �����ֲ����Ž�
function localSolution = localSearch(model_R_best, model_L_best)
    global fragmentNum;
    for q=1:fragmentNum
       for p=q:fragmentNum
           isBreak = 0;
           if q~=p
                new_R = insertSerialnumber(model_R_best,q,p);
                new_L = totaltime(new_R);
                if new_L < model_L_best   %��������бȵ�ǰ��������깤ʱ�����
                    localSolution = [new_L new_R]; %�������·��
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

% �����ֲ����Ž�
function localBestSolution = localBestSearch(model_R_best, model_L_best)
    global fragmentNum;
    prepare_R_best = model_R_best;
    prepare_L_best = model_L_best;
    for q=1:fragmentNum
       for p=q:fragmentNum
           if q~=p
                new_R = swapSerialnumber(model_R_best,q,p);
                new_L = totaltime(new_R);
                if new_L < prepare_L_best   %��������бȵ�ǰ��������깤ʱ�����
                    prepare_R_best = new_R;
                    prepare_L_best = new_L;
                end
           end
       end
    end
    localBestSolution = [prepare_L_best prepare_R_best]; %�������·��
end








