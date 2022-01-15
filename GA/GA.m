function GA()
    clear;clc;
%% 1.������ʼ��
    N=48;                     %���еĸ���
    M=100;                    %��Ⱥ�ĸ���
    C=500;                    %��������
    m=2;                      %��Ӧֵ��һ����̭����ָ��
    Pc=0.5;                   %�������
    Pmutation=0.2;            %�������
    X=[6734	2233	5530	401	3082	7608	7573	7265	6898	1112	5468	5989	4706	4612	6347	6107	7611	7462	7732	5900	4483	6101	5199	1633	4307	675	7555	7541	3177	7352	7545	3245	6426	4608	23	7248	7762	7392	3484	6271	4985	1916	7280	7509	10	6807	5185	3023
];
    Y=[1453	10	1424	841	1644	4458	3716	1268	1885	2049	2606	2873	2674	2035	2683	669	5184	3590	4723	3561	3369	1110	2182	2809	2322	1006	4819	3981	756	4506	2801	3305	3173	1198	2216	3779	4595	2244	2829	2135	140	1569	4899	3239	2676	2993	3258	1942
];
    city=[X',Y'];
    D=zeros(N,N);%���ݳ����������ɵ��ڽӾ���
for   i=1:N
    for j=1:N
         dis=[(city(i,1)-city(j,1))^2+(city(i,2)-city(j,2))^2];
          D(i,j)=dis^(0.5);
           D(j,i)= D(i,j);
    end
end

%% 2.1���ɳ�ʼȺ��
    popm=zeros(M,N);
    for i=1:M
        popm(i,:)=randperm(N);             %��N����������
    end
    %%2.2���ѡ��һ����Ⱥ
    R=popm(1,:);
    
    figure(1);
    scatter(X,Y,'ro'); 
    axis([0 8000 0 6000]);
    figure(2);
    plot_route(city,R);                     %������Ⱥ������֮�������
    axis([0 8000 0 6000]);
    title('48��������������·��ͼ');
    %% 3.��ʼ����Ⱥ������Ӧ����
    fitness=zeros(M,1);
    len=zeros(M,1);
    for i=1:M
        len(i,1)=myLength(D,popm(i,:));
    end
    maxlen=max(len);
    minlen=min(len);
    fitness=fit(len,m,maxlen,minlen);
    rr=find(len==minlen);
    R=popm(rr(1,1),:);
    for i=1:N
        fprintf('%d ',R(i));
    end
    fprintf('\n');
    fitness=fitness/sum(fitness);

    distance_min=zeros(C+1,1);  %%���ε�������С����Ⱥ�ľ���
    while C>=0
        fprintf('������%d��\n',C);
        %%ѡ�����
        nn=0;
        for i=1:size(popm,1)
            len_1(i,1)=myLength(D,popm(i,:));
            jc=rand*0.3;
            for j=1:size(popm,1)
                if fitness(j,1)>=jc
                    nn=nn+1;
                    popm_sel(nn,:)=popm(j,:);
                    break;
                end
            end
        end
        %%ÿ��ѡ�񶼱������ŵ���Ⱥ
        popm_sel=popm_sel(1:nn,:);
        [len_m,len_index]=min(len_1);
        popm_sel=[popm_sel;popm(len_index,:)];

        %%�������
        nnper=randperm(nn);
        A=popm_sel(nnper(1),:);
        B=popm_sel(nnper(2),:);
        for i=1:nn*Pc
            [A,B]=cross(A,B);
            popm_sel(nnper(1),:)=A;
            popm_sel(nnper(2),:)=B;
        end
        %%�������
        for i=1:nn
            pick=rand;
            while pick==0
                pick=rand;
            end
            if pick<=Pmutation
                popm_sel(i,:)=Mutation(popm_sel(i,:));
            end
        end
        %%����Ӧ�Ⱥ���
        NN=size(popm_sel,1);
        len=zeros(NN,1);
        for i=1:NN
            len(i,1)=myLength(D,popm_sel(i,:));
        end
        maxlen=max(len);
        minlen=min(len);
        distance_min(C+1,1)=minlen;
        fitness=fit(len,m,maxlen,minlen);
        rr=find(len==minlen);
        fprintf('minlen=%d\n',minlen);
        R=popm_sel(rr(1,1),:);
        for i=1:N
            fprintf('%d ',R(i));
        end
        fprintf('\n');
        popm=[];
        popm=popm_sel;
        C=C-1;
        %pause(1);
    end
    figure(3)
    plot_route(city,R);
    axis([0 8000 0 6000]);
    title('48����������������Ż�·��ͼ');
    text(1000,5000,str);
    figure(4)
    plot(distance_min);
end
%% ���е������
function plot_route(a,R)
    scatter(a(:,1),a(:,2),'rx');
    hold on;
    plot([a(R(1),1),a(R(length(R)),1)],[a(R(1),2),a(R(length(R)),2)]);
    hold on;
    for i=2:length(R)
        x0=a(R(i-1),1);
        y0=a(R(i-1),2);
        x1=a(R(i),1);
        y1=a(R(i),2);
        xx=[x0,x1];
        yy=[y0,y1];
        plot(xx,yy);
        hold on;
    end
end




