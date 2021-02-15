% author: Eman AL-Hawri
%% define problem variables and constants.
%read network file
readNetFile;

nodes=unique(network);
% Select the gateways nodes: you can use either single or multiple gateways
gateways=[8,20];

links(:,1)=1:numberOfLinks;
netSize=numberOfLinks*numberOfNodes;
networkF=[network links];
[netRows, netCols]=size(network);
numberOfGetways=length(gateways);
numberOfsources=numberOfNodes-numberOfGetways;
decisionVar=zeros(numberOfNodes+numberOfLinks+numberOfLinks+(numberOfsources*numberOfLinks)+numberOfNodes+numberOfsources*numberOfGetways,1);
numberOfVariables=length(decisionVar);
Aeq=zeros(numberOfNodes*numberOfsources,numberOfVariables);
beq=zeros(numberOfNodes*numberOfsources,1);
A=zeros(numberOfNodes+numberOfLinks+numberOfsources*numberOfLinks,numberOfVariables);
b=zeros(numberOfNodes+numberOfLinks+numberOfsources*numberOfLinks,1);
numberOfHops=3;

% Exclude gateways from the nodes.
sources=nodes;
for g=1:length(gateways)
    sources=sources(sources~=gateways(g));
end

%% 1- define the lower and the upper bounds
lb=zeros(numberOfVariables,1);
ub=ones(numberOfVariables,1);

%% Only sources and their links
% if the nodes is selected to be a gateway nodesOutG of this node will 1
% and 0 otherwize.
nodesGateways=logical(nodes==gateways);
nodesOutG=zeros(numberOfNodes,1);
for n=1:numberOfNodes
nodesOutG(n)=sum(nodesGateways(n,:),2);
end

%  exctract all links except the links whose their sources are the gateways
linksOutG=links;
for n=1:numberOfNodes
    for i=1:numberOfLinks
        if networkF(i,1)==n && nodesOutG(n)==1
            linksOutG=linksOutG(linksOutG~=networkF(i,3));
        end
    end
end
%% Draw the graph
gatewayIsSrc=zeros(length(gateways),numberOfLinks);
srcs=network(:,1);
dests=network(:,2);
G=digraph(srcs,dests);
%G=graph(srcs,dests);
 p=plot(G);

% put one when a gateway is a source for links
for i=1:length(gateways)
    for j=1:numberOfLinks
        if(network(j,1)==gateways(i) )
            gatewayIsSrc(gateways(i),j)=1;
        else
            gatewayIsSrc(gateways(i),j)=0;
        end
    end
end
%% This for data to flow between nodes and gateways,using tree-based routing
% costraints: 1- conservation law
isSource=zeros(numberOfNodes,numberOfLinks);
isDest=zeros(numberOfNodes, numberOfLinks);

% check which nodes are source for which links
e=0;
for k=1:numberOfsources
    for i=1:numberOfNodes
        for g=1:numberOfGetways
            for j=1:numberOfLinks
                bool=isAgateway(j,1,gateways,network);
                if(network(j,1)==i )   &&  bool==false
                    isSource(i+e,j)=1;
                else
                    isSource(i+e,j)=0;
                end            
            end
        end
    end
    e=e+numberOfNodes;
end

% check which nodes are distination for which links.
e=0;
for k=1:length(sources)
    for i=1:numberOfNodes
        for g=1:numberOfGetways
            for j=1:numberOfLinks
                bool=isAgateway(j,1,gateways,network);
                if(network(j,2)==i) && network(j,2)~= sources(k) && bool== false
                    isDest(i+e,j)=-1;
                else
                    isDest(i+e,j)=0;
                end
            end
        end
    end
            e=e+numberOfNodes;
end

% put both isSource and IsDestination matrices in one matrix
adjMatrix=zeros(numberOfNodes*numberOfsources, numberOfLinks);
adjMatrix(:,:)=isSource(:,:)+isDest(:,:);

%% fill Aeq matrix with first constraint values
j=numberOfNodes+1;
jump1=numberOfLinks+numberOfNodes;
m=0;
adjM=0;
    for kk=1:numberOfsources
        for i=1:numberOfNodes
            Aeq(i+m,j:jump1)=adjMatrix(i+adjM,:);
        end
        j=jump1+1;
        jump1=jump1+numberOfLinks;
        m=m+numberOfNodes;
        adjM=adjM+numberOfNodes;
    end


% Fill beq for the first constraints if source the sum=1 dest the sum=-1
 % and 0 otherwise(intermediate nodes)
 e=0;
    for mm=1: length(sources)
        for nn= 1:numberOfNodes
            
            if nn==sources(mm)
                beq(nn+e)= 1;
                
            else
                beq(nn+e)=0;
            end
        end
        e=e+numberOfNodes;
    end

% fill Aeq when n= a gateway
jump1=numberOfLinks+numberOfNodes;
betaCols=numberOfNodes+numberOfLinks+numberOfLinks+(numberOfsources*numberOfLinks)+numberOfNodes;
d=1;
a=0;
m=0;
    for kk=1:numberOfsources
        for i=1:numberOfNodes
            if i==sources(kk)
                Aeq(i+m,betaCols+d)=0;
            elseif nodesOutG(i)==1
                Aeq(i+m,betaCols+d+a)=1;
                beq(i+m)=0;
                d=d+1;
            else
                Aeq(i+m,betaCols+d)=0;
            end
        end
        a=a+numberOfGetways;
        d=1;
        m=m+numberOfNodes;
    end

%% every source should select only one gateway.

betaRow=numberOfNodes*numberOfsources;
betaS=betaCols+1;
betacolsEnd=betaS+numberOfGetways-1;

for s=1:numberOfsources
    Aeq(betaRow+s, betaS:betacolsEnd)=1;
    beq(betaRow+s)=1;
    betaS=betacolsEnd+1;
    betacolsEnd=betacolsEnd+numberOfGetways;
    
end
%% calculate how many links are left from each source
startFrom=zeros(1,length(nodes));
for s=1:length(nodes)   
    startFrom(1,s)=sum(isSource(s,:)==1);
end
%%  calculte how many links are entered to each destination
arriveAt=zeros(1,length(nodes));

for d=1:length(nodes)
    arriveAt(1,d)=sum(isDest(d,:)==-1);
end

%% Inequality constraints
% 1- ensure that trees are being built, if a link is performing data forwarding

ineqSl=numberOfNodes+1;
ineqL=numberOfsources * numberOfLinks + numberOfNodes+1;
e=0;
i=0;
g=1;
for lC1=1:numberOfLinks  
    for ineqrows=1:length(sources) 
           bool=isAgateway(lC1,1,gateways,network);
            if  bool ==true
                A(ineqrows+e, ineqSl+i)=0;% sl
                 A(ineqrows+e,ineqL)=0;%yl
            else
                A(ineqrows+e, ineqSl+i)=1;% sl
                A(ineqrows+e,ineqL)=-1;%yl
                ineqSl=ineqSl+numberOfLinks;
             end
    end
    i=i+1;
    ineqL=ineqL+1;
    e=e+numberOfsources;
    ineqSl=numberOfNodes+1;
end

% fill b
b(1:numberOfsources*numberOfLinks,1)=0;

%2- ensure that no other link with the same source can be used, node has a single parent.
rowLoop=(numberOfLinks*numberOfsources)+1;
e=1;
q=0;
tmpf=[];
for i=1:numberOfNodes
    if nodesOutG(i)==1
        q=q+sum(gatewayIsSrc(i,:)==1);
        as=zeros(sum(gatewayIsSrc(i,:)==1),numberOfsources* numberOfLinks);
        tmpf=[tmpf;as];
        continue;
    end
    as=zeros(startFrom(i),numberOfLinks);
    as(e:startFrom(i),e+q:q+startFrom(i))=1;
    
    for k=1:startFrom(i)
        as(k,k+q)=0;
    end
    tmp=repmat(as,1,numberOfsources);
    tmpf=[tmpf;tmp];
    q=q+startFrom(i);
end
colS=1+numberOfNodes;
colE=numberOfsources*numberOfLinks+numberOfNodes;
A(rowLoop:numberOfLinks+rowLoop-1,colS:colE)=tmpf(:,:);

% fill b with Yl variable
c3Start=rowLoop-1;
c3End=c3Start+numberOfLinks-1;
ylStart=numberOfsources * numberOfLinks+(numberOfNodes);

for j=1: numberOfLinks
    bool=isAgateway(j,1,gateways,network);
    if bool==true
        A(c3Start+j,ylStart+j)=0;
        b(c3Start+j)=0;   
    else
        A(c3Start+j,ylStart+j)=netSize;
        b(c3Start+j)=netSize;
    end
end

%3- limit the number of hops for any flow (tree depth is limited)
 rowsC3=rowLoop+numberOfLinks;
 slC3=numberOfNodes+1;
     for i=1:numberOfsources
         A(rowsC3,slC3:slC3+(numberOfLinks)-1)=1;
         rowsC3=rowsC3+1;
         slC3=slC3+numberOfLinks;
     end
 
 %fill b with the number of hops
 rowsC3b=rowLoop+numberOfLinks;
     for i=1:numberOfsources
         b(rowsC3b)=numberOfHops;
         rowsC3b=rowsC3b+1;
     end

%% objective function
% minimizes the total number of nodes performing encoding
f=zeros(numberOfVariables,1);
f(1:numberOfHops,1)=1;
IntCon=1:numberOfVariables;

% Solve the mathematical model
[x] = intlinprog(f,IntCon,A,b,Aeq,beq,lb,ub);

% show the results
if isempty(x)
    return;
else
resu = cell(numberOfVariables,1);
names = cell(numberOfVariables,1);
for i = 1:numberOfNodes
    names{i} = ['Ne' '_' num2str(i)  ];
    resu{i} =    [names{i} ':  ' num2str(x(i))];
end
k=numberOfNodes;
w=0;
for src=1:length(sources)
    for i = numberOfNodes+1:numberOfLinks+numberOfNodes
        
        names{i+w} = ['delta' '_' num2str(sources(src)) '_' num2str(i-k)];
        resu{i+w} =    [names{i+w} ':  ' num2str(x(i)+w)];  
    end
    w=w+numberOfLinks;
end

k=numberOfsources *numberOfLinks+numberOfNodes;
for i=k+1 :k+numberOfLinks
    names{i}=['gamma' '_' num2str(i-k)];
    resu{i} =    [names{i} ':  ' num2str(x(i))];
end

wF=0;
k=numberOfsources *numberOfLinks+numberOfNodes+numberOfLinks;
for i=k+1: k+numberOfLinks
    names{i+wF}=['sigma' '_' num2str(i-k)];
    resu{i+wF}= [names{i+wF} ': ' num2str(x(i)+wF)];
end

k=numberOfsources* numberOfLinks + numberOfNodes+numberOfLinks+numberOfLinks;
for i=k+1:k+numberOfNodes
    names{i}=['eta' '_' num2str(i-k)];
    resu{i}=[names{i} ': ' num2str(x(i))];
end
w=0;
k=numberOfsources* numberOfLinks + numberOfNodes+numberOfLinks+numberOfLinks+numberOfNodes;
for s=1:numberOfsources
    for i=k+1:k+numberOfGetways
        names{i+w}=['beta' '_' num2str(sources(s)) '_'  num2str(gateways(i-k))];
        resu{i+w}=[names{i+w} ':' num2str(x(i)+w)];
        
    end
     w=w+numberOfGetways;
end
ResTable1=table(names,x);
end
%% Optimization tree links gammas
if isempty(x)
    gammaTree=0;
    sigmaTree=0;
    return;
else
    e=0;
    allDeltas=x(numberOfNodes+1:numberOfNodes+numberOfsources*numberOfLinks);
    deltas=zeros(numberOfsources,numberOfLinks);
    for i=1:numberOfsources
    deltas(i,1:numberOfLinks)=x(numberOfNodes+1+e:numberOfNodes+e+numberOfLinks);
    e=e+numberOfLinks;
    end
   
    gammaTree=x(numberOfsources*numberOfLinks+numberOfNodes+1:numberOfsources*numberOfLinks+numberOfNodes+numberOfLinks);
    gammaTree=logical(gammaTree);
    activeLinks=sum(logical(gammaTree));
    Opt1srcOptT=[];%zeros(activeLinks,1);
    Op1destOptT=[];%zeros(activeLinks,1);

    t=1;
    for l=1:numberOfLinks
        if  gammaTree(l)==1
            Opt1srcOptT(t)=network(l,1);
            Op1destOptT(t)=network(l,2);
            
            t=t+1;
            if t>activeLinks
                break;
            end
        end
    end
    
    % Color the graph's links in the solution path
    p.LineWidth=1.7;
    p.ArrowSize=12;
    highlight(p,Opt1srcOptT,Op1destOptT,'EdgeColor','m', 'LineStyle',':');
    highlight(p,gateways,'NodeColor',[0 0.75 0]);
    p.MarkerSize=5;
end
axis off;

%% Functions Section
function bool =isAgateway(i,e,gateways,network)
numberOfGetways=length(gateways);
    for g=1:numberOfGetways
        if network(i,e)==gateways(g)
            bool=true;
            break;
        else
            bool=false;
        end
    end
end

