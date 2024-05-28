clc;clear all;close all;
%    ////////////////////////     INPUT  ///////////////////////////////
printf("Numerical analysis assignment by Nasos Karras(dit20082@go.uop.gr)\nINTEGRAL CALCULATION \n");
printf("Values of x \nx=[1,3]\n");

printf("Choose which function you want to integrate\n")
printf("1. sin(x)\n2. x^3 - 6x^2 + 4x + 12\n3. 2^x \n4.Skip\n")
choice=input("Choose:");

%creating our function via function handles
switch choice
    case 1
        fun=@sin;  % function handle for sin(x) , @(varUsed) function(πx log(x))
    case 2
        fun=@(x) x.^3-6*(x.^2)+4*x+12;
    case  3
        fun=@(x) 2.^x;    
    otherwise
        fun=@sin;
endswitch

%initializing x and y
x=[1,3];y=[fun(1),fun(3)];

%         //////////////////      trapezoid function declaration   ////////////////////////////

function ftrap(x,y,fun)  

%trapezoid formula
%error and relative  error calculations are inside the printf
resultTrapeziou=( x(2)-x(1) ) * ( (y(1)+y(2)) / 2 );
    
    printf("resultTrapeziou %d\n",resultTrapeziou);
    printf("True intergral result %d\n",integral(fun,1,3));
    printf("Error %d\n",error=integral(fun,1,3)-resultTrapeziou);
    printf("Relative error %d %%\n",errorRelative=abs(error/integral(fun,1,3))*100);

endfunction

%         //////////////////    composite  trapezoid function declaration   ////////////////////////////

function ftrapSteps(x,y,fun)  

    stepSum=0;
    step=input("\n THE STEP IS ?\n");
    for i = 1:stepSum-1
        stepSum=stepSum+fun( x(1)+ i * ( (x(2)-x(1))/step ) );
    endfor
                            % h= b-a / n ^^^^ this small misread cost me 1.5 hours >.<

    resultTrapeziou=( ( x(2)-x(1) )/( 2 * step ) )*( y(1)+ 2*stepSum +y(2) );
        
        printf("resultTrapeziou %d\n",resultTrapeziou);
        printf("True intergral result %d\n",integral(fun,1,3));
        printf("Error %d\n",error=integral(fun,1,3)-resultTrapeziou);
        printf("Relative error %d %%\n",errorRelative=abs(error/integral(fun,1,3))*100);

endfunction

%         //////////////////    simpsons 1/3 function declaration   ////////////////////////////

      
function simp1(x,y,fun)  
h=(x(2)-x(1))/2; %  h=B-A/2
%simp1/3 formula
resultSimp=h/3*(y(1)+4*fun((x(2)+x(1))/2)+y(2));
    
    printf("resultSimp1/3 %d\n",resultSimp);
    printf("True intergral result %d\n",integral(fun,x(1),x(2)));
    printf("Error %d\n",error=integral(fun,x(1),x(2))-resultSimp);
    printf("Relative error %d %%\n",errorRelative=abs(error/integral(fun,x(1),x(2)))*100);

endfunction

%/////////////////////////
%         //////////////////    simpsons 3/8 function declaration   ////////////////////////////
    
function simp3(x,y,fun)  

h=(x(2)-x(1))/3; %  h=B-A/3

%simp3/8 formula apo to wikipedia
resultSimp=((3*h)/8)*(y(1)+3*(fun((2*x(1)+x(2))/3)+fun((x(1)+2*x(2))/3))+y(2));

    printf("resultSimp3/8 %d\n",resultSimp);
    printf("True intergral result %d\n",integral(fun,x(1),x(2)));
    printf("Error %d\n",error=abs(integral(fun,x(1),x(2))-resultSimp));
    printf("Relative error %d %%\n",errorRelative=abs(error/integral(fun,x(1),x(2)))*100);

endfunction

%    ///////////////////    function call   ///////////////////

choice=input("Choose method: \n1.Trapezoid (composite included)\n2.Simpsons 1/3\n3.Simpsons 3/8 \n4.Skip\n");
switch choice
    case 1
        ftrap(x,y,fun)
        ftrapSteps(x,y,fun)
    case 2
       simp1(x,y,fun)
    case  3
        simp3(x,y,fun)
    otherwise
        printf("moving on\n");
endswitch

    %///////////////////////// INTERRPLOTATION SECTION ////////////////////////////
printf("\n INTERPLOTATION \n");
choice=input("Choose function to interplot \n1.cos(x)\n2.5*(x^5)+12*(x^4)-(x^3)+9*(x^2)\n3.2^x+c\n4.Skip\n");
c=2022202000082;

choice2=input("Gia logous aplothtas \n allazoume to am apo 2022202000082 se 82? yes=1 no=0\n");
    if choice2==1
        c=82;
    endif

switch choice
    case 1
        fun=@(x) cos(x) + c  ;
    case 2
        fun=@(x) 5*(x.^5)+12*(x.^4)-(x.^3)+9*(x.^2);
    case  3
        fun=@(x) 2.^x+c;   
    otherwise
        fun=@cos;
endswitch

printf("DATA POINTS x,f(x) \n");

   x=0:1:4   % TA 5 gnwsta shmeia
   xMexri10=0:1:9; %ta alla 5 shmeia
   yi=fun(x) %f(x) gia ta 5 shmeia
   

% ///////////////////////////// linear interplotation function ///////////////////////////////

function  linearInter(x,yi,fun)
        plotp=0; %init gia to plot
        plote=0; %init gia to plot tou error
        plotr=0; %init gia to plot tou Relative error
        for xx=0:9 
                    %ousiastika oti einai mesa se auto to for einai h pragmatikh parembolh
                    %0:9 giati trexei to function kai gia ta 5 gnwsta (loop 0 mexri 5) kai ta 5 agnwsta (loop 5 mexri 9)
            sum=0;
            trueResult=fun(xx);
            for i=1:2               %apo ta slides dinotan gia n-th order,ara gia grammikh 1st order (2-1=1)
                product=yi(i);
                for j=1:2
                    if i!=j
                        product=product*( ( xx-x(j) )/( x(i)-x(j)) );
                    endif
                endfor
                sum=sum+product;
            endfor
                printf("Linear interplotation \nx=%d f(x)=%d \nError=%d\nRelative error=%d %%\n",xx,sum,trueResult-sum,abs(((trueResult-sum)/trueResult)*100));
            plotp(end+1)=sum;
            plote(end+1)=trueResult-sum;
            plotr(end+1)=abs(((trueResult-sum)/trueResult)*100);
        endfor
    yy=0:1:9
    plotp(1) = [];  %delete thn init timh tou plotp
    plote(1) = [];  %delete thn init timh tou plote
    plotr(1) = [];  %delete thn init timh tou plote
    subplot(3,1,1)
    plot(yy,plotp,'rs-');
    title('Linear interplotation'),xlabel('x'),ylabel('y'),grid on
    subplot(3,1,2)
    plot(yy,plote,'rs-');
    title('Error of Linear interplotation'),xlabel('x'),ylabel('Error'),grid on
    subplot(3,1,3)
    plot(yy,plotr,'rs-');
    title('Relative error of Linear interplotation'),xlabel('x'),ylabel('Relative Error'),grid on
    plotp
    plote
    plotr
endfunction

% ////////////////////////// LA GRANGE function //////////////////////////////

% x to diastama ,yi=oi times gia to diasthma,n stoixeia sto diasthma, xx mia timh apo to diasthma
%fun exists to ensure the yi stays the same

function  laGrange(x,yi,n,fun)
        plotp=0; %init gia to plot
        plote=0; %init gia to plot tou error
        plotr=0; %init gia to plot tou Relative error
        for xx=0:9
                            %ousiastika oti einai mesa se auto to for einai h pragmatikh parembolh
                            %0:9 giati trexei to function kai gia ta 5 gnwsta (loop 0 mexri 5) kai ta 5 agnwsta (loop 5 mexri 9)
            sum=0;    
            trueResult=fun(xx);
            for i=1:n           %apo ta slides
                product=yi(i);
                for j=1:n
                    if i!=j
                        product=product*( ( xx-x(j) )/( x(i)-x(j)) );
                    endif
                endfor
                sum=sum+product;
            endfor
            printf("Lagrange interplotation\nx=%d f(x)=%d \nError=%d\nRelative error=%d %%\n",xx,sum,trueResult-sum,abs(((trueResult-sum)/trueResult)*100));
            plotp(end+1)=sum;
            plote(end+1)=trueResult-sum;
            plotr(end+1)=abs(((trueResult-sum)/trueResult)*100);
        endfor
    yy=0:1:9;
    plotp(1) = [];  %delete thn init timh tou plotp
    plote(1) = [];  %delete thn init timh tou plote
    plotr(1) = [];  %delete thn init timh tou plote
    subplot(3,1,1)
    plot(yy,plotp,'rs-');
    title('La grange interplotation'),xlabel('x'),ylabel('y'),grid on
    subplot(3,1,2)
    plot(yy,plote,'rs-');
    title('Error of La grange interplotation'),xlabel('x'),ylabel('Error'),grid on
    subplot(3,1,3)
    plot(yy,plotr,'rs-');
    title('Relative error of La grange interplotation'),xlabel('x'),ylabel('Relative Error'),grid on
    plotp
    plote
    plotr
endfunction
% //////////////////////////////// splines ////////////////////////////////////

%           nah

%///////////////////////////////// least squares ///////////////////////////////
            
function leastSquares(x,yi,n,xx)

sumGinomeno=0;sumPower2=0;sumX=0;sumY=0; 

    %calculating the sums
    for i=1:n
        sumGinomeno= sumGinomeno + x(i)*yi(i);
        sumPower2= sumPower2 + x(i)*x(i);
        sumY=sumY+yi(i);
        sumX=sumX+x(i);
    endfor
    %a and b for the y=a*x+b 
    a=  (n*sumGinomeno - sumX*sumY)/(n*sumPower2 - (sumX).^2);
    b= (sumPower2*sumY-sumX*sumGinomeno)/( n*sumPower2 - (sumX).^2);

    fun=@(x) a*x+b; %declaring our function via function handle

    plot(xx,fun(xx),"rs-");
    title('Least squares interplotation'),xlabel('x'),ylabel('y'),grid on
endfunction


% ///////////////////////////// function call ///////////////////////////////


choice=input("Select interplot  method \n1.Linear\n2.La grange\n3.Splines \n4.Least Squares\n5.Skip\n");
switch choice

    case 1
        linearInter(x,yi,fun);
    case 2
        laGrange(x,yi,5,fun); % n=5 giati kserw mono 5 shmeia
    case 3
        printf("oxi gia ta splines\n");
    case 4
        subplot(2,1,1)
            leastSquares(x,yi,5,x) % gia ta gnwsta shmeia
        subplot(2,1,2)
            leastSquares(x,yi,5,xMexri10) % gia ta agnwsta shmeia
    otherwise
        printf("moving on...\n");
endswitch

printf("\n**the end**\n");


%implement  compisite rules, ask people
%compisite trap is ok???


%Ε Σχολιάστε την αποτελεσματικότητα και την ακρίβεια των μεθόδων, την επίδραση της επιλογής της 
%μεθόδου στα σφάλματα και τυχόν διαφορές που παρατηρήθηκαν στην απόδοση των μεθόδων όταν 
%εφαρμόστηκαν στις διάφορες συναρτήσεις
%
% function handle for sin(x) , @(varUsed) function(px log(x))
    %plot(xx,fun(xx),"rs-");
    %title('Least squares interplotation'),xlabel('x'),ylabel('y'),grid on
            %subplot(2,1,1)
            %leastSquares(x,yi,5,x) % gia ta gnwsta shmeia
        %subplot(2,1,2)
            %leastSquares(x,yi,5,xMexri10) % gia ta agnwsta shmeia            