C1=3.7418e8;
C2=1.4388e4;

ts=5780;
E=@(t,T) C1.*t.^(-5)./(exp(C2./t./T)-1);
E1=@(t,T) t.*(C1.*t.^(-5)./(exp(C2./t./T)-1));

% for i=1:1001
% T(i)=650+0.1*(i-1);
% LHS(i) = 0.9*integral(@(t) E(t,T(i)),eps,1.5)+0.1*integral(@(t) E(t,T(i)),1.5,inf); 
% end

up=0.83*integral(@(t) E(t,ts),eps,0.7)+1.0225*integral(@(t) E(t,ts),0.7,1.9)...
   -0.275*integral(@(t) E1(t,ts),0.7,1.9) + 0.5*integral(@(t) E(t,ts),1.9,2.8)+...
    1.816*integral(@(t) E(t,ts),2.8,3.5)-0.47*integral(@(t) E1(t,ts),2.8,3.5)+...
    [0.17*integral(@(t) E(t,ts),3.5,inf);

low=5.67e-8*ts^4;

alpha=up/low


