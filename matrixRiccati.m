function[dPdt] = matrixRiccati(t,A,B,Q,R,N,P,T)


P = reshape(P, size(A));

Q11 = reshape(Q(1,1,:),[length(T) 1]);
Q18 = reshape(Q(1,8,:),[length(T) 1]);
Q19 = reshape(Q(1,9,:),[length(T) 1]);
Q22 = reshape(Q(2,2,:),[length(T) 1]);
Q33 = reshape(Q(3,3,:),[length(T) 1]);
Q38 = reshape(Q(3,8,:),[length(T) 1]);
Q39 = reshape(Q(3,9,:),[length(T) 1]);
Q44 = reshape(Q(4,4,:),[length(T) 1]);
Q46 = reshape(Q(4,6,:),[length(T) 1]);
Q55 = reshape(Q(5,5,:),[length(T) 1]);
Q57 = reshape(Q(5,7,:),[length(T) 1]);
Q58 = reshape(Q(5,8,:),[length(T) 1]);
Q59 = reshape(Q(5,9,:),[length(T) 1]);
Q78 = reshape(Q(7,8,:),[length(T) 1]);
Q79 = reshape(Q(7,9,:),[length(T) 1]);
Q88 = reshape(Q(8,8,:),[length(T) 1]);
Q89 = reshape(Q(8,9,:),[length(T) 1]);
Q99 = reshape(Q(9,9,:),[length(T) 1]);

N21 = reshape(N(2,1,:),[length(T) 1]);
N32 = reshape(N(3,2,:),[length(T) 1]);
N43 = reshape(N(4,3,:),[length(T) 1]);
N54 = reshape(N(5,4,:),[length(T) 1]);
N82 = reshape(N(8,2,:),[length(T) 1]);
N84 = reshape(N(8,4,:),[length(T) 1]);
N92 = reshape(N(9,2,:),[length(T) 1]);
N94 = reshape(N(9,4,:),[length(T) 1]);




method = 'next';
Q11 = interp1(T,Q11,t,method); 
Q18 = interp1(T,Q18,t,method); 
Q19 = interp1(T,Q19,t,method); 
Q22 = interp1(T,Q22,t,method); 
Q33 = interp1(T,Q33,t,method); 
Q38 = interp1(T,Q38,t,method); 
Q39 = interp1(T,Q39,t,method); 
Q44 = interp1(T,Q44,t,method); 
Q46 = interp1(T,Q46,t,method); 
Q55 = interp1(T,Q55,t,method); 
Q57 = interp1(T,Q57,t,method); 
Q58 = interp1(T,Q58,t,method); 
Q59 = interp1(T,Q59,t,method); 
Q78 = interp1(T,Q78,t,method); 
Q79 = interp1(T,Q79,t,method); 
Q88 = interp1(T,Q88,t,method); 
Q89 = interp1(T,Q89,t,method); 
Q99 = interp1(T,Q99,t,method); 

N21 = interp1(T,N21,t,method);
N32 = interp1(T,N32,t,method);
N43 = interp1(T,N43,t,method);
N54 = interp1(T,N54,t,method);
N82 = interp1(T,N82,t,method);
N84 = interp1(T,N84,t,method);
N92 = interp1(T,N92,t,method);
N94 = interp1(T,N94,t,method);


Q = [Q11  0    0    0   0   0   0  Q18  Q19;
     0   Q22   0    0   0   0   0   0    0;
     0    0   Q33   0   0   0   0  Q38  Q39;
     0    0    0   Q44  0  Q46  0   0    0;
     0    0    0    0  Q55  0  Q57 Q58  Q59;
     0    0    0   Q46  0   1   0   0    0;
     0    0    0    0  Q57  0   1  Q78  Q79;
     Q18  0   Q38   0  Q58  0  Q78 Q88  Q89;
     Q19  0   Q39   0  Q59  0  Q79 Q89  Q99];
 

N = [0   0   0   0;
    N21  0   0   0;
     0  N32  0   0;
     0   0  N43  0;
     0   0   0  N54;
     0   0   1   0;
     0   0   0   1;
     0  N82  0  N84;
     0  N92  0  N94];


dPdt = -P*A -A.'*P +P*B*(R\B.')*P -Q +P*B*(R\N.') +N*(R\N.') +N*(R\B.')*P;
dPdt = dPdt(:);















% function[dPdt] = matrixRiccati(t,A,B,Q,R,N,P) %#ok<INUSL>
% P = reshape(P, size(A));
% dPdt = -P*A -A.'*P +P*B*(R\B.')*P -Q +P*B*(R\N.') +N*(R\N.') +N*(R\B.')*P;
% dPdt = dPdt(:);
% %fprintf('Value of Q is [%.2f ,%.2f,%.2f,%.2f,%.2f ,%.2f,%.2f,%.2f,%.2f]\n',Q)


end
