function A = curvspace(Ad,Au,n,m)

% construct a space for A with more grid ponit on the lower end                                                                                                                                          
% Ad = min of A                                                                                                                                                                                          
% Au = max of A                                                                                                                                                                                          
% n =  number of points                                                                                                                                                                                  
% m determines the curvature of A, increase m => A more curvy                                                                                                                                            
% May use m = 2                                                                                                                                                                                            
% output A                                                                                                                                                                                               


x = linspace(0,1,n);                                                                                                                                                             

A0 = x.^m;                                                                                                                                                                                         

A = A0*(Au-Ad)+Ad;

A(1)=Ad; A(end)=Au;       % This is added on Oct.13,2011
