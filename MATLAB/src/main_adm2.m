% Main interface for admission control
%{
% Test Example 1:
Pr_ = [1/3 1/3 1/3; 1/3 1/3 1/3; 1/3 1/3 1/3];
rate_ = [2 1 0];
%groupA = [1 2];
%groupB = 3;
groupA = 3;
groupB = [1 2];
% End of Test Example 1
%}

%{
% Test Example 2:
Pr_ = [1/2 1/2 0; 0 1/2 1/2; 1/2 0 1/2];
rate_ = [2 1 0];
%groupA = [1 2];
%groupB = 3;
groupA = 1;
groupB = [2 3];
% End of Test Example 2
%}

% Test Example 3:
% 
%{
rate_ = [9 8 6 4 3 2 1.5 1 0];
row1_ = [1/9 1/9 1/9 1/9 1/9 1/9 1/9 1/9 1/9];
row2_ = [1/9 1/9 1/9 1/9 1/9 1/9 1/9 1/9 1/9];
row3_ = [1/9 1/9 1/9 1/9 1/9 1/9 1/9 1/9 1/9];
row4_ = [1/9 1/9 1/9 1/9 1/9 1/9 1/9 1/9 1/9];
row5_ = [1/9 1/9 1/9 1/9 1/9 1/9 1/9 1/9 1/9];
row6_ = [1/9 1/9 1/9 1/9 1/9 1/9 1/9 1/9 1/9];
row7_ = [1/9 1/9 1/9 1/9 1/9 1/9 1/9 1/9 1/9];
row8_ = [1/9 1/9 1/9 1/9 1/9 1/9 1/9 1/9 1/9];
row9_ = [1/9 1/9 1/9 1/9 1/9 1/9 1/9 1/9 1/9];
row10_ = [1/9 1/9 1/9 1/9 1/9 1/9 1/9 1/9 1/9];
row11_ = [1/9 1/9 1/9 1/9 1/9 1/9 1/9 1/9 1/9];
row12_ = [1/9 1/9 1/9 1/9 1/9 1/9 1/9 1/9 1/9];
row13_ = [1/9 1/9 1/9 1/9 1/9 1/9 1/9 1/9 1/9];
row14_ = [1/9 1/9 1/9 1/9 1/9 1/9 1/9 1/9 1/9];
row15_ = [1/9 1/9 1/9 1/9 1/9 1/9 1/9 1/9 1/9];
row16_ = [1/9 1/9 1/9 1/9 1/9 1/9 1/9 1/9 1/9];
row17_ = [1/9 1/9 1/9 1/9 1/9 1/9 1/9 1/9 1/9];
row18_ = [1/9 1/9 1/9 1/9 1/9 1/9 1/9 1/9 1/9];
row19_ = [1/9 1/9 1/9 1/9 1/9 1/9 1/9 1/9 1/9];
row20_ = [1/9 1/9 1/9 1/9 1/9 1/9 1/9 1/9 1/9];
%}
%{


row1_ = [0 0.2 0.2 0.2 0.2 0.2];
row2_ = [0 0.2 0.2 0.2 0.2 0.2];
row3_ = [0 0.2 0.2 0.2 0.2 0.2];
row4_ = [0 0.2 0.2 0.2 0.2 0.2];
row5_ = [0 0.2 0.2 0.2 0.2 0.2];
row6_ = [0 0.2 0.2 0.2 0.2 0.2];
row7_ = [0 0.2 0.2 0.2 0.2 0.2];
row8_ = [0 0.2 0.2 0.2 0.2 0.2];
row9_ = [0 0.2 0.2 0.2 0.2 0.2];
row10_ = [0 0.2 0.2 0.2 0.2 0.2];

row1_ = [1/30 1/6 1/6 1/6 1/6 9/30];
row2_ = [1/30 1/6 1/6 1/6 1/6 9/30];
row3_ = [1/30 1/6 1/6 1/6 1/6 9/30];
row4_ = [1/30 1/6 1/6 1/6 1/6 9/30];
row5_ = [1/30 1/6 1/6 1/6 1/6 9/30];
row6_ = [1/30 1/6 1/6 1/6 1/6 9/30];
row7_ = [1/30 1/6 1/6 1/6 1/6 9/30];
row8_ = [1/30 1/6 1/6 1/6 1/6 9/30];
row9_ = [1/30 1/6 1/6 1/6 1/6 9/30];
row10_ = [1/30 1/6 1/6 1/6 1/6 9/30];

row1_ = [1/60 1/6 1/6 1/6 1/6 19/60];
row2_ = [1/60 1/6 1/6 1/6 1/6 19/60];
row3_ = [1/60 1/6 1/6 1/6 1/6 19/60];
row4_ = [1/60 1/6 1/6 1/6 1/6 19/60];
row5_ = [1/60 1/6 1/6 1/6 1/6 19/60];
row6_ = [1/60 1/6 1/6 1/6 1/6 19/60];
row7_ = [1/60 1/6 1/6 1/6 1/6 19/60];
row8_ = [1/60 1/6 1/6 1/6 1/6 19/60];
row9_ = [1/60 1/6 1/6 1/6 1/6 19/60];
row10_ = [1/60 1/6 1/6 1/6 1/6 19/60];

row1_ = [1/36 1/6 1/6 1/6 1/6 11/36];
row2_ = [1/36 1/6 1/6 1/6 1/6 11/36];
row3_ = [1/36 1/6 1/6 1/6 1/6 11/36];
row4_ = [1/36 1/6 1/6 1/6 1/6 11/36];
row5_ = [1/36 1/6 1/6 1/6 1/6 11/36];
row6_ = [1/36 1/6 1/6 1/6 1/6 11/36];
row7_ = [1/36 1/6 1/6 1/6 1/6 11/36];
row8_ = [1/36 1/6 1/6 1/6 1/6 11/36];
row9_ = [1/36 1/6 1/6 1/6 1/6 11/36];
row10_ = [1/36 1/6 1/6 1/6 1/6 11/36];


row1_ = [1/30 1/6 1/6 1/6 1/6 9/30];
row2_ = [1/30 1/6 1/6 1/6 1/6 9/30];
row3_ = [1/30 1/6 1/6 1/6 1/6 9/30];
row4_ = [1/30 1/6 1/6 1/6 1/6 9/30];
row5_ = [1/30 1/6 1/6 1/6 1/6 9/30];
row6_ = [1/30 1/6 1/6 1/6 1/6 9/30];
row7_ = [1/30 1/6 1/6 1/6 1/6 9/30];
row8_ = [1/30 1/6 1/6 1/6 1/6 9/30];
row9_ = [1/30 1/6 1/6 1/6 1/6 9/30];
row10_ = [1/30 1/6 1/6 1/6 1/6 9/30];

row1_ = [1/6 1/6 1/6 1/6 1/6 1/6];
row2_ = [1/6 1/6 1/6 1/6 1/6 1/6];
row3_ = [1/6 1/6 1/6 1/6 1/6 1/6];
row4_ = [1/6 1/6 1/6 1/6 1/6 1/6];
row5_ = [1/6 1/6 1/6 1/6 1/6 1/6];
row6_ = [1/6 1/6 1/6 1/6 1/6 1/6];
row7_ = [1/6 1/6 1/6 1/6 1/6 1/6];
row8_ = [1/6 1/6 1/6 1/6 1/6 1/6];
row9_ = [1/6 1/6 1/6 1/6 1/6 1/6];
row10_ = [1/6 1/6 1/6 1/6 1/6 1/6];

row1_ = [1/25 6/25 6/25 6/25 6/25];
row2_ = [1/25 6/25 6/25 6/25 6/25];
row3_ = [1/25 6/25 6/25 6/25 6/25];
row4_ = [1/25 6/25 6/25 6/25 6/25];
row5_ = [1/25 6/25 6/25 6/25 6/25];
row6_ = [1/25 6/25 6/25 6/25 6/25];
row7_ = [1/25 6/25 6/25 6/25 6/25];
row8_ = [1/25 6/25 6/25 6/25 6/25];
row9_ = [1/25 6/25 6/25 6/25 6/25];
row10_ = [1/25 6/25 6/25 6/25 6/25];

row1_ = [1/10 1/10 1/10 1/10 6/10];
row2_ = [1/10 1/10 1/10 1/10 6/10];
row3_ = [1/10 1/10 1/10 1/10 6/10];
row4_ = [1/10 1/10 1/10 1/10 6/10];
row5_ = [1/10 1/10 1/10 1/10 6/10];
row6_ = [1/10 1/10 1/10 1/10 6/10];
row7_ = [1/10 1/10 1/10 1/10 6/10];
row8_ = [1/10 1/10 1/10 1/10 6/10];
row9_ = [1/10 1/10 1/10 1/10 6/10];
row10_ = [1/10 1/10 1/10 1/10 6/10];
row11_ = [0.5 0.5 0 0 0];
row12_ = [0.5 0.5 0 0 0];
row13_ = [0.5 0.5 0 0 0];
row14_ = [0.5 0.5 0 0 0];
row15_ = [0.5 0.5 0 0 0];
row16_ = [0.5 0.5 0 0 0];
row17_ = [0.5 0.5 0 0 0];
row18_ = [0.5 0.5 0 0 0];
row19_ = [0.5 0.5 0 0 0];
row20_ = [0.5 0.5 0 0 0];

row1_ = [0.02 0.245 0.245 0.245 0.245];
row2_ = [0.02 0.245 0.245 0.245 0.245];
row3_ = [0.02 0.245 0.245 0.245 0.245];
row4_ = [0.02 0.245 0.245 0.245 0.245];
row5_ = [0.02 0.245 0.245 0.245 0.245];
row6_ = [0.02 0.245 0.245 0.245 0.245];
row7_ = [0.02 0.245 0.245 0.245 0.245];
row8_ = [0.02 0.245 0.245 0.245 0.245];
row9_ = [0.02 0.245 0.245 0.245 0.245];
row10_ = [0.02 0.245 0.245 0.245 0.245];
row11_ = [0.2 0.2 0.2 0.2 0.2];
row12_ = [0.2 0.2 0.2 0.2 0.2];
row13_ = [0.2 0.2 0.2 0.2 0.2];
row14_ = [0.2 0.2 0.2 0.2 0.2];
row15_ = [0.2 0.2 0.2 0.2 0.2];
row16_ = [0.2 0.2 0.2 0.2 0.2];
row17_ = [0.2 0.2 0.2 0.2 0.2];
row18_ = [0.2 0.2 0.2 0.2 0.2];
row19_ = [0.2 0.2 0.2 0.2 0.2];
row20_ = [0.2 0.2 0.2 0.2 0.2];

%}
%rate_ = [9 8 6 3 1 0];
%rate_ = [2340 2080 1560 750 220];
rate_ = [2340 0];
row1_ = [0.5 0.5];
row2_ = [0.5 0.5];
row3_ = [0.5 0.5];
row4_ = [0.5 0.5];
row5_ = [0.5 0.5];
row6_ = [0.5 0.5];
row7_ = [0.5 0.5];
row8_ = [0.5 0.5];
row9_ = [0.5 0.5];
row10_ = [0.5 0.5];
row11_ = [0.5 0.5];
row12_ = [0.5 0.5];
row13_ = [0.5 0.5];
row14_ = [0.5 0.5];
row15_ = [0.5 0.5];
row16_ = [0.5 0.5];
row17_ = [0.5 0.5];
row18_ = [0.5 0.5];
row19_ = [0.5 0.5];
row20_ = [0.5 0.5];
Pr_ = [row1_; row2_; row3_; row4_; row5_; row6_; row7_; row8_; row9_; row10_;row11_; row12_; row13_; row14_; row15_; row16_; row17_; row18_; row19_; row20_];
%groupA = [1 2];
%groupB = 3;
%groupA = [1 2 3 4 5 6 7 8 9];
%groupB = 10;
groupA = [1 2 3 4 5 6 7 8 9 10];
groupB = [11 12 13 14 15 16 17 18 19 20];
% End of Test Example 3


[q1, q2, q3, hist_, var_] = admctrl2(Pr_, rate_, groupA, groupB);



