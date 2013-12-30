#!/usr/local/bin/MathematicaScript -script

# SetDirectory["."]
old=ReadList["tmpold"][[1]];
new=ReadList["tmpnew"][[1]];
length=Length[old];

c1=Table[Insert[old,new[[i]][[2]][[1]],{i,2,3}],{i,1,length}][[1]];
c1=Table[Insert[c1,new[[i]][[2]][[2]],{i,2,4}],{i,1,length}][[1]];

mag=c1;

For[n=1,n<=length,n++,
mag[[n]][[2]] =Flatten[mag[[n]][[2]]]
]

For[n=1,n<=length,n++,
mag[[n]][[2]] = Map[Abs,mag[[n]][[2]]]
]

For[n=1,n<=length,n++,
mag[[n]][[2]] = Max[mag[[n]][[2]]]
]

For[n=1,n<=length,n++,
mag[[n]][[1]] = mag[[n]][[1]] [[1]]
]

Export["matho/absoluteMagnitude.dat",mag]



lambdaLength=Length[l[[1]][[2]]];

For[n=1,n<=lambdaLength,n++,
lambda[n]=MapIndexed[If[#2[[2]]==1,#1,#1[[n]]]&,l,{2}]
]

For[n=1,n<=length,n++,
For[m=1,m<=lambdaLength,m++,
ltemp=lambda[m];
ltemp[[n]][[2]]=Map[Abs,ltemp[[n]][[2]]];
lambda[m]=ltemp
]
]

For[n=1,n<=length,n++,
For[m=1,m<=lambdaLength,m++,
ltemp=lambda[m];
ltemp[[n]][[2]] = Max[ltemp[[n]][[2]]];
lambda[m]=ltemp;
]
]

For[n=1,n<=length,n++,
For[m=1,m<=lambdaLength,m++,
ltemp=lambda[m];
ltemp[[n]][[1]] = ltemp[[n]][[1]] [[1]];
lambda[m]=ltemp;
]
]

For[m=1,m<=lambdaLength,m++,
Export["matho/lambda-"<>ToString[m-1]<>".dat",lambda[m]]
]
