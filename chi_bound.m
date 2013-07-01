#!/usr/local/bin/MathematicaScript -script

# SetDirectory["."]
l=ReadList["f20"][[1]];
length=Length[l]

mag=l;

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
Export["matho/plotLambda-"<>ToString[m]<>".dat",lambda[m]]
]



d=l;

length=Length[d]

For[n=1,n<=length,n++,
d[[n]][[2]] =Flatten[d[[n]][[2]]]
]

For[n=1,n<=length,n++,
d[[n]]=Outer[{#1,Chop[Re[#2]],Chop[Im[#2]]}&,d[[n]][[1]],d[[n]][[2]]][[1]]
]

d=Flatten[d,1];

Export["matho/3DPlot.dat",d]



lambdaLength=Length[l[[1]][[2]]]

For[n=1,n<=lambdaLength,n++,
lambdaD[n]=MapIndexed[If[#2[[2]]==1,#1,#1[[n]]]&,l,{2}]
]

For[n=1,n<=length,n++,
For[m=1,m<=lambdaLength,m++,
ltemp=lambdaD[m];
ltemp[[n]]=Outer[{#1,Chop[Re[#2]],Chop[Im[#2]]}&,ltemp[[n]][[1]],ltemp[[n]][[2]]][[1]];
lambdaD[m]=ltemp
]
]

For[m=1,m<=lambdaLength,m++,
ltemp=lambdaD[m];
ltemp=Flatten[ltemp,1];
lambdaD[m]=ltemp
]

For[m=1,m<=lambdaLength,m++,
Export["matho/lambda3D-"<>ToString[m]<>".png",lambdaD[m]]
]
