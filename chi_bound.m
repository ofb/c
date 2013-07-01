#!/usr/local/bin/MathematicaScript -script

# SetDirectory["."]
l=ReadList["catted"][[1]];
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

plotSmall=ListPlot[mag,ImageSize->Large,PlotLabel->"Maximum value sum takes for prime p, over all nontrivial chi and for \[Lambda]\[Element]{0,1,2,3}",PlotRange->{{0,250},{0,1500}}]

plotFull=ListPlot[mag,ImageSize->Large,PlotLabel->"Maximum value sum takes for prime p, over all nontrivial chi and for \[Lambda]\[Element]{0,1,2,3}"]

Export["matho/plotSmall.png",plotSmall]

Export["matho/plotFull.png",plotFull]



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
plots[m]=ListPlot[lambda[m],ImageSize->Large,PlotLabel->"Maximum value sum takes for prime p, over all nontrivial chi and for \[Lambda] index="<>ToString[m]];
Export["matho/plotLambda-"<>ToString[m]<>".png",plots[m]]
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

dPlot=ListPointPlot3D[d,ImageSize->Large,PlotRange->{{0,l[[length]][[1]][[1]]+1},All,All}, PlotLabel->"Character sum values as C-valued function over p, over all nontrivial chi and for \[Lambda]\[Element]{0,1,2,3}",ViewPoint->{.5,0,-1},PreserveImageOptions->False]

Export["matho/3DPlot.png",dPlot]



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
plots[m]=ListPointPlot3D[lambdaD[m],ImageSize->Large,PlotRange->{{0,l[[length]][[1]][[1]]+1},All,All}, PlotLabel->"Character sum values as C-valued function over p, over all nontrivial chi and for \[Lambda] index="<>ToString[m],ViewPoint->{.5,0,-1},PreserveImageOptions->False];
Export["matho/3DPlotLambda-"<>ToString[m]<>".png",plots[m]]
]
