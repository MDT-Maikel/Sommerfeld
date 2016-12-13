(* ::Package:: *)

(* ::Title:: *)
(*Plot Utility*)


(* ::Text:: *)
(*Mathematica package for plotting and exporting plots.*)


(* ::Section:: *)
(*Begin Package*)


(* ::Subsection:: *)
(*Mathematica Version*)


(* ::Input::Initialization:: *)
(* Check for the correct Mathematica version *)
PU$MathematicaVersionNeeded=10;
If[Floor[$VersionNumber]<PU$MathematicaVersionNeeded,
Print[Style["PlotUtility requires at least Mathematica version "<>ToString[PU$MathematicaVersionNeeded]<>" while this is only version "<>ToString[Floor[$VersionNumber]]<>".",Red,Bold]];
Abort[];
];


(* ::Subsection:: *)
(*Begin Package*)


(* ::Input::Initialization:: *)
(* Begin the package and clear the PlotUtility` namespace *)
BeginPackage["PlotUtility`"];
Unprotect@@Names["PlotUtility`*"];
ClearAll@@Names["PlotUtility`*"];


(* ::Subsection:: *)
(*Package Version*)


(* ::Input::Initialization:: *)
(* Usage Descriptions *)
PU$Version::usage="PlotUtility version";

(* Begin Definitions *)
Begin["`Private`"];

PU$Version="0.2";

(* End Definitions *)
End[];

(* Protect & Distribute Definitions *)
DistributeDefinitions[PU$Version];
Protect[PU$Version];


(* ::Section:: *)
(*Utility Functions*)


(* ::Input::Initialization:: *)
(* Usage Descriptions *)
StrFrac::usage="StrFrac[a, b]: Converts the fraction a/b to a string.";
CleanContourPlot::usage="CleanContourPlot[cp]: Cleans up the contour plot cp for exporting.";
EdgeRule::usage="EdgeRule: ";
RestylePlot::usage="RestylePlot: ";

(* Begin Definitions *)
Begin["`Private`"];

StrFrac[a_,b_]:=ToString[FractionBox[a,b,Beveled->True]//DisplayForm,FormatType->StandardForm];

(* http://mathematica.stackexchange.com/questions/3190/saner-alternative-to-contourplot-fill *)
CleanContourPlot[cp_]:=Module[{points,groups,regions,lines},
groups=Cases[cp,{style__,g_GraphicsGroup}:>{{style},g},Infinity];
points=First@Cases[cp,GraphicsComplex[pts_,___]:>pts,Infinity];
regions=Table[
Module[{group,style,polys,edges,cover,graph},
{style,group}=g;
polys=Join@@Cases[group,Polygon[pt_,___]:>pt,Infinity];
edges=Join@@(Partition[#,2,1,1]&/@polys);
cover=Cases[Tally[Sort/@edges],{e_,1}:>e];
graph=Graph[UndirectedEdge@@@cover];
{Sequence@@style,FilledCurve[List/@Line/@First/@Map[First,FindEulerianCycle/@(Subgraph[graph,#]&)/@ConnectedComponents[graph],{3}]]}],
{g,groups}
];
lines=Cases[cp,_Tooltip,Infinity];
Graphics[GraphicsComplex[points,{regions,lines}],Sequence@@Options[cp]]
];

EdgeRule:={EdgeForm[],r_?(MemberQ[{RGBColor,Hue,CMYKColor,GrayLevel},Head[#]]&),i___}:>{EdgeForm[r],r,i};

(* http://mathematica.stackexchange.com/questions/17250/is-it-possible-to-change-the-color-of-plot-in-show *)
RestylePlot[plot_Graphics,styles_List,fillings_List,op:OptionsPattern[Graphics]]:=Module[{x=styles,y=fillings},
Show[MapAt[#/.{{__,ln__Line}:>{Directive@Last[x=RotateLeft@x],ln},{EdgeForm[],__,GraphicsGroup[pts__]}:>{EdgeForm[],Directive@Last[y=RotateLeft@y],GraphicsGroup[pts]}}&,plot,1],op]
];

(* End Definitions *)
End[];

(* Protect & Distribute Definitions *)
DistributeDefinitions[StrFrac,CleanContourPlot,EdgeRule,RestylePlot];
Protect[StrFrac,CleanContourPlot,EdgeRule,RestylePlot];


(* ::Section:: *)
(*Plot Options*)


(* ::Input::Initialization:: *)
(* Some default options for plots *)
Unprotect[PlotOptions];
ClearAll[PlotOptions];
PlotOptions=Sequence[
ImageSize->800,
Axes->False,
Frame->True,
FrameStyle->Directive[20,Thickness[0.002],Background->White],
LabelStyle->Directive[Background->White],
GridLinesStyle->Directive[Dashed,Opacity[0.75],Thickness[0.0015],Gray]
];
(* Protect all variables *)
Protect[PlotOptions];
(* Distribute definition PlotOptions to all kernels *)
DistributeDefinitions[PlotOptions];
(* Set the new plot options *)
SetOptions[Plot,PlotOptions];
SetOptions[LogPlot,PlotOptions];
SetOptions[RegionPlot,PlotOptions];
SetOptions[ContourPlot,PlotOptions];
SetOptions[ListPlot,PlotOptions];
SetOptions[ListDensityPlot,PlotOptions];
SetOptions[ListLogPlot,PlotOptions];
SetOptions[ParametricPlot,PlotOptions];


(* ::Section:: *)
(*Plot Styles*)


(* ::Subsection:: *)
(*Plot Colors, Line Styles & Region Styles*)


(* ::Input::Initialization:: *)
(* Usage Descriptions *)
PlotGridLinesStyle::usage="PlotGridLinesStyle: A predefined grid line style.";
PlotFrameStyleSmall::usage="PlotFrameStyleSmall: A predefined plot frame style (small, for articles or large images).";
PlotFrameStyleNormal::usage="PlotFrameStyleNormal: A predefined plot frame style (normal, for articles or large images).";
PlotFrameStyleLarge::usage="PlotFrameStyleLarge: A predefined plot frame style (large, for presentations or normal size images).";
PlotFrameStyleHuge::usage="PlotFrameStyleHuge: A predefined plot frame style (huge, for presentations, or small images).";

FontSizeSmall::usage="FontSizeSmall: 22";
FontSizeNormal::usage="FontSizeNormal: 24";
FontSizeLarge::usage="FontSizeLarge: 26";
FontSizeHuge::usage="FontSizeHuge: 30";

PlotColor1::usage="PlotColor1: A predefined plot color ("<>ToString[PlotColor1,StandardForm]<>").";
PlotColor2::usage="PlotColor2: A predefined plot color ("<>ToString[PlotColor2,StandardForm]<>").";
PlotColor3::usage="PlotColor3: A predefined plot color ("<>ToString[PlotColor3,StandardForm]<>").";
PlotColor4::usage="PlotColor4: A predefined plot color ("<>ToString[PlotColor4,StandardForm]<>").";
PlotColor5::usage="PlotColor5: A predefined plot color ("<>ToString[PlotColor5,StandardForm]<>").";
PlotColor6::usage="PlotColor6: A predefined plot color ("<>ToString[PlotColor6,StandardForm]<>").";

PlotColorAlt1::usage="PlotColor1: A predefined alternative plot color ("<>ToString[PlotColorAlt1,StandardForm]<>").";
PlotColorAlt2::usage="PlotColor2: A predefined alternative plot color ("<>ToString[PlotColorAlt2,StandardForm]<>").";
PlotColorAlt3::usage="PlotColor3: A predefined alternative plot color ("<>ToString[PlotColorAlt3,StandardForm]<>").";
PlotColorAlt4::usage="PlotColor4: A predefined alternative plot color ("<>ToString[PlotColorAlt4,StandardForm]<>").";
PlotColorAlt5::usage="PlotColor5: A predefined alternative plot color ("<>ToString[PlotColorAlt5,StandardForm]<>").";
PlotColorAlt6::usage="PlotColor6: A predefined alternative plot color ("<>ToString[PlotColorAlt6,StandardForm]<>").";

LineStyle1::usage="LineStyle1: A predefined plot line style.";
LineStyle2::usage="LineStyle2: A predefined plot line style.";
LineStyle3::usage="LineStyle3: A predefined plot line style.";
LineStyle4::usage="LineStyle4: A predefined plot line style.";
LineStyle5::usage="LineStyle5: A predefined plot line style.";
LineStyle6::usage="LineStyle6: A predefined plot line style.";

RegionStyle1::usage="RegionStyle1: A predefined plot region style.";
RegionStyle2::usage="RegionStyle2: A predefined plot region style.";
RegionStyle3::usage="RegionStyle3: A predefined plot region style.";
RegionStyle4::usage="RegionStyle4: A predefined plot region style.";
RegionStyle5::usage="RegionStyle5: A predefined plot region style.";
RegionStyle6::usage="RegionStyle6: A predefined plot region style.";

(* Begin Definitions *)
Begin["`Private`"];

PlotGridLinesStyle=Directive[Dashed, Opacity[0.75],Thickness[0.0015], Gray];
PlotFrameStyleSmall=Directive[22,Thickness[0.002]];
PlotFrameStyleNormal=Directive[24,Thickness[0.00225]];
PlotFrameStyleLarge=Directive[26,Thickness[0.0025]];
PlotFrameStyleHuge=Directive[30,Thickness[0.003]];

FontSizeSmall=22;
FontSizeNormal=24;
FontSizeLarge=26;
FontSizeHuge=30;

PlotColor1=Blend[{ColorData[97,"ColorList"][[1]],Black},0.1];
PlotColor2=Blend[{ColorData[97,"ColorList"][[2]],Black},0.1];
PlotColor3=Blend[{ColorData[97,"ColorList"][[3]],ColorData[97,"ColorList"][[6]],Black},{0.81,0.09,0.10}];
PlotColor4=Blend[{ColorData[97,"ColorList"][[4]],Black},0.05];
PlotColor5=Blend[{ColorData[97,"ColorList"][[5]],Black},0.1];
PlotColor6=Blend[{ColorData[97,"ColorList"][[6]],Black},0.1];

PlotColorAlt1=Blend[{PlotColor1,Black},0.15];
PlotColorAlt2=Blend[{PlotColor2,Black},0.15];
PlotColorAlt3=Blend[{PlotColor3,Black},0.15];
PlotColorAlt4=Blend[{PlotColor4,Black},0.15];
PlotColorAlt5=Blend[{PlotColor5,Black},0.15];
PlotColorAlt6=Blend[{PlotColor6,Black},0.15];

LineStyle1=Directive[Thick,PlotColor1];
LineStyle2=Directive[Thick,PlotColor2];
LineStyle3=Directive[Thick,PlotColor3];
LineStyle4=Directive[Thick,PlotColor4];
LineStyle5=Directive[Thick,PlotColor5];
LineStyle6=Directive[Thick,PlotColor6];

RegionStyle1=Directive[PlotColor1,Opacity[1/3]];
RegionStyle2=Directive[PlotColor2,Opacity[1/3]];
RegionStyle3=Directive[PlotColor3,Opacity[1/3]];
RegionStyle4=Directive[PlotColor4,Opacity[1/3]];
RegionStyle5=Directive[PlotColor5,Opacity[1/3]];
RegionStyle6=Directive[PlotColor6,Opacity[1/3]];

(* End Definitions *)
End[];

(* Protect & Distribute Definitions *)
DistributeDefinitions[PlotGridLinesStyle,PlotFrameStyleSmall,PlotFrameStyleNormal,PlotFrameStyleLarge,PlotFrameStyleHuge,FontSizeSmall,FontSizeNormal,FontSizeLargeFontSizeHuge,PlotColor1,PlotColor2,PlotColor3,PlotColor4,PlotColor5,PlotColor6,PlotColorAlt1,PlotColorAlt2,PlotColorAlt3,PlotColorAlt4,PlotColorAlt5,PlotColorAlt6,LineStyle1,LineStyle2,LineStyle3,LineStyle4,LineStyle5,LineStyle6,RegionStyle1,RegionStyle2,RegionStyle3,RegionStyle4,RegionStyle5,RegionStyle6];
Protect[PlotGridLinesStyle,PlotFrameStyleSmall,PlotFrameStyleNormal,PlotFrameStyleLarge,PlotFrameStyleHuge,FontSizeSmall,FontSizeNormal,FontSizeLargeFontSizeHuge,PlotColor1,PlotColor2,PlotColor3,PlotColor4,PlotColor5,PlotColor6,PlotColorAlt1,PlotColorAlt2,PlotColorAlt3,PlotColorAlt4,PlotColorAlt5,PlotColorAlt6,LineStyle1,LineStyle2,LineStyle3,LineStyle4,LineStyle5,LineStyle6,RegionStyle1,RegionStyle2,RegionStyle3,RegionStyle4,RegionStyle5,RegionStyle6];


(* ::Subsection::Closed:: *)
(*Frame Styles & Grid Lines [OLD]*)


(* ::Input:: *)
(*(* List of nice frame ticks for the coupling constant g *)*)
(*Unprotect[FrameTicksTenth,FrameTicksHalf,FrameTicksOne,FrameTicksPi8,FrameTicksPi4,FrameTicksPi2];*)
(*ClearAll[FrameTicksTenth,FrameTicksHalf,FrameTicksOne,FrameTicksPi8,FrameTicksPi4,FrameTicksPi2];*)
(*FrameTicksTenth={{{0.0,"0.0",{0.01,0}},{0.05,"",{0.005,0}},{0.1,"0.1",{0.01,0}},{0.15,"",{0.005,0}},{0.2,"0.2",{0.01,0}},{0.25,"",{0.005,0}},{0.3,"0.3",{0.01,0}},{0.35,"",{0.005,0}},{0.4,"0.4",{0.01,0}},{0.45,"",{0.005,0}},{0.5,"0.5",{0.01,0}},{0.55,"",{0.005,0}},{0.6,"0.6",{0.01,0}},{0.65,"",{0.005,0}},{0.7,"0.7",{0.01,0}},{0.75,"",{0.005,0}},{0.8,"0.8",{0.01,0}},{0.85,"",{0.005,0}},{0.9,"0.9",{0.01,0}},{0.95,"",{0.005,0}},{1.0,"1.0",{0.01,0}}},None};*)
(*FrameTicksHalf={{{0.0,"0.0",{0.01,0}},{0.25,"",{0.005,0}},{0.5,"0.5",{0.01,0}},{0.75,"",{0.005,0}},{1.0,"1.0",{0.01,0}},{1.25,"",{0.005,0}},{1.5,"1.5",{0.01,0}},{1.75,"",{0.005,0}},{2.0,"2.0",{0.01,0}},{2.25,"",{0.005,0}},{2.5,"2.5",{0.01,0}},{2.75,"",{0.005,0}},{3.0,"3.0",{0.01,0}},{3.25,"",{0.005,0}},{3.5,"3.5",{0.01,0}},{3.75,"",{0.005,0}},{4.0,"4.0",{0.01,0}},{4.25,"",{0.005,0}},{4.5,"4.5",{0.01,0}},{4.75,"",{0.005,0}},{5.0,"5.0",{0.01,0}}},None};*)
(*FrameTicksOne={{{0,"0",{0.01,0}},{0.5,"",{0.005,0}},{1,"1",{0.01,0}},{1.5,"",{0.005,0}},{2,"2",{0.01,0}},{2.5,"",{0.005,0}},{3,"3",{0.01,0}},{3.5,"",{0.005,0}},{4,"4",{0.01,0}},{4.5,"",{0.005,0}},{5,"5",{0.01,0}},{5.5,"",{0.005,0}},{6,"6",{0.01,0}},{6.5,"",{0.005,0}},{7,"7",{0.01,0}},{7.5,"",{0.005,0}},{8,"8",{0.01,0}},{8.5,"",{0.005,0}},{9,"9",{0.01,0}},{9.5,"",{0.005,0}},{10,"10",{0.01,0}}},None};*)
(*FrameTicksPi8={{{0.0,"0",{0.01,0}},{\[Pi]/16,"",{0.005,0}},{\[Pi]/8,"\!\(\*FractionBox[\(\[Pi]\), \(8\)]\)",{0.01,0}},{3\[Pi]/16,"",{0.005,0}},{\[Pi]/4,"\!\(\*FractionBox[\(\[Pi]\), \(4\)]\)",{0.01,0}},{5\[Pi]/16,"",{0.005,0}},{3\[Pi]/8,"\!\(\*FractionBox[\(3  \[Pi]\), \(8\)]\)",{0.01,0}},{7\[Pi]/16,"",{0.005,0}},{\[Pi]/2,"\!\(\*FractionBox[\(\[Pi]\), \(2\)]\)",{0.01,0}},{9\[Pi]/16,"",{0.005,0}},{5\[Pi]/8,"\!\(\*FractionBox[\(5  \[Pi]\), \(8\)]\)",{0.01,0}},{11\[Pi]/16,"",{0.005,0}},{3\[Pi]/4,"\!\(\*FractionBox[\(3  \[Pi]\), \(4\)]\)",{0.01,0}},{13\[Pi]/16,"",{0.005,0}},{7\[Pi]/8,"\!\(\*FractionBox[\(7  \[Pi]\), \(8\)]\)",{0.01,0}},{15\[Pi]/16,"",{0.005,0}},{\[Pi],"\[Pi]",{0.01,0}},{17\[Pi]/16,"",{0.005,0}},{9\[Pi]/8,"\!\(\*FractionBox[\(9  \[Pi]\), \(8\)]\)",{0.01,0}},{19\[Pi]/16,"",{0.005,0}},{5\[Pi]/4,"\!\(\*FractionBox[\(5  \[Pi]\), \(4\)]\)",{0.01,0}},{21\[Pi]/16,"",{0.005,0}},{11\[Pi]/8,"\!\(\*FractionBox[\(11  \[Pi]\), \(8\)]\)",{0.01,0}},{23\[Pi]/16,"",{0.005,0}},{3\[Pi]/2,"\!\(\*FractionBox[\(3  \[Pi]\), \(2\)]\)",{0.01,0}},{25\[Pi]/16,"",{0.005,0}},{13\[Pi]/8,"\!\(\*FractionBox[\(13  \[Pi]\), \(8\)]\)",{0.01,0}},{27\[Pi]/16,"",{0.005,0}},{7\[Pi]/4,"\!\(\*FractionBox[\(7  \[Pi]\), \(4\)]\)",{0.01,0}},{29\[Pi]/16,"",{0.005,0}},{15\[Pi]/8,"\!\(\*FractionBox[\(15  \[Pi]\), \(8\)]\)",{0.01,0}},{31\[Pi]/16,"",{0.005,0}},{2\[Pi],"2\[Pi]",{0.01,0}}},None};*)
(*FrameTicksPi4={{{0.0,"0",{0.01,0}},{\[Pi]/8,"",{0.005,0}},{\[Pi]/4,"\!\(\*FractionBox[\(\[Pi]\), \(4\)]\)",{0.01,0}},{3\[Pi]/8,"",{0.005,0}},{\[Pi]/2,"\!\(\*FractionBox[\(\[Pi]\), \(2\)]\)",{0.01,0}},{5\[Pi]/8,"",{0.005,0}},{3\[Pi]/4,"\!\(\*FractionBox[\(3  \[Pi]\), \(4\)]\)",{0.01,0}},{7\[Pi]/8,"",{0.005,0}},{\[Pi],"\[Pi]",{0.01,0}},{9\[Pi]/8,"",{0.005,0}},{5\[Pi]/4,"\!\(\*FractionBox[\(5  \[Pi]\), \(4\)]\)",{0.01,0}},{11\[Pi]/8,"",{0.005,0}},{3\[Pi]/2,"\!\(\*FractionBox[\(3  \[Pi]\), \(2\)]\)",{0.01,0}},{13\[Pi]/8,"",{0.005,0}},{7\[Pi]/4,"\!\(\*FractionBox[\(7  \[Pi]\), \(4\)]\)",{0.01,0}},{15\[Pi]/8,"",{0.005,0}},{2\[Pi],"2\[Pi]",{0.01,0}},{17\[Pi]/8,"",{0.005,0}},{9\[Pi]/4,"\!\(\*FractionBox[\(9  \[Pi]\), \(4\)]\)",{0.01,0}},{19\[Pi]/8,"",{0.005,0}},{5\[Pi]/2,"\!\(\*FractionBox[\(5  \[Pi]\), \(2\)]\)",{0.01,0}},{21\[Pi]/8,"",{0.005,0}},{11\[Pi]/4,"\!\(\*FractionBox[\(11  \[Pi]\), \(4\)]\)",{0.01,0}},{23\[Pi]/8,"",{0.005,0}},{3\[Pi],"3\[Pi]",{0.01,0}},{25\[Pi]/8,"",{0.005,0}},{13\[Pi]/4,"\!\(\*FractionBox[\(13  \[Pi]\), \(4\)]\)",{0.01,0}},{27\[Pi]/8,"",{0.005,0}},{7\[Pi]/2,"\!\(\*FractionBox[\(7  \[Pi]\), \(2\)]\)",{0.01,0}},{29\[Pi]/8,"",{0.005,0}},{15\[Pi]/4,"\!\(\*FractionBox[\(15  \[Pi]\), \(4\)]\)",{0.01,0}},{31\[Pi]/8,"",{0.005,0}},{4\[Pi],"4\[Pi]",{0.01,0}}},None};*)
(*FrameTicksPi2={{{0.0,"0",{0.01,0}},{\[Pi]/4,"",{0.005,0}},{\[Pi]/2,"\!\(\*FractionBox[\(\[Pi]\), \(2\)]\)",{0.01,0}},{3\[Pi]/4,"",{0.005,0}},{\[Pi],"\[Pi]",{0.01,0}},{5\[Pi]/4,"",{0.005,0}},{3\[Pi]/2,"\!\(\*FractionBox[\(3  \[Pi]\), \(2\)]\)",{0.01,0}},{7\[Pi]/4,"",{0.005,0}},{2\[Pi],"2\[Pi]",{0.01,0}},{9\[Pi]/4,"",{0.005,0}},{5\[Pi]/2,"\!\(\*FractionBox[\(5  \[Pi]\), \(2\)]\)",{0.01,0}},{11\[Pi]/4,"",{0.005,0}},{3\[Pi],"3\[Pi]",{0.01,0}},{13\[Pi]/4,"",{0.005,0}},{7\[Pi]/2,"\!\(\*FractionBox[\(7  \[Pi]\), \(2\)]\)",{0.01,0}},{15\[Pi]/4,"",{0.005,0}},{4\[Pi],"4\[Pi]",{0.01,0}}},None};*)
(*(* Protect all variables *)*)
(*Protect[FrameTicksTenth,FrameTicksHalf,FrameTicksOne,FrameTicksPi8,FrameTicksPi4,FrameTicksPi2];*)
(*(* List of nice grid lines for the coupling constant g *)*)
(*Unprotect[GridLinesTenth,GridLinesHalf,GridLinesOne,GridLinesPi8,GridLinesPi4,GridLinesPi2];*)
(*ClearAll[GridLinesTenth,GridLinesHalf,GridLinesOne,GridLinesPi8,GridLinesPi4,GridLinesPi2];*)
(*GridLinesTenth={0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};*)
(*GridLinesHalf={0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10};*)
(*GridLinesOne={0,1,2,3,4,5,6,7,8,9,10};*)
(*GridLinesPi8={0,\[Pi]/8,\[Pi]/4,3\[Pi]/8,\[Pi]/2,5\[Pi]/8,3\[Pi]/4,7\[Pi]/8,\[Pi],9\[Pi]/8,5\[Pi]/4,11\[Pi]/8,3\[Pi]/2,13\[Pi]/8,7\[Pi]/4,15\[Pi]/8,2\[Pi]};*)
(*GridLinesPi4={0,\[Pi]/4,\[Pi]/2,3\[Pi]/4,\[Pi],5\[Pi]/4,3\[Pi]/2,7\[Pi]/4,2\[Pi],9\[Pi]/4,5\[Pi]/2,11\[Pi]/4,3\[Pi],13\[Pi]/4,7\[Pi]/2,15\[Pi]/4,4\[Pi]};*)
(*GridLinesPi2={0,\[Pi]/2,\[Pi],3\[Pi]/2,2\[Pi],5\[Pi]/2,3\[Pi],7\[Pi]/2,4\[Pi]};*)
(*(* Protect all variables *)*)
(*Protect[GridLinesTenth,GridLinesHalf,GridLinesOne,GridLinesPi8,GridLinesPi4,GridLinesPi2];*)


(* ::Input:: *)
(*(* List of nice frame ticks for the mass m *)*)
(*Unprotect[FrameTicks500,FrameTicks1000,FrameTicks2000];*)
(*ClearAll[FrameTicks500,FrameTicks1000,FrameTicks2000];*)
(*FrameTicks500={{{0,0,{0.01,0}},{250,"",{0.005,0}},{500,500,{0.01,0}},{750,"",{0.005,0}},{1000,1000,{0.01,0}},{1250,"",{0.005,0}},{1500,1500,{0.01,0}},{1750,"",{0.005,0}},{2000,2000,{0.01,0}},{2250,"",{0.005,0}},{2500,2500,{0.01,0}},{2750,"",{0.005,0}},{3000,3000,{0.01,0}},{3750,"",{0.005,0}},{3500,3500,{0.01,0}},{3750,"",{0.005,0}},{4000,4000,{0.01,0}},{4250,"",{0.005,0}},{4500,4500,{0.01,0}},{4750,"",{0.005,0}},{5000,5000,{0.01,0}}},None};*)
(*FrameTicks1000={{{0,0,{0.01,0}},{500,"",{0.005,0}},{1000,1000,{0.01,0}},{1500,"",{0.005,0}},{2000,2000,{0.01,0}},{2500,"",{0.005,0}},{3000,3000,{0.01,0}},{3500,"",{0.005,0}},{4000,4000,{0.01,0}},{4500,"",{0.005,0}},{5000,5000,{0.01,0}},{5500,"",{0.005,0}},{6000,6000,{0.01,0}},{6500,"",{0.005,0}},{7000,7000,{0.01,0}},{7500,"",{0.005,0}},{8000,8000,{0.01,0}},{8500,"",{0.005,0}},{9000,9000,{0.01,0}},{9500,"",{0.005,0}},{10000,10000,{0.01,0}}},None};*)
(*FrameTicks2000={{{0,0,{0.01,0}},{1000,"",{0.005,0}},{2000,2000,{0.01,0}},{3000,"",{0.005,0}},{4000,4000,{0.01,0}},{5000,"",{0.005,0}},{6000,6000,{0.01,0}},{7000,"",{0.005,0}},{8000,8000,{0.01,0}},{9000,"",{0.005,0}},{10000,10000,{0.01,0}}},None};*)
(*(* Protect all variables *)*)
(*Protect[FrameTicks500,FrameTicks1000,FrameTicks2000];*)
(*(* List of nice grid lines for the mass m *)*)
(*Unprotect[GridLines500,GridLines1000,GridLines2000];*)
(*ClearAll[GridLines500,GridLines1000,GridLines2000];*)
(*GridLines500={0,500,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000,6500,7000,7500,8000,8500,9000,9500,10000};*)
(*GridLines1000={0,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000};*)
(*GridLines2000={0,2000,4000,6000,8000,10000};*)
(*(* Protect all variables *)*)
(*Protect[GridLines500,GridLines1000,GridLines2000];*)


(* ::Input:: *)
(*(* List of nice frame ticks for percentages *)*)
(*Unprotect[FrameTicks20Percent];*)
(*ClearAll[FrameTicks20Percent];*)
(*FrameTicks20Percent={{{0,"0%",{0.01,0}},{10,"",{0.005,0}},{20,"20%",{0.01,0}},{30,"",{0.005,0}},{40,"40%",{0.01,0}},{50,"",{0.005,0}},{60,"60%",{0.01,0}},{70,"",{0.005,0}},{80,"80%",{0.01,0}},{90,"",{0.005,0}},{100,"100%",{0.01,0}},{110,"",{0.005,0}},{120,"120%",{0.01,0}},{130,"",{0.005,0}},{140,"140%",{0.01,0}},{150,"",{0.005,0}},{160,"160%",{0.01,0}},{170,"",{0.005,0}},{180,"180%",{0.01,0}},{190,"",{0.005,0}},{200,"200%",{0.01,0}}},None};*)
(*(* Protect all variables *)*)
(*Protect[FrameTicks20Percent];*)
(*(* List of nice grid lines for percentages *)*)
(*Unprotect[GridLines20Percent];*)
(*ClearAll[GridLines20Percent];*)
(*GridLines20Percent={0,20,40,60,80,100,120,140,160,180,200};*)
(*(* Protect all variables *)*)
(*Protect[GridLines20Percent];*)


(* ::Input:: *)
(*(* List of nice frame ticks for \[Chi] distribution plots *)*)
(*Unprotect[FrameTicks\[Chi],FrameTicksTwentieth];*)
(*ClearAll[FrameTicks\[Chi],FrameTicksTwentieth];*)
(*FrameTicks\[Chi]={{{0,"0",{0.01,0}},{1,"",{0.005,0}},{2,"",{0.005,0}},{3,"",{0.005,0}},{4,"",{0.005,0}},{5,"5",{0.01,0}},{6,"",{0.005,0}},{7,"",{0.005,0}},{8,"",{0.005,0}},{9,"",{0.005,0}},{10,"10",{0.01,0}},{11,"",{0.005,0}},{12,"",{0.005,0}},{13,"",{0.005,0}},{14,"",{0.005,0}},{15,"15",{0.01,0}},{16,"",{0.005,0}},{17,"",{0.005,0}},{18,"",{0.005,0}},{19,"",{0.005,0}},{20,"20",{0.01,0}},{21,"",{0.005,0}},{22,"",{0.005,0}},{23,"",{0.005,0}},{24,"",{0.005,0}},{25,"25",{0.01,0}},{26,"",{0.005,0}},{27,"",{0.005,0}},{28,"",{0.005,0}},{29,"",{0.005,0}},{30,"30",{0.01,0}}},None};*)
(*FrameTicksTwentieth={{{0.0,"0.00",{0.01,0}},{0.025,"",{0.005,0}},{0.05,"0.05",{0.01,0}},{0.075,"",{0.005,0}},{0.1,"0.10",{0.01,0}},{0.125,"",{0.005,0}},{0.15,"0.15",{0.01,0}},{0.175,"",{0.005,0}},{0.2,"0.20",{0.01,0}},{0.225,"",{0.005,0}},{0.25,"0.25",{0.01,0}},{0.275,"",{0.005,0}},{0.3,"0.30",{0.01,0}},{0.325,"",{0.005,0}},{0.35,"0.35",{0.01,0}},{0.375,"",{0.005,0}},{0.4,"0.40",{0.01,0}},{0.425,"",{0.005,0}},{0.45,"0.45",{0.01,0}},{0.475,"",{0.005,0}},{0.5,"0.50",{0.01,0}}},None};*)
(*(* Protect all variables *)*)
(*Protect[FrameTicks\[Chi],FrameTicksTwentieth];*)
(*(* List of nice grid lines for percentages *)*)
(*Unprotect[GridLines\[Chi],GridLinesTwentieth];*)
(*ClearAll[GridLines\[Chi],GridLinesTwentieth];*)
(*GridLines\[Chi]={0,5,10,15,20,25,30};*)
(*GridLinesTwentieth={0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5};*)
(*(* Protect all variables *)*)
(*Protect[GridLines\[Chi],GridLinesTwentieth];*)


(* ::Subsection:: *)
(*Frame Styles & Grid Lines*)


(* ::Input::Initialization:: *)
(* Usage Descriptions *)
GridLinesStep::usage="GridLinesStep[step, min, max]: Use to create custum grid lines, to be used with GridLines plot setting.";
FrameTicksStep::usage="FrameTicksStep[step, min, max, nolabels(optional]: Use to create custom frame ticks, to be used with FrameTicks plot setting.";
GridLinesStepLog::usage="GridLinesStepLog[multiplystep, min, max]: Use to create custum grid lines on a log scale, to be used with GridLines plot setting.";
FrameTicksStepLog::usage="FrameTicksStepLog[multiplystep, min, max, nolabels(optional)]: Use to create custom frame ticks on a log scale, to be used with FrameTicks plot setting.";

(* Begin Definitions *)
Begin["`Private`"];

GridLinesStep[step_,min_,max_]:=Table[gridline,{gridline,min,max,step}];

Default[FrameTicksStep,4]=False;
FrameTicksStep[step_,min_,max_,nolabels_.]:={Prepend[Delete[#,0]&/@Table[{{tick,If[!nolabels,If[Abs[Mod[tick,1]]<=$MachineEpsilon,NumberForm[tick,{Infinity,1}],tick],"",""],{0.01,0}},{tick+step/2,"",{0.005,0}}},{tick,min,max,step}],{min-step/2,"",{0.005,0}}],None};

GridLinesStepLog[multiplystep_,min_,max_]:=Module[{steps},
steps=Ceiling[Log[multiplystep,max/min]];
Table[min*multiplystep^step,{step,0,steps,1}]
];

Default[FrameTicksStepLog,4]=False;
FrameTicksStepLog[multiplystep_,min_,max_,nolabels_.]:=Module[{steps},
steps=Ceiling[Log[multiplystep,max/min]];
{Delete[#,0]&/@Table[FlattenAt[{{min*multiplystep^step,If[!nolabels,TextString[min*multiplystep^step],"",Superscript[10,Log10[min*multiplystep^step]]],{0.01,0}},Table[{min*multiplystep^step*(1+(multiplystep-1)/9 substep),"",{0.005,0}},{substep,1,8,1}]},2],{step,0,steps,1}],None}
];

(* End Definitions *)
End[];

(* Protect & Distribute Definitions *)
DistributeDefinitions[GridLinesStep,FrameTicksStep,GridLinesStepLog,FrameTicksStepLog];
Protect[GridLinesStep,FrameTicksStep,GridLinesStepLog,FrameTicksStepLog];


(* ::Subsection::Closed:: *)
(*Legends*)


(* ::Input::Initialization:: *)
(* Usage Descriptions *)
LineStyle::usage="LineStyle[color]: Creates a line style out of the given color.";
MakeLegend::usage="MakeLegend[color, style, area(optional)]: Creates a plot legend entry for the given color, style and area type (area == 0,1,2,3,4,5).";

(* Begin Definitions *)
Begin["`Private`"];

LineStyle[color_]:=Directive[Thick,color];

Default[MakeLegend,3]=0;
MakeLegend[color_,style_,area_.]:=Module[{},
(* Enclosed area legend *)
If[area==5,
Return[
Plot[{1/6,5/6},
{x,0,1},
Filling->{1->{2}},
PlotRange->{{0,1},{0,1}},
PlotStyle->{Directive[style,LineStyle[color]],Directive[style,LineStyle[color]]},
FillingStyle->Directive[Opacity[1/3],color],
Frame->False,Axes->False]
];
];
(* Enclosed area legend with middle line *)
If[area==4,
Return[
Plot[{1/2,1/6,5/6},
{x,0,1},
Filling->{2->{3}},
PlotRange->{{0,1},{0,1}},
PlotStyle->{Directive[style,LineStyle[color]],Directive[Dashed,LineStyle[color]],Directive[Dashed,LineStyle[color]]},
FillingStyle->Directive[Opacity[1/3],color],
Frame->False,Axes->False]
];
];
(* High area legend *)
If[area==3,
Return[
Plot[{1/6,1/6,5/6},
{x,0,1},
Filling->{2->{3}},
PlotRange->{{0,1},{0,1}},
PlotStyle->{Directive[style,LineStyle[color]],None,None},
FillingStyle->Directive[Opacity[1/3],color],
Frame->False,Axes->False]
];
];
(* Double area legend *)
If[area==2,
Return[
Plot[{1/2,1/3,2/3,1/6,5/6},
{x,0,1},
Filling->{2->{3},4->{5}},
PlotRange->{{0,1},{0,1}},
PlotStyle->{Directive[style,LineStyle[color]],None,None,None,None},
FillingStyle->Directive[Opacity[1/6],color],
Frame->False,Axes->False]
];
];
(* Single area legend *)
If[area==1,
Return[
Plot[{1/2,1/6,5/6},
{x,0,1},
Filling->{2->{3}},
PlotRange->{{0,1},{0,1}},
PlotStyle->{Directive[style,LineStyle[color]],None,None},
FillingStyle->Directive[Opacity[1/3],color],
Frame->False,Axes->False]
];
];
(* No area legend *)
Return[
Plot[1/2,
{x,0,1},
PlotRange->{{0,1},{0,1}},
PlotStyle->Directive[style,LineStyle[color]],
Frame->False,Axes->False]
];
];

(* End Definitions *)
End[];

(* Protect & Distribute Definitions *)
DistributeDefinitions[LineStyle,MakeLegend];
Protect[LineStyle,MakeLegend];


(* ::Section:: *)
(*End Package*)


(* ::Input::Initialization:: *)
(* End the package, functions are protected at their respective locations *)
EndPackage[];


(* ::Input::Initialization:: *)
(* Package information displayed on end of loading *)
If[PU$PackageLogging===True,
Print["PlotUtility version "<>PU$Version<>" loaded."];
];
