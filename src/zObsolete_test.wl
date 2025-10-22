
(* --- Clear all variables --- *)    
ClearAll["Global`*"];


(* --- Parameters --- *)
R = 1;
r = 0.5;

(* --- Reaction functions --- *)
reactionTrp[Trp_, IPA_imine_, IPA_imine_dimer_, CPA_, PDVA_, DVA_, PDV_, DV_, PVA_, PV_, VA_, Violacein_, x_] := -0.75 * VioA[x] * Trp / (0.000125 + Trp);
reactionIPA_imine[Trp_, IPA_imine_, IPA_imine_dimer_, CPA_, PDVA_, DVA_, PDV_, DV_, PVA_, PV_, VA_, Violacein_, x_] := -0.75 * VioB[x] * IPA_imine / (0.000125 + IPA_imine)+0.75 * VioA[x] * Trp / (0.000125 + Trp);
reactionIPA_imine_dimer[Trp_, IPA_imine_, IPA_imine_dimer_, CPA_, PDVA_, DVA_, PDV_, DV_, PVA_, PV_, VA_, Violacein_, x_] := -0.75 * VioE[x] * IPA_imine_dimer / (0.000125 + IPA_imine_dimer)-0.001 * IPA_imine_dimer+0.75 * VioB[x] * IPA_imine / (0.000125 + IPA_imine);
reactionCPA[Trp_, IPA_imine_, IPA_imine_dimer_, CPA_, PDVA_, DVA_, PDV_, DV_, PVA_, PV_, VA_, Violacein_, x_] := +0.001 * IPA_imine_dimer;
reactionPDVA[Trp_, IPA_imine_, IPA_imine_dimer_, CPA_, PDVA_, DVA_, PDV_, DV_, PVA_, PV_, VA_, Violacein_, x_] := -0.75 * VioC[x] * PDVA / (0.000125 + PDVA)-0.75 * VioD[x] * PDVA / (0.000125 + PDVA)-0.001 * PDVA+0.75 * VioE[x] * IPA_imine_dimer / (0.000125 + IPA_imine_dimer);
reactionDVA[Trp_, IPA_imine_, IPA_imine_dimer_, CPA_, PDVA_, DVA_, PDV_, DV_, PVA_, PV_, VA_, Violacein_, x_] := -0.001 * DVA+0.75 * VioC[x] * PDVA / (0.000125 + PDVA);
reactionPDV[Trp_, IPA_imine_, IPA_imine_dimer_, CPA_, PDVA_, DVA_, PDV_, DV_, PVA_, PV_, VA_, Violacein_, x_] := +0.001 * PDVA;
reactionDV[Trp_, IPA_imine_, IPA_imine_dimer_, CPA_, PDVA_, DVA_, PDV_, DV_, PVA_, PV_, VA_, Violacein_, x_] := +0.001 * DVA;
reactionPVA[Trp_, IPA_imine_, IPA_imine_dimer_, CPA_, PDVA_, DVA_, PDV_, DV_, PVA_, PV_, VA_, Violacein_, x_] := -0.75 * VioC[x] * PVA / (0.000125 + PVA)-0.001 * PVA+0.75 * VioD[x] * PDVA / (0.000125 + PDVA);
reactionPV[Trp_, IPA_imine_, IPA_imine_dimer_, CPA_, PDVA_, DVA_, PDV_, DV_, PVA_, PV_, VA_, Violacein_, x_] := +0.001 * PVA;
reactionVA[Trp_, IPA_imine_, IPA_imine_dimer_, CPA_, PDVA_, DVA_, PDV_, DV_, PVA_, PV_, VA_, Violacein_, x_] := -0.001 * VA+0.75 * VioC[x] * PVA / (0.000125 + PVA);
reactionViolacein[Trp_, IPA_imine_, IPA_imine_dimer_, CPA_, PDVA_, DVA_, PDV_, DV_, PVA_, PV_, VA_, Violacein_, x_] := +0.001 * VA;
(* --- Enzyme concentration and location --- *)
VioA[x_] := Piecewise[{{ 0.25, 0.1 < x < 0.2 }, { 0.25, 0.8 < x < 0.9 }}]
VioB[x_] := Piecewise[{{ 0.25, 0.6 < x < 0.7 }}]
VioC[x_] := Piecewise[{{ 0.025, 0.05 < x < 0.15 }}]
VioD[x_] := Piecewise[{{ 0.025, 0.2 < x < 0.3 }}]
VioE[x_] := Piecewise[{{ 0.25, 0.4 < x < 0.5 }}]

(* --- Left domain solver: x in [0,r] *)
leftSolver = ParametricNDSolveValue[
{
-6.6e-10*Trp''[x] + reactionTrp[Trp[x], IPA_imine[x], IPA_imine_dimer[x], CPA[x], PDVA[x], DVA[x], PDV[x], DV[x], PVA[x], PV[x], VA[x], Violacein[x], x] == 0,
Trp[0] == Trp0,
Trp'[0] == 0,
-6.6e-10*IPA_imine''[x] + reactionIPA_imine[Trp[x], IPA_imine[x], IPA_imine_dimer[x], CPA[x], PDVA[x], DVA[x], PDV[x], DV[x], PVA[x], PV[x], VA[x], Violacein[x], x] == 0,
IPA_imine[0] == IPA_imine0,
IPA_imine'[0] == 0,
-6.6e-10*IPA_imine_dimer''[x] + reactionIPA_imine_dimer[Trp[x], IPA_imine[x], IPA_imine_dimer[x], CPA[x], PDVA[x], DVA[x], PDV[x], DV[x], PVA[x], PV[x], VA[x], Violacein[x], x] == 0,
IPA_imine_dimer[0] == IPA_imine_dimer0,
IPA_imine_dimer'[0] == 0,
-6.6e-10*CPA''[x] + reactionCPA[Trp[x], IPA_imine[x], IPA_imine_dimer[x], CPA[x], PDVA[x], DVA[x], PDV[x], DV[x], PVA[x], PV[x], VA[x], Violacein[x], x] == 0,
CPA[0] == CPA0,
CPA'[0] == 0,
-6.6e-10*PDVA''[x] + reactionPDVA[Trp[x], IPA_imine[x], IPA_imine_dimer[x], CPA[x], PDVA[x], DVA[x], PDV[x], DV[x], PVA[x], PV[x], VA[x], Violacein[x], x] == 0,
PDVA[0] == PDVA0,
PDVA'[0] == 0,
-6.6e-10*DVA''[x] + reactionDVA[Trp[x], IPA_imine[x], IPA_imine_dimer[x], CPA[x], PDVA[x], DVA[x], PDV[x], DV[x], PVA[x], PV[x], VA[x], Violacein[x], x] == 0,
DVA[0] == DVA0,
DVA'[0] == 0,
-6.6e-10*PDV''[x] + reactionPDV[Trp[x], IPA_imine[x], IPA_imine_dimer[x], CPA[x], PDVA[x], DVA[x], PDV[x], DV[x], PVA[x], PV[x], VA[x], Violacein[x], x] == 0,
PDV[0] == PDV0,
PDV'[0] == 0,
-6.6e-10*DV''[x] + reactionDV[Trp[x], IPA_imine[x], IPA_imine_dimer[x], CPA[x], PDVA[x], DVA[x], PDV[x], DV[x], PVA[x], PV[x], VA[x], Violacein[x], x] == 0,
DV[0] == DV0,
DV'[0] == 0,
-6.6e-10*PVA''[x] + reactionPVA[Trp[x], IPA_imine[x], IPA_imine_dimer[x], CPA[x], PDVA[x], DVA[x], PDV[x], DV[x], PVA[x], PV[x], VA[x], Violacein[x], x] == 0,
PVA[0] == PVA0,
PVA'[0] == 0,
-6.6e-10*PV''[x] + reactionPV[Trp[x], IPA_imine[x], IPA_imine_dimer[x], CPA[x], PDVA[x], DVA[x], PDV[x], DV[x], PVA[x], PV[x], VA[x], Violacein[x], x] == 0,
PV[0] == PV0,
PV'[0] == 0,
-6.6e-10*VA''[x] + reactionVA[Trp[x], IPA_imine[x], IPA_imine_dimer[x], CPA[x], PDVA[x], DVA[x], PDV[x], DV[x], PVA[x], PV[x], VA[x], Violacein[x], x] == 0,
VA[0] == VA0,
VA'[0] == 0,
-6.6e-10*Violacein''[x] + reactionViolacein[Trp[x], IPA_imine[x], IPA_imine_dimer[x], CPA[x], PDVA[x], DVA[x], PDV[x], DV[x], PVA[x], PV[x], VA[x], Violacein[x], x] == 0,
Violacein[0] == Violacein0,
Violacein'[0] == 0,
},
{ Trp, IPA_imine, IPA_imine_dimer, CPA, PDVA, DVA, PDV, DV, PVA, PV, VA, Violacein },
{ x, 0, r },
{ Trp0, IPA_imine0, IPA_imine_dimer0, CPA0, PDVA0, DVA0, PDV0, DV0, PVA0, PV0, VA0, Violacein0}
],


(* --- Right domain solver: x in [r,R] *)
rightSolver = ParametricNDSolveValue[
{
-6.6e-10*Trp''[x] + reactionTrp[Trp[x], IPA_imine[x], IPA_imine_dimer[x], CPA[x], PDVA[x], DVA[x], PDV[x], DV[x], PVA[x], PV[x], VA[x], Violacein[x], x] == 0,
Trp[R] == TrpR,
Trp'[R] == 9e-05/6.6e-10 * (2.5e-05 - Trp[R]),
-6.6e-10*IPA_imine''[x] + reactionIPA_imine[Trp[x], IPA_imine[x], IPA_imine_dimer[x], CPA[x], PDVA[x], DVA[x], PDV[x], DV[x], PVA[x], PV[x], VA[x], Violacein[x], x] == 0,
IPA_imine[R] == IPA_imineR,
IPA_imine'[R] == 9e-05/6.6e-10 * (0.0 - IPA_imine[R]),
-6.6e-10*IPA_imine_dimer''[x] + reactionIPA_imine_dimer[Trp[x], IPA_imine[x], IPA_imine_dimer[x], CPA[x], PDVA[x], DVA[x], PDV[x], DV[x], PVA[x], PV[x], VA[x], Violacein[x], x] == 0,
IPA_imine_dimer[R] == IPA_imine_dimerR,
IPA_imine_dimer'[R] == 9e-05/6.6e-10 * (0.0 - IPA_imine_dimer[R]),
-6.6e-10*CPA''[x] + reactionCPA[Trp[x], IPA_imine[x], IPA_imine_dimer[x], CPA[x], PDVA[x], DVA[x], PDV[x], DV[x], PVA[x], PV[x], VA[x], Violacein[x], x] == 0,
CPA[R] == CPAR,
CPA'[R] == 9e-05/6.6e-10 * (0.0 - CPA[R]),
-6.6e-10*PDVA''[x] + reactionPDVA[Trp[x], IPA_imine[x], IPA_imine_dimer[x], CPA[x], PDVA[x], DVA[x], PDV[x], DV[x], PVA[x], PV[x], VA[x], Violacein[x], x] == 0,
PDVA[R] == PDVAR,
PDVA'[R] == 9e-05/6.6e-10 * (0.0 - PDVA[R]),
-6.6e-10*DVA''[x] + reactionDVA[Trp[x], IPA_imine[x], IPA_imine_dimer[x], CPA[x], PDVA[x], DVA[x], PDV[x], DV[x], PVA[x], PV[x], VA[x], Violacein[x], x] == 0,
DVA[R] == DVAR,
DVA'[R] == 9e-05/6.6e-10 * (0.0 - DVA[R]),
-6.6e-10*PDV''[x] + reactionPDV[Trp[x], IPA_imine[x], IPA_imine_dimer[x], CPA[x], PDVA[x], DVA[x], PDV[x], DV[x], PVA[x], PV[x], VA[x], Violacein[x], x] == 0,
PDV[R] == PDVR,
PDV'[R] == 9e-05/6.6e-10 * (0.0 - PDV[R]),
-6.6e-10*DV''[x] + reactionDV[Trp[x], IPA_imine[x], IPA_imine_dimer[x], CPA[x], PDVA[x], DVA[x], PDV[x], DV[x], PVA[x], PV[x], VA[x], Violacein[x], x] == 0,
DV[R] == DVR,
DV'[R] == 9e-05/6.6e-10 * (0.0 - DV[R]),
-6.6e-10*PVA''[x] + reactionPVA[Trp[x], IPA_imine[x], IPA_imine_dimer[x], CPA[x], PDVA[x], DVA[x], PDV[x], DV[x], PVA[x], PV[x], VA[x], Violacein[x], x] == 0,
PVA[R] == PVAR,
PVA'[R] == 9e-05/6.6e-10 * (0.0 - PVA[R]),
-6.6e-10*PV''[x] + reactionPV[Trp[x], IPA_imine[x], IPA_imine_dimer[x], CPA[x], PDVA[x], DVA[x], PDV[x], DV[x], PVA[x], PV[x], VA[x], Violacein[x], x] == 0,
PV[R] == PVR,
PV'[R] == 9e-05/6.6e-10 * (0.0 - PV[R]),
-6.6e-10*VA''[x] + reactionVA[Trp[x], IPA_imine[x], IPA_imine_dimer[x], CPA[x], PDVA[x], DVA[x], PDV[x], DV[x], PVA[x], PV[x], VA[x], Violacein[x], x] == 0,
VA[R] == VAR,
VA'[R] == 9e-05/6.6e-10 * (0.0 - VA[R]),
-6.6e-10*Violacein''[x] + reactionViolacein[Trp[x], IPA_imine[x], IPA_imine_dimer[x], CPA[x], PDVA[x], DVA[x], PDV[x], DV[x], PVA[x], PV[x], VA[x], Violacein[x], x] == 0,
Violacein[R] == ViolaceinR,
Violacein'[R] == 9e-05/6.6e-10 * (0.0 - Violacein[R]),
},
{ Trp, IPA_imine, IPA_imine_dimer, CPA, PDVA, DVA, PDV, DV, PVA, PV, VA, Violacein },
{ x, r, R },
{ TrpR, IPA_imineR, IPA_imine_dimerR, CPAR, PDVAR, DVAR, PDVR, DVR, PVAR, PVR, VAR, ViolaceinR}
],
