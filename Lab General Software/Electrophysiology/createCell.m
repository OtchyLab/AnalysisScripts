function cell = createCell(exper, cellNum, locationDesc, electrodeType, PosAP, PosLM, PosDV, sutterDepth, numLatTurns, signalNum, chanNum, bAntidromic)

cell.exper = exper;
cell.num = cellNum;
cell.locDesc = locationDesc;
cell.electrodeType = electrodeType;
cell.signalNum = signalNum;
cell.chanNum = chanNum;
cell.AP = PosAP;
cell.LM = PosLM;
cell.DV = PosDV;
cell.latPos = numLatTurns;
cell.sutterDepth = sutterDepth;
cell.bAntidromic = bAntidromic;