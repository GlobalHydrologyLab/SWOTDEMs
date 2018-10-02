function [comp] = relativeSteepnessComp(skm,z,nSlopes,steepMaskTrue)

[Rd,~,sId] = relativeSteepness(skm,z,nSlopes);
simCutoff = nanmean(Rd) + nanstd(Rd);
comp.steepMask = Rd>simCutoff;
comp.steepId = sId(comp.steepMask);

comp.correct = steepMaskTrue + comp.steepMask == 2;
comp.falseNeg = steepMaskTrue - comp.steepMask == 1;
comp.falsePos = steepMaskTrue - comp.steepMask == -1;

comp.Stats = [sum(comp.correct) sum(comp.falseNeg) sum(comp.falsePos)];

end

