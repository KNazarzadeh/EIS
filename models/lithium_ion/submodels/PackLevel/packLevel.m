% function packLevel (battery, AgingParams,  parallelCellNum, seriesCellNum, cycleTotalNum)

load("CellLevel.mat")

cycleTotalNum = 11;

PackLevelParams.applied_current = vertcat(CellLevel.applied_current{2:cycleTotalNum, 1});
PackLevelParams.Voltage = vertcat(CellLevel.Voltage{2:cycleTotalNum, 1});
PackLevelParams.total_capacity_loss = vertcat(CellLevel.total_capacity_loss{2:cycleTotalNum, 1});

PackLevelParams.applied_current = PackLevelParams.applied_current * parallelCellNum;
PackLevelParams.Voltage = PackLevelParams.Voltage * seriesCellNum;
PackLevelParams.total_capacity_loss = PackLevelParams.total_capacity_loss * parallelCellNum;


% end
