Simulation <- Simulate(LifeHistory, Exploitation, nsim=2)
Sampled_Data <- Generate_Data(Simulation, Sampling = Sampling)
Data <- Import_Data(Sampled_Data, sim=1)
Write_Data2XL(Simulation, xlfile='inst/Data_Example_Binned.xlsx', binned=TRUE, dev=T)
Write_Data2XL(Simulation, xlfile='inst/Data_Example_Raw.xlsx', dev=T)

