devtools::install_github('blue-matter/SLAM')

library(SLAM)

Data <- Import(Example_Data()[1])

# Report(Data)

MyAssess <- Assess(Data)


Report(MyAssess)
