## We use DATA4 created in ReadData.R


DATA4 <- combine_DataFiles(L1 = DATA2, L2 = DATA1)
str(DATA4)
DATA4$ddot_env

P = Palaeodose_Computation(DATA4, Nb_sample = 2, Iter = 1000, SampleNames = DATA4$SampleNames)

create_MeasuresDataFrame(P, DATA4, .02, .3)
