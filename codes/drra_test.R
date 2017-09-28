# ------------------------------------------------------------------------------
### Test code
#
setwd("C:\\Users\\Hyowon\\Research\\Dose_response_modelling\\Data")

data.org <- read.csv("Dirk_data_2.csv")

drra.org <- drra(dose = data.org$Dose, response = data.org$Response, n.init = 100)
drra.form <- drra(Response ~ Dose, data = data.org, n.init = 100)

