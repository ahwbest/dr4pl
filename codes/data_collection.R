
k <- 6
data.no <- 32

conc.ref=conc.ref.sorted[k]

df.each.conc=df.plot[df.plot$RefDose==conc.ref,]

data.err <- df.each.conc
data.err <- subset(data.err, select = c(Dose, Response))
colnames(data.err) <- c("Dose", "Response")

setwd("C:\\Users\\Hyowon\\Research\\Dose_response_modelling\\Data")
write.csv(data.err, paste("Dirk_data_", data.no, ".csv", sep = ""), row.names = FALSE)

drm.drc <- drm(Response ~ Dose,
               data = data.err,
               fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "IC50")))

plot(drm.drc,
     type = "all",
     pch = 16,
     lwd = 2,
     bty = "l")