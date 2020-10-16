library(readxl)
library(vars)
afp <- read_excel("G:/Mi unidad/10mo Semestre/Series de Tiempo II/Tareas/afp.xlsx")
View(afp)

attach(afp)
data <- data.frame(Cuprum.B,Habitat.B,PlanVital.B,ProVida.B)

afpts <- ts(data,frequency = 12,start=c(2005,8))

#A
plot.ts(afpts)

#B
VARselect(afpts, lag.max=5, type="const")

#C
Acoef(afpfit)

afpfit <- VAR(afpts, lag.max=10, type="const")
A_df <- data.frame(Acoef(afpfit))
A_est <- cbind(A_df[,1],A_df[,2],A_df[,3],A_df[,4])
A_est

#D
stability(afpfit)

#E
newdata <- data[-c(174,175,176,177),]

newafpts <- ts(newdata,frequency = 12,start=c(2005,8))

newafpfit <- VAR(newafpts,lag.max=10, type="const")

Acoef(newafpfit)

predict(newafpfit,n.ahead=4,ci=0.95)

ecmc <- ((2.99-0.29)**2+(-3.83-0.62)**2+(-13.74-0.45)**2+(7.68-0.41)**2)/4
ecmh <- ((3.22-0.33)**2+(-3.52-0.64)**2+(-13.13-0.48)**2+(7.55-0.45)**2)/4
ecmpr <- ((3.09-0.21)**2+(-3.78-0.60)**2+(-13.39-0.44)**2+(7.91-0.41)**2)/4
ecmpl <- ((2.97-0.29)**2+(-3.79-0.60)**2+(-13.70-0.44)**2+(8.19-0.40)**2)/4

#F
plot(predict(newafpfit,n.ahead=4,ci=0.95))

#G
normality.test(afpfit)
