library(astsa)
library(rugarch)
library(tseries)

ts.plot(sales)

sales.dif1 <- diff(sales)
sales.dif1.sq <- sales.dif1**2

ts.plot(sales.dif1)

spec <- ugarchspec()
fit <- ugarchfit(spec = spec, data = sales.dif1)

show(fit)
summary(fit)

acf(sales.dif1)
acf(sales.dif1.sq)

pacf(sales.dif1)
pacf(sales.dif1.sq)

hist(sales.dif1)

boxplot(sales.dif1)

qqnorm(sales.dif1)

Box.test(sales.dif1, type="Ljung")
Box.test(sales.dif1.sq, type="Ljung")
