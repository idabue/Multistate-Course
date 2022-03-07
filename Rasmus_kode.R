library(survival)
library(Epi)
library(popEpi)

data(cancer)
data(lung)
lung$sex <- factor(lung$sex,levels=1:2,labels=c("M","W"))
lung$time <- lung$time/(365.25/12)
head(lung)
km <- survfit(Surv(time,status==2)~1,data=lung)
km
summary(km)
plot(km)
km_sex <- survfit(Surv(time,status==2)~sex,data=lung)
plot(km_sex,col=1:2)
plot(km_sex, col = c("blue", "red"), lwd = 1, conf.int = TRUE)
lines(km_sex, col = c("blue", "red"), lwd = 3)
with(lung,tapply(age,sex,mean))
survdiff(Surv(time,status==2)~sex,data=lung)

c0 <- coxph(Surv(time,status==2)~sex,data=lung)
c1 <- coxph(Surv(time,status==2)~sex+age,data=lung)
c2 <- coxph(Surv(time,status==2)~sex+I(age/10),data=lung)
ci.exp(c0)
ci.exp(c1)
ci.exp(c2)
summary(c1)


p1 <- glm(cbind(status==2,time)~sex+age,
          family=poisreg,
          data=lung)
px <- glm(status == 2 ~ sex + age + offset(log(time)),
          family = poisson,
          data = lung)
ci.exp(p1)
ci.exp(c2)
ci.exp(px)

L1 <-  Lexis(exit = list(tfl = time),
              exit.status = factor(status,
                                     levels = 1:2,
                                     labels = c("Alive","Dead")),
              data = lung)
L1_cox <- coxph(Surv(tfl,tfl + lex.dur,lex.Xst == "Dead") ~ sex + age,
             data = L1)
boxes(L1, boxpos = TRUE)

cL <- coxph.Lexis(L1, tfl ~ sex + age)
pL <- glm.Lexis(L1, ~ sex + age)

cL
pL

S1 <- splitMulti(L1,tfl=0:36)
summary(L1)
summary(S1)
wh <- names(L1)[1:10] # names of variables in some order
subset(L1, lex.id == 10)[,wh]
subset(S1,lex.id==10)[,wh]
ps <- glm(cbind(lex.Xst == "Dead", lex.dur) ~ Ns(tfl, knots = seq(0, 36, 12)) + sex + age,
           family = poisreg,
           data = S1)
pc <- glm(cbind(lex.Xst == "Dead", lex.dur) ~ sex + age,
           family = poisreg,
           data = L1)

ci.exp(ps)
summary(ps)

round(cbind(ci.exp(cl),
            ci.exp(ps, subset = c("sex","age")),
            ci.exp(pc, subset = c("sex","age"))), 3)

plot(ps)
data(DMlate)
# str(DMlate)
set.seed(1952)
DMlate <- DMlate[sample(1:nrow(DMlate), 2000),]
str(DMlate)
Ldm <- Lexis(entry = list(per = dodm,
                          age = dodm - dobth,
                          tfd = 0),
              exit = list(per = dox),
              exit.status = factor(!is.na(dodth),
                                     labels = c("DM","Dead")),
              data = DMlate)
summary(Ldm)

Cdm <- cutLexis(Ldm,
                 cut = Ldm$dooad,
                 timescale = "per",
                 new.state = "OAD")
summary(Cdm)
Adm <- subset(Cdm, lex.Cst == "DM")
summary(Adm)
boxes(Adm,
      boxpos = TRUE,
      scale.R = 100,
      show.BE = T)
boxes(Cdm,
      boxpos = TRUE,
      scale.R = 100,
      show.BE = T)

m3 <- survfit(Surv(tfd,
                    tfd + lex.dur,
                    lex.Xst) ~ 1,
               id = lex.id,
               data = Adm)

head(cbind(time = m3$time, m3$pstate))
tail(cbind(time = m3$time, m3$pstate))
names(m3)
m3$transitions
summary(Adm)
m3$pstate
par( mfrow=c(1,2) )
 matplot(m3$time, m3$pstate,
           type="s", lty=1, lwd=4,
           col=c("ForestGreen","red","black"),
           xlim=c(0,15), xaxs="i",
           ylim=c(0,1), yaxs="i" )
stackedCIF(m3, lwd=3, xlim=c(0,15), xaxs="i", yaxs="i" )
text( rep(12,3), c(0.9,0.3,0.6), levels(Cdm) )
box()

 
 
Sdm <- splitMulti(Adm, tfd = seq(0,20,0.1) )
summary(Adm) 
round(cbind(
with(subset(Sdm, lex.Xst == "OAD" ), quantile(tfd + lex.dur, 0:10/10)),
with(subset(Sdm, lex.Xst == "Dead"), quantile(tfd + lex.dur, 0:10/10))),3)
okn <- c(0,0.5,3,6)
dkn <- c(0,2.0,5,9)
OAD.glm <- glm.Lexis(Sdm, ~ Ns(tfd, knots = okn), from = "DM", to = "OAD" )
Dead.glm <- glm.Lexis(Sdm, ~ Ns(tfd, knots = dkn), from = "DM", to = "Dead")

int <- 0.01
nd <- data.frame(tfd = seq(0, 15, int))
l.glm <- ci.pred( OAD.glm, nd)
m.glm <- ci.pred(Dead.glm, nd)
matshade(nd$tfd,
          cbind(l.glm, m.glm) * 100,
          plot = TRUE,
          log = "y", ylim = c(2, 20),
          col = rep(c("red","black"), 2), lwd = 3)

cR <- ci.Crisk(mods = list(OAD = OAD.glm,
                           Dead = Dead.glm), nd = nd)

matshade(as.numeric(dimnames(cR$Crisk)[[1]]),
         cbind(cR$Crisk[,1,],
               cR$Crisk[,2,],
               cR$Crisk[,3,]), plot = TRUE,
         lwd = 2, col = c("limegreen","red","black"))

mat2pol(cR$Crisk[,3:1,1], col = c("forestgreen","red","black")[3:1])
matshade(as.numeric(dimnames(cR$Srisk)[[1]]),
          cbind(cR$Srisk[,1,],
                cR$Srisk[,2,]), plot = TRUE,
         lwd = 2, col = c("black","red"),
         ylim = 0:1, yaxs = "i")

