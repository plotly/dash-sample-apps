data("Smarket")

Smarket$Direction = as.factor(Smarket$Direction)
names(Smarket)
glm.fit <- glm(Direction ~ Lag1 + Lag2 + Lag3 + Lag4 + 
               Lag5 + Volume, data = Smarket, 
               family = binomial)


train = Smarket[Smarket$Year < 2005,]
test = Smarket[!(Smarket$Year < 2005),]

glm.fit <- glm(Direction ~ Lag1 + Lag2 + Lag3 + Lag4 + Lag5 + Volume, 
               data = Smarket, 
               family = binomial, subset = Smarket$Year<2005)

summary(glm.fit)

glm.probs <- predict(glm.fit, 
                     newdata = test, 
                     type = "response")

glm.pred <- ifelse(glm.probs > 0.5, "Up", "Down")

Direction.2005 = test$Direction

table(glm.pred, Direction.2005)

mean(glm.pred == Direction.2005)


#we dont have a good accuracy. lets improve by using a smaller
#model since lag 4 and 5 dont have good p values.


glm.fit = glm(formula = Direction ~ Lag1 + Lag2 + Lag3 + Volume + Year, family = binomial, 
    data = Smarket)

summary(glm.fit)

# 
# glm.probs <- predict(glm.fit, 
#                      newdata = test, 
#                      type = "response")
# 
# glm.pred <- ifelse(glm.probs > 0.5, "Up", "Down")
# 
# Direction.2005 = test$Direction
# 
# table(glm.pred, Direction.2005)
# 
# mean(glm.pred == Direction.2005)
# 

newdata = data.frame(Lag1 = 0.381, Lag2 = -0.192, Lag3 = -2.624, 
                     Volume = 2, Year = 2020)

predict(glm.fit, newdata, type = 'response')

