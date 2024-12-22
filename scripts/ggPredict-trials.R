
library(tidyverse)
library(ggiraphExtra)

# hp = x
# drat = colors
model1 <- lm(mpg ~ hp*drat, data = mtcars)

ggPredict(model1, interactive = TRUE)

# disp = x
# drat = facets
# hp = colors
model2 <- lm(mpg ~ disp + hp + drat + hp*drat, data = mtcars)

ggPredict(model2, interactive = TRUE)
