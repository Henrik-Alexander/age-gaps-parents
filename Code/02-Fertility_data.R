#######################################
# Purpose: Loading the birth data     #
# Author: Henrik-Alexander Schubert   #
# Date: 14.06.2023                    #
# E-Mail: schubert@demogr.mpg.de      #
# Pre-Requisites: full repository     #
#######################################

# Load the functions
source("Functions/packages.R")
source("Functions/functions.R")
source("Functions/graphics.R")


# Load the survival package
library(survival)
library(eha)
library(casebase)


### Clean the bio-birth data  ----------------------------------------

# Load the birth data
fert <- read_stata("SOEP_V36/Stata/biobirth.dta")

# Remove respondents that where no asked the question
fert <- fert |> filter(bioyear != -1 & gebjahr != -1)

# Filter men
fert <- fert |> filter(sex == 1)

# Remove unimportant variables
fert <- fert |> select(!starts_with("kidsex"))

# Make everything as double
fert <- fert |> mutate(across(where(is.factor), as.double))

# Make missing, where values are either -2 or -1
fert <- fert |> replace_with_na_all(condition = ~.x %in% c(-2, -1))

# Clean the names
names(fert) <- sub("(.*)(\\d{2})$", "\\1_\\2", names(fert))

# Make a life-course perspective
fert2 <- fert |> pivot_longer(cols = starts_with("kid"), names_pattern = "([a-z]*)_([0-9]*)", values_to = "Value", names_to = c("Variable", "Number"))

# Filter first births
fert2 <- fert2 |> filter(Number == "01")

# Pivot wider
fert2 <- fert2 |> pivot_wider(names_from = c(Variable, Number), values_from = Value)

# Create cohorts - split by 5 year groups
fert2 <- fert2 |> mutate(cohort = cut(gebjahr, breaks = seq(1900, 2020, by = 10), dig.lab = 4))

# Double check
fert2 <- fert2 |> filter(!is.na(gebjahr) & !is.na(bioyear))

### Prepare the survival data ------------------------------------

# Create an event and censoring variable
fert2 <- fert2 |> mutate(Event = if_else(is.na(kidgeb_01), 0, 1),
                         Censoring = if_else(Event == 0, bioyear - gebjahr, kidgeb_01 - gebjahr))

# Look at the survival times
with(fert2, Surv(Censoring, Event))

# Make the Kaplan-Meier
km <- survfit(Surv(Censoring, Event) ~ 1, conf.type = "log",
              conf.int = 0.95, type = "kaplan-meier", error = "greenwood",
              data = fert2)

# Print the output of the kaplan-meier
summary(km)

# Plot the kaplan meier
with(km, data.frame(time, n.risk, n.event, surv, n.censor, cumhaz, std.chaz, lower, upper)) |> 
  filter(time <= 50) |> 
  ggplot(aes(time, y=  surv, ymin = lower, ymax = upper)) +
  geom_line() +
  geom_ribbon(alpha = .3)

# Fit by cohort
km_coh <- survfit(Surv(Censoring, Event) ~ cohort, data = fert2, conf.int = 0.95, type = "kaplan-meier", error = "greenwood") 

# Plot
plot(km_coh)

### Fit the smoothed hazards -------------------------------------

# Select the central cohorts
cohorts <- c("(1940,1950]","(1950,1960]","(1960,1970]")

# Filter the data
fert2 <- fert2 |> filter(cohort %in% cohorts)

# Run the regression model
mod_cb <- fitSmoothHazard(Event ~ ns(log(Censoring), df = 5) + cohort,
                          data = fert2, 
                          time = "Censoring")

# Print the result
summary(mod_cb)

# Plot the result
plot_results <- plot(mod_cb,
                     hazard.params = list(xvar = "Censoring",
                                          by = "cohort",
                                          alpha = 0.10,
                                          ylab = "Hazard",
                                          plot = FALSE))

# Plot the result
plot_results$fit |>
  filter(Censoring >= 15 & Censoring <= 50 & cohort %in% cohorts) |> 
  ggplot(aes(Censoring, visregFit, group = cohort, colour = cohort)) +
  geom_line()

### Make parametric hazard models --------------------------------

# Exponential
exp <- par_surv(distribution = "exponential")

# Weibull
weib <- par_surv(distribution = "weibull")

# Gompertz
#gomp <- par_surv(distribution = "gompertz")

# Gaussian
gauss <- par_surv(distribution = "gaussian")

# Lognormal
lognor <- par_surv(distribution = "lognormal")

# Log-logistic
loglog <- par_surv(distribution = "loglogistic")

# Generalized gamma


# Make a result table
stargazer(exp, weib, gauss, lognor, loglog)


### Predictions --------------------------------------------------

# Create an example dataset
pred_data <- expand.grid(Censoring = 18:55, cohort = unique(fert2$cohort))

# Predict the result
predict(lognor, type = "uquantile")

### Discrete time model ------------------------------------------

# Create the prediction data
pred_data <- expand.grid(Censoring = 18:55, cohort = unique(fert2$cohort))

# Estimate a logistic regression
logist <- glm(Event ~ ns(Censoring, df = 5), data = fert2)

# Predict the results
pred_data$prediction <- predict(logist, pred_data)

# Plot the result
ggplot(pred_data, aes(Censoring, prediction)) +
  geom_line() #+
#  geom_ribbon(aes(ymin = lower, ymax = upper))

### A non-parametric approach ------------------------------------

# Create a start variable
fert2 <- fert2 |> mutate()

# Split the data
spell_data <- survSplit(fert2, cut = 15:55, end = "Censoring", event = "Event", start = "start")

# Estimate the exposures
exposures <- spell_data |> group_by(start, cohort) |> count()

# Count the events
births <- spell_data |> group_by(start, cohort) |> summarise(birth = sum(Event))

# Combine
unparametric <- inner_join(exposures, births) |> mutate(rate = birth / n)


# Plot the result
plot_raw <- unparametric |> filter(cohort %in% cohorts & start >= 18) |> 
  ggplot(aes(start, rate, colour = cohort, group = cohort, shape = cohort)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(expand = c(0, 0)) +
  ggtitle("Age-specific first-birth rates for men") +
  labs(caption = "Data = SOEP Wave 36") +
  ylab("Age-specific fertility rate (Parity 1)") +
  xalb("Age")

# Plot interpolated
plot_interpol <- unparametric |> filter(cohort %in% cohorts & start >= 18) |> 
  ggplot(aes(start, rate, colour = cohort, group = cohort, linetype = cohort)) +
  geom_smooth() +
  scale_x_continuous(expand = c(0, 0)) +
  ggtitle("Age-specific first-birth rates for men (smoothed)") +
  labs(caption = "Data = SOEP Wave 36") +
  ylab("Age-specific fertility rate (Parity 1)") +
  xalb("Age")

## Save the graphs
ggsave(plot_raw, filename = "Figures/raw_firstbirth_soep.pdf")
ggsave(plot_interpol, filename = "Figures/smooth_firstbirth_soep.pdf")

### Combine with background variables ----------------------------

# Join with the individual data
#data <- inner_join(df, fert, by = "persnr")


# Filter men
data <- data %>% filter(sex == 2)

# Plot the data
p <- ggplot(data, aes(x = sumkids)) +
  geom_bar(position = "dodge", mapping = aes(y = ..prop.., group = treatment)) +
  scale_x_continuous(limits = c(-0.5, 5))+
  facet_wrap(~ treatment)
