#####
# Initial model
# Poisson GLM
#####
data <- mydataY %>% 
  select(ANNEE = ADEC, SEXE = SEXE, AGE = AGE2, obs = dcf, OFFSET = POPf)
res <- glm(obs ~ ANNEE * factor(AGE) * SEXE, 
           offset = log(OFFSET), 
           family = poisson, 
           data = data %>% filter(ANNEE %in% 2010:2019))
newdata <- data %>% 
  filter(ANNEE %in% 2010:2023) %>% 
  select(ANNEE, AGE, SEXE, obs, OFFSET)
estimations_glm <- predict(res, newdata = newdata, type = "response", se.fit = TRUE)
pred2 <- cbind(newdata, val_pred = estimations_glm$fit, se = estimations_glm$se.fit) %>% 
  mutate(diff = obs - val_pred)
# Adding Years of Life Lost (YLL)
data_pred <- left_join(
  pred2, 
  esp2019 %>% select(AGE, SEXE, ESP19 = ESP)
) %>% 
  mutate(yll = diff * ESP19)
# Bootstrap for confidence intervals on annual predictions
set.seed(1)
num_iterations <- 10000
results <- vector("list", num_iterations)
sim <- i <- table <- test <- NULL
for (i in 1:num_iterations) {
  test <- data_pred %>% filter(ANNEE %in% 2020:2023)
  sim <- rnorm(nrow(test), mean = test$val_pred, sd = test$se)
  test$boot <- sim
  test <- test %>%
    mutate(esp = (obs - boot) * ESP19) %>%
    group_by(ANNEE, SEXE) %>%
    summarise(
      time = i,
      diff = sum(obs - val_pred),
      diff_e = sum(obs - boot),
      yll = sum(esp),
      .groups = 'drop'
    )
  results[[i]] <- test
  if (i %% 1000 == 0) {print(i)}
}
table <- bind_rows(results)
# Total excess deaths (ED) and YLL over 4 years
table %>% 
  ungroup() %>% 
  group_by(time) %>% 
  summarise(
    diff = sum(diff),
    diff_e = sum(diff_e),
    yll = sum(yll)
  ) %>% 
  summarise(
    ed = mean(diff),
    LCL = quantile(diff_e, 0.025),
    UCL = quantile(diff_e, 0.975),
    YLL = mean(yll),
    lcl = quantile(yll, 0.025),
    ucl = quantile(yll, 0.975)
  )
# General table with YLL and ED
left_join(
  table %>% 
    group_by(ANNEE, time) %>% 
    summarise(
      diff = sum(diff),
      diff_e = sum(diff_e),
      yll = sum(yll)
    ) %>% 
    group_by(ANNEE = ANNEE) %>% 
    summarise(
      obs_mean = mean(diff),
      boot_mean = mean(diff_e),
      LCL = quantile(diff_e, 0.025),
      UCL = quantile(diff_e, 0.975),
      YLL = mean(yll),
      lcl = quantile(yll, 0.025),
      ucl = quantile(yll, 0.975)
    ),
  data_pred %>% 
    filter(ANNEE %in% 2020:2024) %>% 
    group_by(ANNEE = ANNEE) %>% 
    summarise(dct = sum(val_pred))
) %>% 
  mutate(prop = obs_mean / dct * 100)
# Median YLL per individual
data_pred %>% 
  filter(ANNEE %in% 2020:2023 & diff >= 0) %>% 
  group_by(ANNEE) %>% 
  summarise(median = weighted_median_iqr(ESP19, diff))
# Bootstrap for proportion under 60 years
set.seed(1)
num_iterations <- 10000
results <- vector("list", num_iterations)
sim <- i <- table <- test <- NULL
for (i in 1:num_iterations) {
  test <- data_pred %>% filter(ANNEE %in% 2020:2023)
  sim <- rnorm(nrow(test), mean = test$val_pred, sd = test$se)
  test$boot <- sim
  test <- test %>%
    mutate(esp = (obs - boot) * ESP19)
  test <- left_join(
    test %>% 
      filter(AGE < 60) %>% 
      group_by(ANNEE, SEXE) %>% 
      summarise(yll59 = sum(esp), .groups = 'drop'),
    test %>% 
      group_by(ANNEE, SEXE) %>% 
      summarise(yll = sum(esp), .groups = 'drop'),
    by = c("ANNEE", "SEXE")
  ) %>% 
    mutate(prop = yll59 / yll * 100, time = i)
  results[[i]] <- test
  if (i %% 1000 == 0) {print(i)}
}
table_yll <- bind_rows(results)
# Proportion of YLL under 60 with CI
table_yll %>% 
  group_by(ANNEE, time) %>% 
  summarise(
    yll59 = sum(yll59),
    yll = sum(yll)
  ) %>% 
  mutate(prop = yll59 / yll * 100) %>% 
  group_by(ANNEE) %>% 
  summarise(
    pr60 = round(mean(prop), 1),
    lcl = round(quantile(prop, 0.025), 1),
    ucl = round(quantile(prop, 0.975), 1)
  ) %>% 
  mutate(ic = paste0(pr60, "% [", lcl, "%; ", ucl, "%]"))
#####
# Spline Model
#####
# Fit a spline-based Poisson regression model
res <- glm(
  dcf ~ SEXE * factor(AGE2) + ns(WEXE, df = 8) + ns(date, df = 3), 
  offset = log(POPf), 
  family = poisson, 
  data = mydataW %>% filter(ADEC %in% 2010:2019)
)
# Predictions using the spline model with standard errors
newdata <- mydataW %>%
  ungroup() %>%
  filter(ADEC %in% 2010:2023) %>%
  select(-dcm, -POPm)
estimations_glm <- predict(res, newdata = newdata, type = "response", se.fit = TRUE)
# Add predictions and residuals to the dataset
newdata$pred <- estimations_glm$fit
newdata$se <- estimations_glm$se
newdata <- newdata %>%
  mutate(diff = dcf - pred)
# Calculate excess mortality
newdata %>%
  filter(ADEC %in% 2020:2023) %>%
  group_by(ADEC) %>%
  summarise(
    diff = sum(diff), 
    pred = sum(pred), 
    prop = diff / pred * 100, 
    obs = sum(dcf)
  )
# Plot observed and predicted values using the spline model
newdata %>%
  group_by(ADEC) %>%
  summarise(
    diff = sum(diff, na.rm = TRUE), 
    pred = sum(pred), 
    prop = diff / pred * 100, 
    obs = sum(dcf)
  ) %>%
  ggplot(aes(x = ADEC)) +
  # Line for observed values
  geom_line(aes(y = obs), color = "black", linetype = "dotted", size = 1.2) +
  # Line for predicted values
  geom_line(aes(y = pred), color = "green", size = 1) +
  # Labels and themes
  labs(
    title = "Observed and Predicted Values",
    x = "Year",
    y = "Values"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12)
  )
#####
# Loop over all periods to calculate indicators with the spline model
#####
# Prepare data for the model
data <- mydataW %>%
  select(obs = dcf, AGE = AGE2, SEXE, ANNEE = ADEC, WEXE, date, OFFSET = POPf)
esp2019 <- read_excel("~/data/population/esp2019_2024_05.xlsx")
result_all <- data.frame()
result2_all <- data.frame()
# Loop through different reference periods
for (train in 7:13) {
  # Define reference years
  start_year <- 2020 - train
  reference_years <- start_year:2019
  
  # Fit the spline model for the reference period
  res <- glm(
    obs ~ factor(AGE) * SEXE + ns(WEXE, df = 8) + ns(date, df = 3), 
    offset = log(OFFSET), 
    family = poisson, 
    data = data %>% filter(ANNEE %in% reference_years)
  )
  
  # Adjust the prediction period based on the reference period
  prediction_years <- ifelse(start_year > 2010, start_year, 2010):2023
  newdata <- data %>%
    filter(ANNEE %in% prediction_years) %>%
    select(WEXE, AGE, SEXE, date, OFFSET, obs, ANNEE)
  
  # Generate predictions
  estimations_glm <- predict(res, newdata = newdata, type = "response", se.fit = TRUE)
  pred2 <- cbind(newdata, val_pred = estimations_glm$fit, se = estimations_glm$se.fit) %>%
    mutate(diff = obs - val_pred)
  
  # Add Years of Life Lost (YLL) to predictions
  data_pred <- left_join(pred2, esp2019 %>% select(AGE, SEXE, ESP19 = ESP)) %>%
    mutate(yll = diff * ESP19)
  
  # Results by sex and total
  result_sexe <- data_pred %>%
    group_by(ANNEE, SEXE) %>%
    summarise(
      difft = sum(diff), 
      obst = sum(obs), 
      prop = difft / sum(val_pred) * 100, 
      yllt = sum(yll), 
      .groups = "drop"
    ) %>%
    mutate(ylli = yllt / difft, train = train)
  
  result_tot <- data_pred %>%
    group_by(ANNEE) %>%
    summarise(
      difft = sum(diff), 
      obst = sum(obs), 
      prop = difft / sum(val_pred) * 100, 
      yllt = sum(yll), 
      .groups = "drop"
    ) %>%
    mutate(SEXE = "tot", ylli = yllt / difft, train = train)
  
  result <- bind_rows(result_sexe, result_tot)
  result_all <- bind_rows(result_all, result)
  
  # Calculate contribution of individuals under 60
  result2 <- left_join(
    data_pred %>%
      filter(ANNEE %in% 2020:2023 & AGE < 60) %>%
      group_by(ANNEE) %>%
      summarise(yll60 = sum(yll), .groups = "drop"),
    data_pred %>%
      filter(ANNEE %in% 2020:2023) %>%
      group_by(ANNEE) %>%
      summarise(yllt = sum(yll), .groups = "drop")
  ) %>%
    mutate(prop = yll60 / yllt * 100, train = train)
  
  result2_all <- bind_rows(result2_all, result2)
  print(train)
}
result_all
result2_all