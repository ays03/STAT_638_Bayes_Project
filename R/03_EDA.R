library(naniar)
library(tsibble)
library(feasts)
library(zoo)

### Missingness
# Count missing values per column
missing_summary <- colSums(is.na(df))

# Visualize missingness (deterministic)
gg_miss_var(df)

# Missingness over time
df %>%
  mutate(missing_gap = is.na(Global_active_power)) %>%
  ggplot(aes(x = DateTime, y = missing_gap)) +
  geom_point(alpha = 0.1) +
  labs(y = "Missing (TRUE/FALSE)")



### Descriptive Stats
stats <- df %>%
  summarise(across(
    where(is.numeric),
    list(
      mean = ~mean(.x, na.rm = TRUE),
      median = ~median(.x, na.rm = TRUE),
      sd = ~sd(.x, na.rm = TRUE),
      min = ~min(.x, na.rm = TRUE),
      max = ~max(.x, na.rm = TRUE)
    )
  ))


### raw time series visualization
df %>%
  ggplot(aes(x = DateTime, y = Global_active_power)) +
  geom_line(alpha = 0.4) +
  labs(y = "Global Active Power (kW)")


### Seasonal Decomposition
df_ts <- df %>%
  as_tsibble(index = DateTime)

df_ts %>%
  model(STL(Global_active_power ~ season(window = "periodic"))) %>%
  components() %>%
  autoplot()


### Distributions
# Histogram
df %>%
  ggplot(aes(x = Global_active_power)) +
  geom_histogram(bins = 100, fill = "steelblue", color = "white")

# Log histogram (only positive values)
df %>%
  filter(Global_active_power > 0) %>%
  ggplot(aes(x = log(Global_active_power))) +
  geom_histogram(bins = 100, fill = "tomato", color = "white")

### Autocorrelation
df_ts %>%
  ACF(Global_active_power) %>%
  autoplot()

df_ts %>%
  PACF(Global_active_power) %>%
  autoplot()

### Outlier detection
df <- df %>%
  mutate(
    gap_z = (Global_active_power - mean(Global_active_power, na.rm = TRUE)) /
      sd(Global_active_power, na.rm = TRUE),
    is_outlier = abs(gap_z) > 3
  )

df %>%
  filter(is_outlier) %>%
  ggplot(aes(x = DateTime, y = Global_active_power)) +
  geom_point(color = "red")

### Need to feature engineer an energy column, we might also consider adding a 
# column (dummy var) indicating a holiday season to see if it is correalted
# with energy consumption and can help improve the model.

