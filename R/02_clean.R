
# Load packages explicitly
library(tidyverse)
library(lubridate)


# Set reproducible locale for dates
Sys.setlocale("LC_TIME", "C")

# Set reproducible timezone
Sys.setenv(TZ = "UTC")

# Define path explicitly
path <- "/data/household_power_consumption.txt"

#import data with explicit column types tp aovid R from guessing
df <- read_delim(
  file = path,
  delim = ";",
  na = "?",
  col_types = cols(
    Date = col_character(),
    Time = col_character(),
    Global_active_power = col_double(),
    Global_reactive_power = col_double(),
    Voltage = col_double(),
    Global_intensity = col_double(),
    Sub_metering_1 = col_double(),
    Sub_metering_2 = col_double(),
    Sub_metering_3 = col_double()
  )
)

# parse data-time deterministically (no locale guessing)
df <- df %>%
  mutate(
    DateTime = dmy_hms(
      paste(Date, Time),
      tz = "UTC"   # reproducible timezone
    )
  ) %>%
  arrange(DateTime)



