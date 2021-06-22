# Load and fix dates on the historical tick data
# This does not need to be run, as the datasets have been saved as an R object

library(data.table)

EURUSD <- fread(file.path("pairs", "EURUSD.csv.gz"))
USDCHF <- fread(file.path("pairs", "USDCHF.csv.gz"))
USDJPY <- fread(file.path("pairs", "USDJPY.csv.gz"))
GBPUSD <- fread(file.path("pairs", "GBPUSD.csv.gz"))

USDCHF[, Date := strptime(Date, format = "%Y%m%d %H:%M:%OS")]
USDJPY[, Date := strptime(Date, format = "%Y%m%d %H:%M:%OS")]
GBPUSD[, Date := strptime(Date, format = "%Y%m%d %H:%M:%OS")]

save(EURUSD, USDCHF, USDJPY, GBPUSD, file = "pairs/currencyPairs.RData")
