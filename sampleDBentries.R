# Database examples
library(RSQLite)
library(ggplot2)

falkorDb <- dbConnect(drv = RSQLite::SQLite(), "falkor.db")
dbListTables(falkorDb)
bet_df <- dbGetQuery(falkorDb, "SELECT * FROM raw_data_smol WHERE mz>? AND mz<?", 
                     params=pmppm(117.078979+1.007276))
dbDisconnect(falkorDb)

ggplot(bet_df, aes(x=rt, y=int, group=file)) + geom_line()
