#' @examples
#' # Test for R Algorithm
#'
#' # usage test
#' data.usage <- system.file("extdata", "4465.csv", package = "nilm15hz")
#' app.usage  <- app.usage.main(data.file = data.usage,
#'                            app.meta.json = meta.list.jason.4465)
#'
#' myplot(    app.usage)
#' test.usage(app.usage, file.for.usage = data.usage) # An output of TRUE means OK
#'
#'
#' # detect and usage tests (takes 30 ~ 40 mins)
#' mytest <- TestNILM()
#' myplot(    mytest$app.usage)
#' test.usage(mytest$app.usage) # An output of TRUE means OK
#'
#'
#' # Test for system (from edev)
#' system("python CloudNilm.py -s 10002813 --REMOVE=all --DEBUG")
#' system("python CloudNilm.py -s 10002813 -j app_search -p 2015-12-20_00:00:00 2015-12-29_23:59:59 --DEBUG --PROCESS 1")
#' system("python CloudNilm.py -s 10002813 -j app_usage  -p 2016-01-11_00:00:00 2016-01-11_23:59:59 --DEBUG --PROCESS 1")
#'
#' # Read usage (cassandra)
#' setPAT = "640e2849d4b51ab2de08486e28dbc2985008ee8b"
#' devtools::install_github("EncoredTech/JediETL", auth_token = setPAT)
#' library(JediETL)
#' myresult <- dumpNilmApplianceUsageEncoredAPI(site_id = 10002813,
#'                                              power_unit = "W",
#'                                              start_date = "2016-01-11",
#'                                              end_date   = "2016-01-12")
#library(jsonlite)

test.usage <- function(app.usage,
                       file.for.usage = "/home/sulgik/data/15hz/7138/1450278000000-1450364399999.10005266_7138_3.csv"){

  out <- c()
  mysum <-
    colSums(ldply(app.usage, function(y){
      sapply(as.character(0:23), function(x) y[['hourly']][[x]])})[,-1])

  gold <- PreprocessNHz(file.for.usage)
  gold$timestamp <-as.POSIXct(gold$timestamp, tz = "Asia/Seoul", origin = "1970-01-01")

  goldHour <- hourlyF(gold, col = "active_power")

  test <- all(ldply(as.character(0:23), function(x) abs((mysum[[x]] - goldHour[[x]])/goldHour[[x]])) < .05)
  out <- c(out, test)
  return(out)
}

TestNILM <- function(file.for.meta  = "/home/sulgik/data/15hz/7138/1449068400000-1450277999999.10005266_7138_3.csv",
                     file.for.usage = "/home/sulgik/data/15hz/7138/1450278000000-1450364399999.10005266_7138_3.csv",
                     show.plot = FALSE){
  ptm = proc.time()
  #  if (usage.only){
  #    load(mymeta)
  #  } else {
  # "/home/sulgik/data/15hz/1450537200000-1451401199999.10004465_6299_3.csv"
  meta.list.json <- app.detect.main(data.file = file.for.meta, find.ac = F)
  #  }
  time1 <- proc.time() - ptm

  ptm = proc.time()
  app.usage <- app.usage.main(data.file = file.for.usage, app.meta.json = meta.list.json)
  proc.time() - ptm
  time2 <- proc.time() - ptm

  if (show.plot) myplot(app.usage)

  out <- list(meta.list.json = meta.list.json,
              app.usage      = app.usage,
              time.meta  = time1,
              time.usage = time2)
  class(out) <- "TestNILM"
  return(out)
}

print.TestNILM <- function(x, ...) {
  if(!inherits(x, "TestNILM")) stop("ERROR: This function ONLY supports the TestNILM object.")
  cat(x$app.usage)
  return(invisible(x))
}


myplot <- function(app.usage){
  mydata <- as.data.frame(sapply(app.usage, function(x) sapply(0:23, function(y) x$hourly[[as.character(y)]])))
  names(mydata) <- c(sapply(app.usage, function(x) x$shape_type))
  mydata$hour <- 0:23

  mydata.melt <- melt(mydata, id.vars = "hour")
  mydata.melt %>% ggplot(aes(x = hour, y = value, fill = variable)) +
    geom_bar(stat = "identity")
}
