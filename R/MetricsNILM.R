#' @examples
#' # Evaluation from cassandra
#' library(encoredLog)
#' myvalidOut <- validateNilmApplianceUsage(key = c("백창현", "최재필"),
#'                                          start_date = "2016-01-23",
#'                                          end_date = "2016-01-24")
#' (mymet <- MetricsNILM(myvalidOut))
#' summaryF(mymet)
#'



MetricsNILM <-
  function(validate.out, threshold = 5)
  {
    # key <- "백창현"; start_date = "2016-01-19"; end_date = "2016-01-20"; threshold = 10
    # validate.out <- validateNilmApplianceUsage(key = key, start_date = start_date, end_date = end_date)
    if (is.data.frame(validate.out) & nrow(validate.out) > 0) {
      out <- validate.out %>%
        mutate(scanner_on = factor(scanner_usage > threshold, levels = c(T, F)),
               nilm_on    = factor(   nilm_usage > threshold, levels = c(T, F)) ) %>%
        plyr::dlply(.(siteApp),
                    function(x) c(rawdata = x, with(x, list(appliance = unique(appliance),
                                                            confusion.matrix = table(scanner_on, nilm_on)))))
      metricsF <- function(x){
        data.frame(appliance = x$appliance,
                   daily.nilm = sum(x$rawdata.nilm_usage,    na.rm = T),
                   daily.scnr = sum(x$rawdata.scanner_usage, na.rm = T),
                   accuracy.amt = max(0, 1 - abs(sum(x$rawdata.nilm_usage, na.rm = T) - sum(x$rawdata.scanner_usage, na.rm = T)) /
                     sum(x$rawdata.scanner_usage, na.rm = T)),
                   accuracy  = (x$confusion.matrix["TRUE", "TRUE"] + x$confusion.matrix["FALSE", "FALSE"]) /
                     sum(x$confusion.matrix),
                   recall    =  x$confusion.matrix["TRUE", "TRUE"] / sum(x$confusion.matrix["TRUE", ]),
                   precision =  x$confusion.matrix["TRUE", "TRUE"] / sum(x$confusion.matrix[, "TRUE"]),
                   scannerOn_nilmOn   = x$confusion.matrix["TRUE" , "TRUE" ],
                   scannerOn_nilmOff  = x$confusion.matrix["TRUE" , "FALSE"],
                   scannerOff_nilmOn  = x$confusion.matrix["FALSE", "TRUE" ],
                   scannerOff_nilmOff = x$confusion.matrix["FALSE", "FALSE"])}

      plyr::ldply(out, metricsF)
    } else {stop("nodata")}
  }

summaryF <- function(metricsFout, max.par = "accuracy.amt") {
  # myvalidOut  <- validateNilmApplianceUsage(key = c("백창현"), start_date = "2016-01-11", end_date = "2016-01-12")
  # metricsFout <- MetricsNILM(myvalidOut)

  metricsFout %>%
    ddply(.(appliance), function(x) x[which.max(x[[max.par]])[1], ])
}



# scanner <- myscanner("최재필", start_date = "2015-01-23", end_date = "2015-01-24")

#validateNilmApplianceUsage <-
#  function(app.usage,
#           key = NULL, start_date, end_date,
#           time_unit  = getTimeUnitSet('nilm'),
#           power_unit = c('mW', 'W', 'kW'),
#           sampling = 15, shape_type1 = "washingMachine")
#  {
    # app.usage = s$app.usage; key = c("백창현"); start_date = '2015-12-17'; end_date = '2015-12-18'; time_unit = "hourly"; power_unit='W'; sampling = 15; shape_type1 = "washingMachine", shape_type2 ="세탁기"
#
#     require(encoredDataETL)
#     time.unit     <- time_unit
#     power.unit    <- power_unit
#
#
#     meta.test.set <-
#       searchMetaData(key, type = 'nilm_test') %>%
#       select(key, site_id, user, scanner_serial_hex, appliance, powerConsumption)
#     stopifnot(nrow(meta.test.set) > 0)
#
#     compare.table.list <- list()
#
#     scanner <- llply(1:nrow(meta.test.set), function(x){
#         return( dumpScannerPeriodicUsage(
#           meta.test.set[x,]$user,
#           meta.test.set[x,]$appliance,
#           start_date, end_date,
#           time_unit = "hourly", power_unit = "W"))
#     })
#
#     # remove null from q
#
#     outer(1:length(NILM_output), 1:length(app.usage), function(x, y){
#
#       scanner_usage <- scanner[[x]]
#       nilm_usage    <- unlist(s$app.usage[[y]][['hourly']])
#
#       scanner_usage$on <- scanner_usage$usage > 10
#       nilm_usage$on    <- nilm_usage > 10
#
#     })
#
#
#
#
#         comparison <- function(given_appliance){
#
#        qappliance == s
#     }
#
#         mapply(q, app.usage, function())
#
#
#
#
#     for (i in 1:nrow(meta.test.set)) {
#       # i = 2
#       print(paste('Key:', i))
#       user      <- meta.test.set$user[i]
#       site.id   <- meta.test.set$site_id[i]
#       appliance <- meta.test.set$appliance[i]
#
#       if (is.na(meta.test.set$scanner_serial_hex[i])) {
#         print(paste('No scanner data: index is', i))
#         next()
#       }
#
#       tmp.nilm <- s$app.usage
#
#       target.app.nilm <- subset(tmp.nilm, typeName == appliance)
#
#
#
#
#       tmp.scanner
#
#       compare.table <-
#         tmp.scanner %>%
#         left_join(target.app.nilm, by = 'date') %>%
#         select(user, appliance, appliance_id, date, scanner_usage, nilm_usage) %>%
#         mutate(
#           delta      = scanner_usage - nilm_usage,
#           error_rate = ifelse(
#             scanner_usage == 0 | nilm_usage == 0,
#             0,
#             ifelse(
#               scanner_usage > nilm_usage,
#               abs(delta) / nilm_usage * 100,
#               abs(delta) / scanner_usage * 100
#             )
#           )
#         )
#       compare.table.list[[i]] <- compare.table
#     }
#
#     result <- ldply(compare.table.list)
#     return(result)
#   }
#


myscanner <- function(key = NULL, start_date, end_date){
  meta.test.set <-
    searchMetaData(key, type = 'nilm_test') %>%
    select(key, site_id, user, scanner_serial_hex, appliance, powerConsumption)
  stopifnot(nrow(meta.test.set) > 0)

  meta.test.set <- meta.test.set[!is.na(meta.test.set$scanner_serial_hex), , drop = F]

  output <-
    llply(1:nrow(meta.test.set), function(x){
      dumpScannerPeriodicUsage(
        meta.test.set[x, ]$user,
        meta.test.set[x, ]$appliance,
        start_date = start_date,
        end_date   = end_date,
        time_unit = "hourly", power_unit = "W")
    })
  return(output)
}


