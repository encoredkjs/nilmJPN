
app.usage.test <- function(
  test.apps=c("kimchi-cyclic","kimchi-pattern","microwave-hslee-1", "ac-hslee-1","ac-knbae-1"), show.fig=T){
  
  #
  #
  # Library
  #
  library(NILM1Hz)
  library(plyr)
  library(reshape)
  library(lubridate)
  library(evaluateNILM)
  library(ggplot2)
  library(RJSONIO)
  
  # from AC
  library( lubridate )
  library( plyr )
  library( ggplot2 )
  library( depmixS4 )
  library( hsmm )
  
  
  ################################################################
  #
  # --------------------  사용량 코드 테스트 -----------------------
  #
  #               사용 데이터셋: 클로즈 베타 그룹
  #
  ################################################################
  
  correct.usage.list = list(
    'kimchi-cyclic'=list('usage'=442.2242, 'file'='/home/hslee/test_data/meta_test_data/test_data_hslee_20150906.csv', 'meta'='
{
  "kimchi-cyclic" : {
  "meta-version" : 1,
               "shape_type" : "cyclic_box",
               "rising_edge" : {
               "ap_height" : 161.56,
               "ap_sigma" : 2.9721,
               "rp_height" : 55.853,
               "rp_sigma" : 1.9206
               },
               "falling_edge" : {
               "ap_height" : -134.1,
               "ap_sigma" : 1.8373,
               "rp_height" : -49.893,
               "rp_sigma" : 0.58444
               },
               "cycle" : {
               "working_time" : 84.938,
               "working_time_sigma" : 5.5674,
               "duty_cycle" : 0.12012
               },
               "box.no" : 18,
               "generation_info" : {
               "data_used" : {
               "start" : "2015-10-06 00:00:02",
               "end" : "2015-10-07 00:00:00",
               "sampling" : 1
               },
               "computed" : "2015-10-07 17:41:42"
               }
  }
}
'),
  'kimchi-pattern'=list('usage'=442.2242, 'file'='/home/hslee/test_data/meta_test_data/test_data_hslee_20150906.csv', 'meta'='
{        
  "kimchi-pattern" : {
  	"meta-version" : 1,
		"shape_type" : "pattern_scan",
                        "rising_edge" : {
                        "summit_flag" : 1,
                        "cluster" : 3,
                        "ap_min" : 132,
                        "ap_max" : 187.44,
                        "ap_med" : 161.75,
                        "rp_min" : 47.499,
                        "rp_max" : 60.756,
                        "rp_med" : 53.769,
                        "sum" : 107,
                        "min.t" : 25.023,
                        "med.t" : 687.96,
                        "min.med.rate" : 0.036373,
                        "lost.sig.num" : 24,
                        "lost.sig.rate" : 18.321
                        },
                        "falling_edge" : {
                        "cluster" : 2,
                        "ap_min" : -161.2,
                        "ap_max" : -74.805,
                        "ap_med" : -131.94,
                        "rp_min" : -53.705,
                        "rp_max" : -43.315,
                        "rp_med" : -48.338,
                        "slotNum.zero" : 7,
                        "slotNum.one" : 85,
                        "ZerotoOneratio" : 0.082353,
                        "EffTimeOn.med" : 89.96,
                        "EffTimeOn.min" : 80.021,
                        "EffTimeOn.max" : 123.95,
                        "EffTimeOn.sd" : 5.0891,
                        "EffAP_Drop.med" : -132.35,
                        "EffAP_Drop.min" : -152.06,
                        "EffAP_Drop.max" : -116.11,
                        "EffAP_Drop.sd" : 5.6348,
                        "EffRP_Drop.med" : -48.33,
                        "EffRP_Drop.min" : -53.525,
                        "EffRP_Drop.max" : -44.564,
                        "EffRP_Drop.sd" : 1.5692
                        },
                        "supportRatio" : 0.99261,
                        "generation_info" : {
                        "data_used" : {
                        "start" : "2015-09-06 00:00:02",
                        "end" : "2015-09-07 00:00:00",
                        "sampling" : 1
                        },
                        "computed" : "2015-10-27 01:36:52"
                        },
                        "parameters" : {
                        "genResolution.n" : 20,
                        "genEffSize.n" : 15,
                        "staPeriodicity.p" : 0.1,
                        "endEffSlot.p" : 0.1,
                        "endConsistency.n" : 1,
                        "min_sigmag.active.n" : 20,
                        "min_sigmag.reactive.n" : 3,
                        "min_lossrate.p" : 0.8,
                        "clustering.method" : 3
                        }
}
}
'),
  'microwave-hslee-1'=list('usage'=271.9535, 'file'='/home/hslee/test_data/meta_test_data/test_data_hslee_20150815-20150819.csv', 'meta'='
{
  "microwave-hslee-1" :  {
  	"meta-version" : 1,
		"shape_type" : "pattern_scan_heavy",
                           "summit_flag" : 1,
                           "rising_edge" : {
                           "ap_min" : 1156,
                           "ap_median" : 1180.1,
                           "ap_max" : 1212.9,
                           "rp_min" : 322.37,
                           "rp_median" : 348.82,
                           "rp_max" : 362.37
                           },
                           "falling_edge" : {
                           "ap_min" : -1174.8,
                           "ap_median" : -1153.1,
                           "ap_max" : -1146.5,
                           "rp_min" : -338.99,
                           "rp_median" : -307.71,
                           "rp_max" : -306.39
                           },
                           "duration" : {
                           "duration_min" : 21.001,
                           "duration_median" : 51.489,
                           "duration_max" : 92.058
                           },
                           "n.pt" : 6,
                           "generation_info" : {
                           "data_used" : {
                           "start" : "2015-08-15 00:00:02",
                           "end" : "2015-08-20 00:00:00",
                           "sampling" : 1
                           },
                           "computed" : "2015-10-20 17:29:00"
                           }
}
}'),
  'ac-hslee-1'=list('usage'=16342.54, 'file'='/home/hslee/test_data/meta_test_data/test_data_hslee_20150815-20150819.csv', 'meta'='
{        
  "ac-hslee-1" : {
  "meta-version" : 1,
  "shape_type" : "HMM"
  }
}'),
  'ac-knbae-1'=list('usage'=7091.619, 'file'='/home/hslee/test_data/meta_test_data/test_data_knbae_20150815-20150819.csv', 'meta'='
{        
  "ac-knbae-1" : {
  "meta-version" : 1,
  "shape_type" : "HMM"
  }
}')
)
  

  usage.answer = data.frame()
  for( f in (test.apps) ){
    test.info = correct.usage.list[[f]]
    test.file <-  test.info[['file']]
    test.meta <-  test.info[['meta']]
    test.usage <-  test.info[['usage']]
    
    print(paste('calculate for ',f))
    
    #  
    # 
    #source('nilm_app_usage_main.R')
    usage.out <- app.usage.main(test.file, test.meta, find.heavy=T, find.pattern=T, find.rice=T, find.cyclic=T, find.ac=T, find.standby=F, show.fig=show.fig)
    #
    #
    usage.out[[f]][['daily']]
    
    #
    # 데이터 저장
    #
    usage.answer <- rbind(usage.answer,data.frame(app_id=f,usage.out=usage.out[[f]][['daily']],usage.correct=test.usage))
  }

  print('---------------------- Answer Sheet -----------------------')
  print(usage.answer)
  print('---------------------- End of Answer -----------------------')

  print('---------------------- Incorrect Answers -----------------------')
  print(usage.answer[ abs(usage.answer$usage.correct-usage.answer$usage.out)/usage.answer$usage.correct > 0.20, ])
  print('---------------------- End of Incorrect Answers -----------------------')

  #
  # 메타 코드 테스트 끝!
  ###################################################################################
}
