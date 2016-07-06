
app.detect.test <- function(
  test.sites=c("2813","2822","2880","2882","2894","2893","3700","3795","4033","4465"),
  test.heavy=F, test.patternHigh=F, test.pattern=F, test.pattern_extend=F,
  test.cyclic=F, test.rice=F, test.standby=F, test.ac=F){

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
  # --------------------  메타 코드 테스트 -----------------------
  #
  #               사용 데이터셋: 클로즈 베타 그룹
  #
  ################################################################
  correct.meta.cyclic.list = list(
    '2813'=data.frame(
      'meta_shape'='cyclic_box',
      'rising_ap'=c(122),
      'rising_rp'=c(33),
      'falling_ap'=c(118),
      'falling_rp'=c(27),
      'time_on'=c(0),
      stringsAsFactors=FALSE
    ), 
    '2822'=data.frame(
      'meta_shape'='cyclic_box',
      'rising_ap'=c(222),
      'rising_rp'=c(65),
      'falling_ap'=c(189),
      'falling_rp'=c(57),
      'time_on'=c(0),
      stringsAsFactors=FALSE
    ), 
    '2880'=data.frame(
      'meta_shape'='cyclic_box',
      'rising_ap'=c(85),
      'rising_rp'=c(0),
      'falling_ap'=c(866),
      'falling_rp'=c(0),
      'time_on'=c(0),
      stringsAsFactors=FALSE
    ), 
    '2882'=data.frame(
      'meta_shape'='cyclic_box',
      'rising_ap'=c(143),
      'rising_rp'=c(88),
      'falling_ap'=c(129),
      'falling_rp'=c(86),
      'time_on'=c(0),
      stringsAsFactors=FALSE
    ),   
    '2894'=data.frame(
      'meta_shape'='cyclic_box',
      'rising_ap'=c(155),
      'rising_rp'=c(66),
      'falling_ap'=c(141),
      'falling_rp'=c(63),
      'time_on'=c(0),
      stringsAsFactors=FALSE
    ),
    '2893'=data.frame(
      'meta_shape'='cyclic_box',
      'rising_ap'=c(194),
      'rising_rp'=c(138),
      'falling_ap'=c(148),
      'falling_rp'=c(129),
      'time_on'=c(0),
      stringsAsFactors=FALSE
    ),
    '3700'=data.frame(
      'meta_shape'='cyclic_box',
      'rising_ap'=c(109),
      'rising_rp'=c(5),
      'falling_ap'=c(96),
      'falling_rp'=c(2),
      'time_on'=c(0),
      stringsAsFactors=FALSE
    ),
    '4033'=data.frame(
      'meta_shape'='cyclic_box',
      'rising_ap'=c(103),
      'rising_rp'=c(26),
      'falling_ap'=c(91),
      'falling_rp'=c(19),
      'time_on'=c(0),
      stringsAsFactors=FALSE
    ),
    '4465'=data.frame(
      'meta_shape'='cyclic_box',
      'rising_ap'=c(180),
      'rising_rp'=c(19),
      'falling_ap'=c(175),
      'falling_rp'=c(11),
      'time_on'=c(0),
      stringsAsFactors=FALSE
    )
  )
  
  correct.meta.pattern.list = list(
    '2813'=data.frame(
      'meta_shape'='pattern_scan',
      'rising_ap'=c(134,538,125),
      'rising_rp'=c(40,5,24),
      'falling_ap'=c(126,540,119),
      'falling_rp'=c(39,5,23),
      'time_on'=c(0,0,0),
      stringsAsFactors=FALSE
      ), 
    '2822'=data.frame(
      'meta_shape'='pattern_scan',
      'rising_ap'=c(224,115),
      'rising_rp'=c(66,109),
      'falling_ap'=c(191,101),
      'falling_rp'=c(59,106),
      'time_on'=c(0,0),
      stringsAsFactors=FALSE
    ), 
    '2880'=data.frame(
      'meta_shape'='pattern_scan',
      'rising_ap'=c(50),
      'rising_rp'=c(6),
      'falling_ap'=c(56),
      'falling_rp'=c(6),
      'time_on'=c(0),
      stringsAsFactors=FALSE
    ), 
    '2882'=data.frame(
      'meta_shape'='pattern_scan',
      'rising_ap'=c(152,2060),
      'rising_rp'=c(94,7),
      'falling_ap'=c(142,2060),
      'falling_rp'=c(93,8),
      'time_on'=c(0,0),
      stringsAsFactors=FALSE
    ),   
    '2894'=data.frame(
      'meta_shape'='pattern_scan',
      'rising_ap'=c(82),
      'rising_rp'=c(33),
      'falling_ap'=c(78),
      'falling_rp'=c(30),
      'time_on'=c(0),
      stringsAsFactors=FALSE
    ),
    '2893'=data.frame(
      'meta_shape'='pattern_scan',
      'rising_ap'=c(190),
      'rising_rp'=c(136),
      'falling_ap'=c(156),
      'falling_rp'=c(130),
      'time_on'=c(0),
      stringsAsFactors=FALSE
    ),
    '3700'=data.frame(
      'meta_shape'='pattern_scan',
      'rising_ap'=c(124),
      'rising_rp'=c(6),
      'falling_ap'=c(116),
      'falling_rp'=c(4.5),
      'time_on'=c(0),
      stringsAsFactors=FALSE
    ),
    '3795'=data.frame(
      'meta_shape'='pattern_scan',
      'rising_ap'=c(68),
      'rising_rp'=c(6),
      'falling_ap'=c(64),
      'falling_rp'=c(7),
      'time_on'=c(0),
      stringsAsFactors=FALSE
    ),
    '4033'=data.frame(
      'meta_shape'='pattern_scan',
      'rising_ap'=c(102),
      'rising_rp'=c(25),
      'falling_ap'=c(93),
      'falling_rp'=c(23),
      'time_on'=c(0),
      stringsAsFactors=FALSE
    ),
    '4465'=data.frame(
      'meta_shape'='pattern_scan',
      'rising_ap'=c(111,178),
      'rising_rp'=c(87,19),
      'falling_ap'=c(85,144),
      'falling_rp'=c(83,13),
      'time_on'=c(0),
      stringsAsFactors=FALSE
    )
  )
  
  file.list=list(
    '2813'='/home/hslee/test_data/meta_test_data/1442415600000-1443625199999.10002813_4637_3.csv.gz',
    '2822'='/home/hslee/test_data/meta_test_data/1442415600000-1443625199999.10002822_4646_3.csv.gz',
    '2880'='/home/hslee/test_data/meta_test_data/1442415600000-1443625199999.10002880_4704_3.csv.gz',
    '2882'='/home/hslee/test_data/meta_test_data/1442415600000-1443625199999.10002882_4706_3.csv.gz',
    #'2894'='/home/hslee/test_data/meta_test_data/1442415600000-1443625199999.10002894_4718_3.csv.gz',
    '2894'='/home/hslee/test_data/meta_test_data/1442242800000-1443365998999.10002894_4718_3.csv.gz',
    #'2894'='/home/hslee/test_data/meta_test_data/1442415600000-1443711599999.10002894_4718_3.csv.gz',
    '2893'='/home/hslee/test_data/meta_test_data/1442415600000-1443625199999.10002893_4717_3.csv.gz',
    #'2893'='/home/hslee/test_data/ClosedBeta/1442242800000-1443365998999.10002893_4717_3.csv',
    '3700'='/home/hslee/test_data/meta_test_data/1442415600000-1443625199999.10003700_5532_3.csv.gz',
    '3795'='/home/hslee/test_data/meta_test_data/1442415600000-1443625199999.10003795_5627_3.csv.gz',
    '4033'='/home/hslee/test_data/meta_test_data/1442415600000-1443625199999.10004033_5866_3.csv.gz',
    #'4465'='/home/hslee/test_data/meta_test_data/1442415600000-1443625199999.10004465_6299_3.csv.gz'
    '4465'='/home/hslee/test_data/meta_test_data/1442242800000-1443365998999.10004465_6299_3.csv.gz'
  )
  
  meta.json.list = list()
  meta.df.list = list()
  for( f in test.sites ) { #names(file.list) ){
    print(file.list[[f]])
    data.file <-  file.list[[f]]
  
    #  
    # import source file
    meta.list.json = app.detect.main(
      data.file = data.file, find.heavy=test.heavy, find.patternHigh=test.patternHigh,  find.pattern=test.pattern, find.pattern_extend=test.pattern_extend,
      find.rice=test.rice, find.standby=test.standby, find.cyclic=test.cyclic, find.ac=test.ac, 
      check.data=F, debug.mode=F)
      
    if( length(meta.list.json) > 0 ) meta.list = fromJSON(meta.list.json)
    if( length(meta.list)==0 ) meta.list = list() # == fromJSON(I(toJSON(''))) ) meta.list = list()
    #
    #
    
    #
    # 데이터 확인
    #
    meta.json.df = data.frame()
    for(k in names(meta.list)){
      meta <- meta.list[[k]]
      print(meta)
      meta_shape <- meta[['shape_type']] 
      if( meta_shape == 'pattern_scan' ){
        rrp = meta[['rising_edge']][['rp_med']]
        frp = meta[['falling_edge']][['EffRP_Drop.med']]
        if( is.null(rrp) | class(rrp) == 'character' ) rrp = 0
        if( is.null(frp) | class(frp) == 'character' ) frp = 0
        
        meta.json.df <- rbind(meta.json.df, 
                              data.frame(meta_shape=meta_shape,
                                         rising_ap=meta[['rising_edge']][['ap_med']],
                                         rising_rp=rrp,
                                         falling_ap=meta[['falling_edge']][['EffAP_Drop.med']], 
                                         falling_rp=frp,
                                         time_on=meta[['falling_edge']][['EffTimeOn.med']]))
      } else if( meta_shape == 'high_power' ){
        meta.json.df <- rbind(meta.json.df, 
                              data.frame(meta_shape=meta_shape,
                                         rising_ap=meta[['rising_edge']][['ap_height']],rising_rp=meta[['rising_edge']][['rp_height']],
                                         falling_ap=meta[['falling_edge']][['ap_height']],falling_rp=meta[['falling_edge']][['rp_height']],
                                         time_on=0))
      } else if( meta_shape == 'ricecooker_pattern_scan' ){
        meta.json.df <- rbind(meta.json.df, 
                              data.frame(meta_shape=meta_shape,
                                         rising_ap=meta[['rising_edge']][['ap_med']],rising_rp=meta[['rising_edge']][['rp_med']],
                                         falling_ap=0,falling_rp=0,
                                         time_on=meta[['falling_edge']][['cooking.num.count']]))
      } else if( meta_shape == 'cyclic_box' ){
        meta.json.df <- rbind(meta.json.df, 
                              data.frame(meta_shape=meta_shape,
                                         rising_ap=meta[['rising_edge']][['ap_height']],rising_rp=meta[['rising_edge']][['rp_height']],
                                         falling_ap=meta[['falling_edge']][['ap_height']],falling_rp=meta[['falling_edge']][['rp_height']],
                                         time_on=0))
      } 
    }
    meta.df.list[[f]] = meta.json.df
    meta.json.list[[f]] = meta.list.json
  }
  

  #
  # 주어진 답지와의 비교
  #
  correct.meta.list <- list()
  if( test.pattern_extend ){
    correct.meta.pattern.list<-lapply(correct.meta.pattern.list, function(x) {
      x$find <- T
      return(x)
    })
    correct.meta.list <- append(correct.meta.list,correct.meta.pattern.list)
  }
  if( test.cyclic ){
    correct.meta.cyclic.list<-lapply(correct.meta.cyclic.list, function(x) {
      x$find <- T
      return(x)
    })
    correct.meta.list <- append(correct.meta.list,correct.meta.cyclic.list)
  }
  
  
  
  for( k in names(correct.meta.list)){
    if( (k %in% test.sites) ){
      correct.meta.list[[k]]$find = F
    }
  }
  
  for( k in test.sites ) { #names(correct.meta.list)) {
    correct_meta = correct.meta.list[[k]]
    answer_meta = meta.df.list[[k]]
    if( length(answer_meta) < 1 ) 
      next
    
    for( i in seq(nrow(correct_meta))) {
      s = correct_meta[i,'meta_shape']
      rap = correct_meta[i,'rising_ap']
      rrp = correct_meta[i,'rising_rp']
      fap = correct_meta[i,'falling_ap']
      frp = correct_meta[i,'falling_rp']
      
      for( j in seq(nrow(answer_meta))) {
        c_s = answer_meta[j,'meta_shape']
        c_rap = abs(answer_meta[j,'rising_ap'])
        c_rrp = abs(answer_meta[j,'rising_rp'])
        c_fap = abs(answer_meta[j,'falling_ap'])
        c_frp = abs(answer_meta[j,'falling_rp'])
        
        if( s==c_s 
            & abs(rap-c_rap)/rap<0.15 
            & abs(rrp-c_rrp)/rrp<0.15 
            & abs(fap-c_fap)/fap<0.15 
            & abs(frp-c_frp)/frp<0.15 ) {
          correct.meta.list[[k]][i,'find'] = T
          break
        }
      }
    }
  }
  
  
  print('---------------------------- Output from Tests  --------------------------')
  print(meta.df.list)
  print('---------------------------- End of Output  --------------------------')
  
  #
  # print results
  print('---------------------------- Test Results --------------------------')
  print(paste('Test Site:',test.sites))
  print('Test failed list:')
  for( k in names(correct.meta.list)) {
    for( i in seq(nrow(correct.meta.list[[k]]))){
      if( correct.meta.list[[k]][i,'find'] == F ) {
        print(paste('site id:',k))
        print(correct.meta.list[[k]][i,])
      }
    }
  }
  print('---------------------------- End of Test  --------------------------')
  
  #
  # 메타 코드 테스트 끝!
  ###################################################################################

}



