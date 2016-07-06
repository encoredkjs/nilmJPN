meta2cyclic_box <- function( meta.json, data, postprocessing = T, show.fig = F, filename = "meta2cyclic" ){
  
  if(meta.json[['shape_type']] != 'cyclic_box') stop("Meta is not correct : meta2cyclic_box")
  
  if( show.fig ) png( paste0(filename,"-%d.png"), width=3000)
  
  # step 1
  box.shape <- find.box.shape(meta.json, data, show.fig); pre.box.shape <- box.shape
  
  # step 2
  if( postprocessing ) box.shape <- post.processing(meta.json, pre.box.shape, show.fig)

  if( show.fig ) dev.off()
  return( box.shape )
}

meta2heavy_load <- function( meta.json, data, show.fig = F, filename = "meta2heavy" ){
  
  if(meta.json[['shape_type']] != 'high_power') stop("Meta is not correct : meta2heavy_load")
  
  if( show.fig ) png( paste0(filename,"-%d.png"), width=3000)
  
  # step 1
  heavy.load <- find.heavy.load( meta.json=meta.json, data=data, debug.mode=show.fig)
  
  if( show.fig ) dev.off()
  return( heavy.load )
}