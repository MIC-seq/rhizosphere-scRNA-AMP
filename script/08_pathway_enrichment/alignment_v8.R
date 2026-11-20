.libPaths(c('/public3/labmember/qianqh/soil_amp_analysis/14_metabolic_functions/software/Rlib',.libPaths()))
library(data.table);library(jsonlite)
b<-'/public3/labmember/qianqh/soil_amp_analysis/14_metabolic_functions/fg_MES_20261027/main_supp_v8';g<-fread(file.path(b,'panel_geometry_mm.tsv'));k<-72/25.4
pan<-lapply(1:nrow(g),function(i){r<-g[i];list(id=r$id,bbox_pt=c(15+r$x*k,910-582-112*k+r$y*k,15+(r$x+r$w)*k,910-582-112*k+(r$y+r$h)*k))})
rows<-lapply(1:3,function(i)list(id=paste0('crop_row_',i),panels=paste0('p',i,1:3)));cols<-lapply(1:3,function(j)list(id=paste0('crop_col_',j),panels=paste0('p',1:3,j)));cols[[4]]<-list(id='umap_col',panels=paste0('u',1:3))
write_json(list(schema_version=1,figure=list(width_pt=595.276,height_pt=910),panels=pan,row_groups=rows,column_groups=cols),file.path(b,'main_alignment.json'),auto_unbox=TRUE,pretty=TRUE)
