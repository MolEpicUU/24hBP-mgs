rm(list=ls())

library(rio)
library(ggplot2)
library(ggpubr)
library(ggforce)
library(ggnewscale)
library(data.table)

#functions
scale_override <- function(which, scale) {
  if(!is.numeric(which) || (length(which) != 1) || (which %% 1 != 0)) {
    stop("which must be an integer of length 1")
  }

  if(is.null(scale$aesthetics) || !any(c("x", "y") %in% scale$aesthetics)) {
    stop("scale must be an x or y position scale")
  }

  structure(list(which = which, scale = scale), class = "scale_override")
}

CustomFacetWrap <- ggproto(
  "CustomFacetWrap", FacetWrap,
  init_scales = function(self, layout, x_scale = NULL, y_scale = NULL, params) {
    # make the initial x, y scales list
    scales <- ggproto_parent(FacetWrap, self)$init_scales(layout, x_scale, y_scale, params)

    if(is.null(params$scale_overrides)) return(scales)

    max_scale_x <- length(scales$x)
    max_scale_y <- length(scales$y)

    # ... do some modification of the scales$x and scales$y here based on params$scale_overrides
    for(scale_override in params$scale_overrides) {
      which <- scale_override$which
      scale <- scale_override$scale

      if("x" %in% scale$aesthetics) {
        if(!is.null(scales$x)) {
          if(which < 0 || which > max_scale_x) stop("Invalid index of x scale: ", which)
          scales$x[[which]] <- scale$clone()
        }
      } else if("y" %in% scale$aesthetics) {
        if(!is.null(scales$y)) {
          if(which < 0 || which > max_scale_y) stop("Invalid index of y scale: ", which)
          scales$y[[which]] <- scale$clone()
        }
      } else {
        stop("Invalid scale")
      }
    }

    # return scales
    scales
  }
)

facet_wrap_custom <- function(..., scale_overrides = NULL) {
  # take advantage of the sanitizing that happens in facet_wrap
  facet_super <- facet_wrap(...)

  # sanitize scale overrides
  if(inherits(scale_overrides, "scale_override")) {
    scale_overrides <- list(scale_overrides)
  } else if(!is.list(scale_overrides) ||
            !all(vapply(scale_overrides, inherits, "scale_override", FUN.VALUE = logical(1)))) {
    stop("scale_overrides must be a scale_override object or a list of scale_override objects")
  }

  facet_super$params$scale_overrides <- scale_overrides

  ggproto(NULL, CustomFacetWrap,
    shrink = facet_super$shrink,
    params = facet_super$params
  )
}


info=import("raw/scapis_metagenomics_mgs_annotations_v1.0.tsv")
names(info)[1]="var.x"

#import results(main results)
sbp=import("results/res_fastlm.24SBP_multiBMI_model.tsv")
dbp=import("results/res_fastlm.24DBP_multiBMI_model.tsv")


sbp=sbp[which(sbp$q.value<0.05),]
dbp=dbp[which(dbp$q.value<0.05),]

sbp$var.y="ABPM_SBP"
dbp$var.y="ABPM_DBP"

sbp=sbp[,which(names(sbp)%in%grep(".infl",names(sbp),value=T)==F)]
dbp=dbp[,which(names(dbp)%in%grep(".infl",names(dbp),value=T)==F)]

#import results(office - main popualtion)
sbp_val=import("results/res_fastlm.SBP_multiBMI_model_validation.tsv")
dbp_val=import("results/res_fastlm.DBP_multiBMI_model_validation.tsv")

sbp_val=sbp_val[which(sbp_val$var.x %in% sbp$var.x),]
dbp_val=dbp_val[which(dbp_val$var.x %in% dbp$var.x),]

#ABPM model2
sbp2=import("results/res_fastlm.24SBP_multi_model.tsv")
dbp2=import("results/res_fastlm.24DBP_multi_model.tsv")

sbp2=sbp2[which(sbp2$var.x %in% sbp$var.x),]
dbp2=dbp2[which(dbp2$var.x %in% dbp$var.x),]

sbp_com=rbind(sbp,sbp2[,grep("infl",names(sbp2),value=T,invert=T)],sbp_val)
dbp_com=rbind(dbp,dbp2[,grep("infl",names(dbp2),value=T,invert=T)],dbp_val)

####################################
#forestplot example 1
sbp_com$var.y<-ifelse(sbp_com$var.y=="bbps","Office BP (non-ABPM cohort)",
ifelse(sbp_com$var.y=="ABPM_SBP","24-hours BP Model 1 (ABPM cohort)","24-hours BP Model 2 (ABPM cohort)"))
sbp_com$var.y <- factor(sbp_com$var.y,levels=c("24-hours BP Model 1 (ABPM cohort)","24-hours BP Model 2 (ABPM cohort)","Office BP (non-ABPM cohort)"))
cols <- c("24-hours BP Model 1 (ABPM cohort)","24-hours BP Model 2 (ABPM cohort)","Office BP (non-ABPM cohort)")
#model <- rev(cols)[sbp_com$var.y]
colnames(sbp_com)[colnames(sbp_com) == "var.y"] ="Model"

##############################
dbp_com$var.y<-ifelse(dbp_com$var.y=="bbpd","Office BP (non-ABPM cohort)",
ifelse(dbp_com$var.y=="ABPM_DBP","24-hours BP Model 1 (ABPM cohort)","24-hours BP Model 2 (ABPM cohort)"))
dbp_com$var.y <- factor(dbp_com$var.y,levels=c("24-hours BP Model 1 (ABPM cohort)","24-hours BP Model 2 (ABPM cohort)","Office BP (non-ABPM cohort)"))
cols <- c("24-hours BP Model 1 (ABPM cohort)","24-hours BP Model 2 (ABPM cohort)","Office BP (non-ABPM cohort)")
#model <- rev(cols)[sbp_com$var.y]
colnames(dbp_com)[colnames(dbp_com) == "var.y"] ="Model"

dbp_com$var="DBP"
sbp_com$var="SBP"


xxx=data.frame(rbind(dbp_com,sbp_com),stringsAsFactors=F)
xxx$var<-ifelse(xxx$var=="DBP","DBP"," SBP")
xxx$var <- factor(xxx$var, ordered = T, levels=c(" SBP", "DBP"))

xxx=merge(xxx,info[,c("var.x","maintax")],by="var.x")
xxx$var.x=paste(xxx$maintax, "\n(",xxx$var.x,")",sep="")


t1=c("Streptococcus sp001556435\n(hMGS.01022)", "Intestinimonas massiliensis\n(hMGS.01559)")

ann_text<-data.frame(
        var.x=t1,
        estimate=rep(15,length(t1)),
        var=rep("DBP",length(t1)),
        Label=rep("N.S.",length(t1)),
        Model=rep("Office (ABPM cohort)",length(t1)))

#x=c("[Ruminococcus] torques (HG3A.0034)","Oscillospiraceae sp. (HG3A.0130)","Gemmiger formicilis (HG3A.0027)","Faecalibacterium prausnitzii (HG3A.0025)","Eubacteriales sp. (HG3A.0069)")
#ann_text[ann_text$var.x%in%x,"var"]=" SBP"


xxx$var.x=as.factor(xxx$var.x)

x1=xxx[which(xxx$Model%in%grep("non-ABPM",xxx$Model,value=T)==F),]
x2=xxx[which(xxx$Model%in%grep("non-ABPM",xxx$Model,value=T)==T),]


f3=ggplot(data=xxx) +
geom_point(data=xxx,aes(reorder(var.x, estimate), y = -500, color=Model),size = 2,position=position_dodge2(0.9),show.legend = F)+
 geom_rect(aes(xmin = 1-0.5 , xmax =1+0.5 , ymin = -Inf, ymax = Inf), fill="grey90",color=NA,alpha = 0.1)+
 geom_rect(aes(xmin = 3-0.5 , xmax =3+0.5 , ymin = -Inf, ymax = Inf), fill="grey90",color=NA,alpha = 0.1)+
# geom_rect(aes(xmin = 53-0.5 , xmax =53+0.5 , ymin = -Inf, ymax = Inf), fill="grey90",color=NA,alpha = 0.1)+

   geom_hline(yintercept=0, lty=2) +
geom_text(data=ann_text,aes(x=reorder(var.x, estimate),y=-0.3),label=ann_text$Label,size=7,color="black",fontface ="bold")+
  theme_bw()+
theme(
panel.background = element_rect(fill = "white"),
 panel.grid.major.x = element_line(color = 'gray80', linetype = 3),
        panel.grid.major.y = element_blank(),
 #panel.grid.minor.x = element_line(color = 'gray80', linetype = 3),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=20),
        axis.text.y = element_text(size=15, face="italic"), #element_blank(),# this makes the horizontal Strep labels disappear. If you want to allow them back in remove the element_blank() and add this element_text(size=9)
        axis.text.x = element_text(size=20),
        axis.ticks.length=unit(0,"cm"),
legend.title=element_text(size=20),

        legend.position="top",
        legend.key.size = unit(0.5,"cm"),
        legend.text = element_text(size=20),

        strip.placement = "inside",  # outside
        strip.text =  element_text(size=20, face="bold"), #element_text(size=9, face = 'bold'),  # replace this with element_blank() if you turn axis.text.y into element_text

        panel.border = element_rect(color = "gray95", fill = NA, size = 1),
        panel.spacing = unit(0, "lines"), # space between facets

        strip.background =element_blank(),

        panel.grid.minor.y = element_blank(),
panel.spacing.x = unit(2, "lines")
)+
scale_color_manual(name="b",
values = c("#0077BB",  "#33BBEE",  "red4"),
labels=c("a",  "b",  "c"),
guides(color="none"))+
  new_scale_color()+
 geom_point(data=x1,aes(x=reorder(var.x, estimate), y = -500, color=Model),size = 2,position=position_dodge2(1))+
scale_color_manual(name="ABPM cohort:",
values = c("red4","#0077BB"),
labels=c("24-hours BP Model 1","24-hours BP Model 2"))+
new_scale_color()+
geom_point(data=x2,aes(x=reorder(var.x, estimate), y = -500, color=Model),size = 2,position=position_dodge2(0.9))+
scale_color_manual(name="non-ABPM cohort:",
values = "#33BBEE",
labels="Office BP")+

  coord_flip()+

facet_wrap_custom(~var,ncol=2,scales = "free_x", scale_overrides = list(
    scale_override(1, scale_y_continuous(limits=c(round(min(xxx$lower[xxx$var%in%" SBP"])), round(max(xxx$upper[xxx$var%in%" SBP"]))+0.25))),
    scale_override(2, scale_y_continuous(limits = c(round(min(xxx$lower[xxx$var%in%"DBP"]))-0.4, round(max(xxx$upper[xxx$var%in%"DBP"]))+0.1)))))+

xlab("MGS") + ylab("Estimate")+
#ylim(round(min(xxx$lower))-1,round(max(xxx$upper))+1)+
 new_scale_color()+
geom_point(data=xxx,aes(reorder(var.x,estimate), y = estimate, color=Model),size = 1.8,position=position_dodge2(0.75),show.legend = F)+
geom_errorbar(data=xxx,aes(x=reorder(var.x,estimate),y=estimate,ymax = upper, ymin = lower,color=Model),position=position_dodge2(1),width=0.75,show.legend = F)+
scale_color_manual(name="a",
values = c("#0077BB","#33BBEE","red4" ),
labels=c("x","y","z"),
guides(color="none"))+
 guides(color = guide_legend(reverse=T))

ggsave(f3,file="/proj/nobackup/sens2019512/wharf/ssayols/ssayols-sens2019512/test.pdf",width=17)



ggsave(f3,file="results/figures/figure2.pdf",width=17)


