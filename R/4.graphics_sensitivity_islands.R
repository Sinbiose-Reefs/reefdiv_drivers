
## --------------------
# FIGURES

## chamar os pacotes
source("R/packages.R")
source("R/functions.R")


#### -------------------------
# Fig 2
# accumulation of FD relative to SR

load (here ("output", "complete_results_for_fig2_minus_islands.RData"))

# accumulation of diversity

## plotting
## Minimum Sample Size (MSS)

p1 <- ggplot (complete_results_for_fig2 [[1]], aes (x=EST.rich,y=FD, fill=Organism, 
                                              colour=Organism)) + 
  geom_point() +
  scale_color_manual(values=c("#f2a154", "#0e49b5")) +
  geom_smooth(method="lm", formula = y ~x) + 
  theme_classic() + xlab ("") + ylab ("Functional diversity")+
  theme (legend.position = "none") + 
  ggtitle ("Minimum Sample Size (MSS)")

## Precision-based Sample Size (PSSi)

p2 <- ggplot (complete_results_for_fig2 [[2]], aes (x=EST.rich,y=FD, fill=Organism, 
                                              colour=Organism)) + 
  geom_point() +
  scale_color_manual(values=c("#f2a154", "#0e49b5")) +
  geom_smooth(method="lm", formula = y ~x) + 
  theme_classic() + xlab ("Estimated Richness") + ylab ("")+
  theme (legend.position = "none") + 
  ggtitle ("Precision-based Sample Size (PSSi)")

## Model-based Sample Size (Lomolino's model)

p3 <- ggplot (complete_results_for_fig2 [[3]], aes (x=EST.rich,y=FD, fill=Organism, 
                                              colour=Organism)) + 
  geom_point() +
  scale_color_manual(values=c("#f2a154", "#0e49b5")) +
  geom_smooth(method="lm", formula = y ~x) + 
  theme_classic() + xlab ("") + ylab ("")+
  theme (legend.position = "none") + 
  ggtitle ("Lomolino's model")

pdf(file=here("output","vectorized","Fig2_minus_islands.pdf"),height=4,width=11)

grid.arrange(p1,p2,p3,
             ncol=9,nrow=6,
             layout_matrix = rbind (c(1,1,1,2,2,2,3,3,3),
                                    c(1,1,1,2,2,2,3,3,3),
                                    c(1,1,1,2,2,2,3,3,3),
                                    c(1,1,1,2,2,2,3,3,3),
                                    c(1,1,1,2,2,2,3,3,3),
                                    c(1,1,1,2,2,2,3,3,3)))
dev.off()

### -------------------------
#  FIG 3
# coefficients of drivers

load (here ("output", "complete_results_for_fig3_minus_islands.RData"))

# plot
dodge <- c(0.4,0.4)
pd <- position_dodge(dodge)
pdf_pt <- position_dodge(dodge)

a <- ggplot (complete_results[which(complete_results$Parameter != "Intercept"),], 
             
             aes  (y=Parameter, x=Estimate, 
                          colour=Algorithm, fill=Algorithm)) + 
  
  geom_errorbar(aes(xmin=lower,xmax=upper),width = 0.2,size=0.5,
                position=pd) + theme_classic() + 
  
  geom_point(position=(pdf_pt), 
             size=1.5)+ 
  
  geom_vline(xintercept = 0, linetype="dashed", 
             color = "gray50", size=0.1)+
  
  facet_wrap(~Organism+Index,scale="free",ncol=4) + 

  scale_color_manual(values=c("gray90", "gray60", "gray40")) + 
  
  xlab("Standardized effect size") + 
  
  ylab ("Parameter") + 
  
  #xlim(-0.5,0.5) +
  
  theme(axis.text.x = element_text(angle = 45,size=7))

a

ggsave (file=here("output","vectorized","Fig3_minus_islands.pdf"),width=11,height=6)


