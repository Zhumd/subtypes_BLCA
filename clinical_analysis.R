    
data01 <- data.frame(x=rep(1,12), y=1:12,ID= paste0('P', seq(1,12)))  
dat02 <- data.frame(x=rep(4,12), y=1:12,Gender=c('Male','Female', 'Male','Male','Male','Male','Male','Female','Female','Male','Male','Male'))  
dat03 <- data.frame(x=rep(5,12), y=1:12,LVI=c('Positive','Negative', 'Positive','Negative','Positive','Negative','Positive','Positive','Positive','NA','NA','Negative'))  
dat04 <- data.frame(x=rep(6,12), y=1:12,TNM=c('T1N1M0','T1N0M0', 'T1N1M0','T1N0M0','T1N1M0','T1N1M0', 'T1N0M0','T1N0M0','T1N0M0','NA','NA','T1N0M0'))

data07 <- data.frame(x=rep(2,12), y=1:12,Age=c(57,69,65,60,74,62,62,64,68,58,58,28))  
dat06 <- data.frame(x=rep(7,12), y=1:12,Grade=c('High', 'High', 'High', 'High', 'High', 'High', 'High', 'High', 'High', 'High', 'High', 'Low'))
data08 <- data.frame(x=rep(3,12), y=1:12,Tissues=rep('BCa',12))
colnames(data01) <- c('x','y','label')
colnames(data07) <- c('x','y','label')
colnames(data08) <- c('x','y','label')
data01 <- rbind(data01,data07,data08)

dat05 <- data.frame(ymin=seq(1.5,11.5,by=2), ymax=seq(2.5,12.5,by=2))
dat09 <- data.frame(x=rep(8,12), y=1:12,Recurrence=c('No', 'Yes', 'Yes', 'No', 'No', 'No', 'Yes', 'No', 'Yes', 'No', 'Yes', 'No'))


x_labels<-c("Patients ID","Age","Tissues","Gender","LVI","TNM",#'Grade'
            "Grade","Recurrence")

ggplot()+
  scale_x_continuous(limits = c(0.5,8),
                     breaks = 1:8,
                     labels = x_labels,
                     position = "top")+
  geom_text(data=data01,aes(x=x,y=rev(y),label=label))+
  geom_point(data=dat02,aes(x=x,y=rev(y),color=`Gender`),
             size=5)+
  scale_color_manual(values = c("Male"="#d38e91",
                                "Female"="#1f639a"),
                     labels = c("Male" = "Male", "Female" = "Female")) +
  ggnewscale::new_scale_color()+
  geom_point(data=dat03,aes(x=x,y=rev(y),color=`LVI`),
             size=5)+
  scale_color_manual(values = c("Positive"="#f6d65b",
                                "Negative"="#3eada2",
                                "NA"="#eb553a"))+
  ggnewscale::new_scale_color()+
  geom_point(data=dat04,aes(x=x,y=rev(y),color=TNM),
             size=5,shape=15)+
  scale_color_manual(values = c("T1N1M0"="#f49b5e",
                                "T1N0M0"="#ecc6f6",
                                "NA"="#eb553a"))+
  ggnewscale::new_scale_color()+
  geom_point(data=dat06,aes(x=x,y=rev(y),color=Grade),
             size=5,shape=15)+
  scale_color_manual(values = c( "High"="#4848AB",
                                 "Low"="#89520e"))+
  ggnewscale::new_scale_color()+
  geom_point(data=dat09,aes(x=x,y=rev(y),color=Recurrence),
             size=5,shape=15)+
  scale_color_manual(values = c( "No"="#33FF57",
                                 "Yes"="#8B0000"))+
  ggnewscale::new_scale_color()+
  theme_bw()+
  theme(axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(hjust=0.5,vjust=0.5))+
  geom_rect(data=dat05,
            aes(xmin=-Inf,xmax=Inf,ymin=ymin,ymax=ymax),
            fill="gray",alpha=0.2)
ggsave('./results/P_LHBC/clinical2_BLCA.pdf',width = 8, height = 8)
