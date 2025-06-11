##### Plots for the Publication #####

tany_integrated1 <- subset(x=tany_integrated, subset = Diet == "Chow" & Age == "6+ weeks") #Figures were plotted using Chow diet fed 6+ weeks old adult mice.

#Plot for the Fig 1c.

desired_order <- c("beta 2","beta 1","alpha 2","alpha 1") 
tany_integrated1$label <- factor(tany_integrated1$label, levels = desired_order)

features2 <- c("Homer1", "Dlg4", "Shank2", "Nlgn2", "Gphn")

plot <- DotPlot(
  tany_integrated1,
  features = features2,
  assay = "RNA",
  cols = c("lightgrey", "blue"),
  dot.min = 0,
  dot.scale = 6,
  group.by = "label",
  scale = TRUE,
  scale.by = "size"
) +  theme(plot.background = element_rect(fill = "white", colour = NA), 
        axis.title = element_text(size = 10),
        axis.text.x = element_text(angle = 45,face = 'italic', hjust = 1, vjust = 1),  # x-axis labels in italics
        axis.text.y = element_text(face = 'plain',hjust = 1, vjust = 1) # ,plot.margin = margin(t = 1, r = 5, b = 1, l = 1)
        ) + 
  xlab("Transcripts") + 
  ylab("Tanycyte Subtypes")  # y-axis label

ggsave(filename = "Fig_1c.pdf", height = 4, width = 5, plot = plot) 

#Plot for the Supplementary Fig 6a.
tany_integrated1 <- subset(x=tany_integrated, subset = Diet == "Chow" & Age == "6+ weeks") #Figures were plotted using Chow diet fed 6+ weeks old adult mice.

features1 <- c("Gria3", "Gabbr1", "Gabbr2", "Chrm1", "Trhr", "Kiss1r", "Crhr1", "Crhr2", "Oxtr", "Avpr1a", "Sstr2", "Adra1a","Adra2a","Adrb1","Adrb2","Adrb3","Drd1", "Drd2", "Drd5")

plot <- DotPlot(
  tany_integrated1,
  features = features1,
  assay = "RNA",
  cols = c("lightgrey", "blue"),
  dot.min = 0,
  dot.scale = 6,
  group.by = "label",
  scale = TRUE,
  scale.by = "radius",
  scale.min = NA,
  scale.max = 15
) +  theme(plot.background = element_rect(fill = "white", colour = NA), 
        axis.title = element_text(size = 14),
        axis.text.x = element_text(face = 'plain', hjust = 1, vjust = 1),  # x-axis labels in italics
        axis.text.y = element_text(face = 'italic',hjust = 1, vjust = 1), # y-axis labels 
        plot.margin = margin(t = 1, r = 5, b = 1, l = 1), legend.position = "bottom") + 
  ylab("Transcripts") + 
  xlab("Tanycyte Subtypes")  # y-axis label

# To flip the coordinates
plot <- plot + coord_flip()  # Flip the axes

ggsave(filename = "SuppFig6a.pdf", height = 8, width = 3, plot = plot)


#Plot for the Supplementary Fig 8a.
tany_integrated1 <- subset(x=tany_integrated, subset = Diet == "Chow" & Age == "6+ weeks") #Figures were plotted using Chow diet fed 6+ weeks old adult mice.

features <- c("Cd59a", "Lyz2", "Slc17a8", "Pygm", "Ephb1", "Vcan", "Crlf3", "Crym", "Frzb", "Pttg1", "Sprr1a", "A2m", "Scn7a", "Adm", "Col25a1", "Trhr", "Abhd11", "Trhde", "Cers6", "Dio2")

plot <- DotPlot(
  tany_integrated1,
  features = features,
  assay = "RNA",
  cols = c("lightgrey", "blue"),
  dot.min = 0,
  dot.scale = 6,
  group.by = "label",
  scale = TRUE,
  scale.by = "radius"
) +  theme(plot.background = element_rect(fill = "white", colour = NA), 
        axis.title = element_text(size = 14),
        axis.text.x = element_text(angle = 45,face = 'italic', hjust = 1, vjust = 1),  # x-axis labels in italics
        axis.text.y = element_text(face = 'plain',hjust = 1, vjust = 1), # y-axis labels 
        plot.margin = margin(t = 1, r = 5, b = 1, l = 1)) + 
  xlab("Transcripts") + 
  ylab("Tanycyte Subtypes")  # y-axis label
  
ggsave(filename = "SuppFig8a.pdf", height = 3.5, width =8, plot = plot) 