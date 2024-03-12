################################################################################
####              Figures homomer duplication                               ####
################################################################################

# Load libraries
library(tidyverse)
library(magrittr)
library(ggplot2)
library(cowplot)
library(Cairo)
library(ggpubr)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(viridis)
library(agricolae)

theme_set(theme_cowplot() +
            theme(panel.background = element_rect(fill = 'white'),
                  plot.background = element_rect(fill = 'white'),
                  axis.text = element_text(size = 10),
                  axis.title = element_text(size = 12, face = 'bold'),
                  strip.text = element_text(size = 12, face = 'bold'),
                  legend.text = element_text(size = 10),
                  legend.title = element_text(size = 12),
                  axis.line = element_blank(),
                  strip.background = element_blank(),
                  panel.grid.major = element_line(colour="#8C8C8C", linetype = 'dashed'),
                  panel.grid.minor = element_blank(),
                  panel.border = element_rect(colour = 'black', linewidth = 1) 
            )
)

panel_label_size = 14

### Set the path to the main directory
setwd('/path/to/Homomer_duplication_2024/')

# A function to draw heatmaps
draw_CHeatmap <- function(in_heatmap){
  return(grid.grabExpr(draw(in_heatmap)))
}

## Variables for heatmaps ##
row_title_size = 8
column_title_size = 8
row_name_size = 7
column_name_size = 7
legend_title_size = 7
legend_label_size = 7
heatmap_width = 2.5
heatmap_height = 2.5
legend_height = 1.5
legend_width = 0.15
ht_opt$TITLE_PADDING = unit(c(0.5, 0.5), "points")

#### Figure 1 ####

#### Figure 1A ####
## Load the diagram for the simulations

fig_1A <- ggdraw() + 
  draw_image('Figures/Diagrams/Diagram1A.png')

fig_1A

#### Figure 1B ####

fig_1B <- ggdraw() +
  draw_image('Figures/Diagrams/Diagram1B.png')
fig_1B

#### Figure 1C ####

#### Effects of folding energy ####

# Load the data on protein stability
data_stab <- read_delim('Results_solution_space/solutionSpace_subunitStab.tsv', delim = '\t')

# Add columns for the percentages of complexes, without considering monomers
data_stab %<>% mutate(pct_AA = 100 * cAA / (cAA + cBB + cAB),
                      pct_AB = 100 * cAB / (cAA + cBB + cAB),
                      pct_BB = 100 * cBB / (cAA + cBB + cAB))

data_stab_wide <- data_stab %>% select(new_stab_A, new_stab_B, pct_AB) %>%
  arrange(desc(new_stab_A)) %>%
  pivot_wider(names_from = new_stab_B, values_from = pct_AB)

# Move the first column to the index
first_column <- data_stab_wide$new_stab_A
data_stab_wide <- as.matrix(data_stab_wide[, 2:ncol(data_stab_wide)])
rownames(data_stab_wide) <- first_column

#### Draw the heatmap ####
# Define column labels
column_labels = seq(from = -10, to = 9.99, by = 0.01)
new_column_labels = column_labels
for(i in 1:length(column_labels)){
  
  cur_val = column_labels[i]
  new_val = ifelse(mod(cur_val, 5) == 0 && cur_val != -10,
                   toString(cur_val), '')
  
  new_column_labels[i] = new_val
  
}

# Define row labels
row_labels = rev(seq(from = -10, to = 9.99, by = 0.01))
new_row_labels = row_labels
for(i in 1:length(row_labels)){
  
  cur_val = row_labels[i]
  new_val = ifelse(mod(cur_val, 5) == 0 && cur_val != -10,
                   toString(cur_val), '')
  
  new_row_labels[i] = new_val
  
}


p_stab_HET <- Heatmap(
  data_stab_wide[seq(from = 4, to = 2000, by = 4), seq(from = 4, to = 2000, by = 4)],
  cluster_columns = F, cluster_rows = F, 
  col = colorRamp2(
    breaks = c(0, 20, 40, 60, 80, 100),
    colors = viridis(6)),
  show_column_names = T, row_names_side = 'left',
  width=unit(heatmap_width, 'cm'), height = unit(heatmap_height, 'cm'),
  border = T,
  row_title = expression(paste(bold('\u0394'),
                               bold(G['fold,A']), 
                               sep = '')),
  column_title = expression(paste(bold('\u0394'),
                                  bold(G['fold,B']), sep = '')),
  row_title_gp = gpar(fontsize=row_title_size, fontface = 'bold'),
  column_title_gp = gpar(fontsize=column_title_size, fontface = 'bold'),
  column_title_side = 'bottom',
  row_names_rot = 90, 
  row_names_centered = T,
  row_names_gp = gpar(fontsize=row_name_size, fontface = 'bold', align = 'left'),
  column_names_gp = gpar(fontsize=column_name_size, fontface='bold', align = 'left'),
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(
    at = c(0, 20, 40, 60, 80, 100),
    # title = "Percentage of HET AB (%)", 
    title = '',
    title_gp = gpar(fontsize = legend_title_size),
    legend_height = unit(legend_height, "cm"),
    grid_width = unit(legend_width, "cm"),
    border='black',
    lwd=1.7,
    labels_gp = gpar(fontsize = legend_label_size),
    title_position = "leftcenter-rot"
  ), 
  row_labels = new_row_labels[seq(from  = 4, to = 2000, by = 4)],
  column_labels = new_column_labels[seq(from = 1, to = 2000, by = 4)], 
  rect_gp = gpar(col = 'transparent',
                 lwd = 0
  ), use_raster = TRUE
)
p_stab_HET

#### Figure 1D ####
#### Repeat for HM BB ####
data_stab_wide <- data_stab %>% select(new_stab_A, new_stab_B, pct_BB) %>%
  arrange(desc(new_stab_A)) %>%
  pivot_wider(names_from = new_stab_B, values_from = pct_BB)

# Move the first column to the index
first_column <- data_stab_wide$new_stab_A
data_stab_wide <- as.matrix(data_stab_wide[, 2:ncol(data_stab_wide)])
rownames(data_stab_wide) <- first_column

#### Draw the heatmap ####
# Define column labels
column_labels = seq(from = -10, to = 9.99, by = 0.01)
new_column_labels = column_labels
for(i in 1:length(column_labels)){
  
  cur_val = column_labels[i]
  new_val = ifelse(mod(cur_val, 5) == 0 && cur_val != -10,
                   toString(cur_val), '')
  
  new_column_labels[i] = new_val
  
}

# Define row labels
row_labels = rev(seq(from = -10, to = 9.99, by = 0.01))
new_row_labels = row_labels
for(i in 1:length(row_labels)){
  
  cur_val = row_labels[i]
  new_val = ifelse(mod(cur_val, 5) == 0 && cur_val != -10, 
                   toString(cur_val), '')
  
  new_row_labels[i] = new_val
  
}

p_stab_HM_BB <- Heatmap(
  data_stab_wide[seq(from = 4, to = 2000, by = 4), seq(from = 4, to = 2000, by = 4)],
  cluster_columns = F, cluster_rows = F, 
  col = colorRamp2(
    breaks = c(0, 20, 40, 60, 80, 100),
    colors = viridis(6)),
  show_column_names = T, row_names_side = 'left',
  width=unit(heatmap_width, 'cm'), height = unit(heatmap_height, 'cm'),
  border = T,
  row_title = expression(paste(bold('\u0394'),
                               bold(G['fold,A']), 
                               sep = '')),
  column_title = expression(paste(bold('\u0394'),
                                  bold(G['fold,B']), 
                                  sep = '')),
  row_title_gp = gpar(fontsize=row_title_size, fontface = 'bold'),
  column_title_gp = gpar(fontsize=column_title_size, fontface = 'bold'),
  column_title_side = 'bottom',
  row_names_rot = 90, 
  row_names_centered = T,
  row_names_gp = gpar(fontsize=row_name_size, fontface = 'bold', align = 'left'),
  column_names_gp = gpar(fontsize=column_name_size, fontface='bold', align = 'left'),
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(
    at = c(0, 20, 40, 60, 80, 100),
    # title = "Percentage of HM BB (%)", 
    title = '', 
    title_gp = gpar(fontsize = legend_title_size),
    legend_height = unit(legend_height, "cm"),
    grid_width = unit(legend_width, "cm"),
    border='black',
    lwd=1.7,
    labels_gp = gpar(fontsize = legend_label_size),
    title_position = "leftcenter-rot"
  ), 
  row_labels = new_row_labels[seq(from = 4, to = 2000, by = 4)],
  column_labels = new_column_labels[seq(from = 1, to = 2000, by = 4)], 
  rect_gp = gpar(col = 'transparent',
                 lwd = 0
  ), use_raster = TRUE
)
p_stab_HM_BB

#### Figure 1E ####
#### A figure for total complexes ####

data_stab_wide <- data_stab %>% 
  mutate(total_complexes = cAA + cAB + cBB, 
         total_activity = 0.1*cA + 0.1*cB + 1*cAA + 1*cBB + 1*cAB) %>%
  select(new_stab_A, new_stab_B, total_activity) %>%
  arrange(desc(new_stab_A)) %>%
  pivot_wider(names_from = new_stab_B, values_from = total_activity)

# Move the first column to the index
first_column <- data_stab_wide$new_stab_A
data_stab_wide <- as.matrix(data_stab_wide[, 2:ncol(data_stab_wide)])
rownames(data_stab_wide) <- first_column

#### Draw the heatmap ####
# Define column labels
column_labels = seq(from = -10, to = 9.99, by = 0.01)
new_column_labels = column_labels
for(i in 1:length(column_labels)){
  
  cur_val = column_labels[i]
  new_val = ifelse(mod(cur_val, 5) == 0 && cur_val != -10,
                   toString(cur_val), '')
  
  new_column_labels[i] = new_val
  
}

# Define row labels
row_labels = rev(seq(from = -10, to = 9.99, by = 0.01))
new_row_labels = row_labels
for(i in 1:length(row_labels)){
  
  cur_val = row_labels[i]
  new_val = ifelse(mod(cur_val, 5) == 0 && cur_val != -10,
                   toString(cur_val), '')
  
  new_row_labels[i] = new_val
  
}

p_stab_total_complexes <- Heatmap(
  data_stab_wide[seq(from = 4, to = 2000, by = 4), seq(from = 4, to = 2000, by = 4)],
  cluster_columns = F, cluster_rows = F, 
  col = colorRamp2(
    breaks = c(0, 20, 40, 60, 80, 100),
    colors = magma(6)),
  show_column_names = T, row_names_side = 'left',
  width=unit(heatmap_width, 'cm'), height = unit(heatmap_height, 'cm'),
  border = T,
  row_title = expression(paste(bold('\u0394'),
                               bold(G['fold,A']), 
                               sep = '')),
  column_title = expression(paste(bold('\u0394'),
                                  bold(G['fold,B']), 
                                  sep = '')),
  row_title_gp = gpar(fontsize=row_title_size, fontface = 'bold'),
  column_title_gp = gpar(fontsize=column_title_size, fontface = 'bold'),
  column_title_side = 'bottom',
  row_names_rot = 90, 
  row_names_centered = T,
  row_names_gp = gpar(fontsize=row_name_size, fontface = 'bold', align = 'left'),
  column_names_gp = gpar(fontsize=column_name_size, fontface='bold', align = 'left'),
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(
    at = c(0, 20, 40, 60, 80, 100),
    # title = 'Total activity',
    title = '',
    title_gp = gpar(fontsize = legend_title_size),
    legend_height = unit(legend_height, "cm"),
    grid_width = unit(legend_width, "cm"),
    border='black',
    lwd=1.7,
    labels_gp = gpar(fontsize = legend_label_size),
    title_position = "leftcenter-rot"
  ), 
  row_labels = new_row_labels[seq(from = 4, to = 2000, by = 4)],
  column_labels = new_column_labels[seq(from = 1, to = 2000, by = 4)], 
  rect_gp = gpar(col = 'transparent',
                 lwd = 0
  ), use_raster = TRUE
)
p_stab_total_complexes

#### Figure 1F ####

alpha = 80
beta = 0.5

data_stab_wide <- data_stab %>% mutate(total_complexes = cAA + cAB + cBB) %>%
  mutate(total_activity = cAA + cAB + cBB + (0.1*cA) + (0.1*cB)) %>%
  mutate(logw = round(log2(beta^(log2(total_activity / alpha)^2)), 5)) %>%
  select(new_stab_A, new_stab_B, logw) %>%
  arrange(desc(new_stab_A)) %>%
  pivot_wider(names_from = new_stab_B, values_from = logw)

# Move the first column to the index
first_column <- data_stab_wide$new_stab_A
data_stab_wide <- as.matrix(data_stab_wide[, 2:ncol(data_stab_wide)])
rownames(data_stab_wide) <- first_column

#### Draw the heatmap ####
# Define column labels
column_labels = seq(from = -10, to = 9.99, by = 0.01)
new_column_labels = column_labels
for(i in 1:length(column_labels)){
  
  cur_val = column_labels[i]
  new_val = ifelse(mod(cur_val, 5) == 0 && cur_val != -10,
                   toString(cur_val), '')
  
  new_column_labels[i] = new_val
  
}

# Define row labels
row_labels = rev(seq(from = -10, to = 9.99, by = 0.01))
new_row_labels = row_labels
for(i in 1:length(row_labels)){
  
  cur_val = row_labels[i]
  new_val = ifelse(mod(cur_val, 5) == 0 && cur_val != -10,
                   toString(cur_val), '')
  
  new_row_labels[i] = new_val
  
}

p_stab_logw <- Heatmap(
  data_stab_wide[seq(from = 4, to = 2000, by = 4), seq(from = 4, to = 2000, by = 4)],
  cluster_columns = F, cluster_rows = F, 
  col = colorRamp2(
    breaks = c(-1, -0.75, -0.5, -0.25, 0),
    colors = cividis(5)
  ),
  show_column_names = T, row_names_side = 'left',
  width=unit(heatmap_width, 'cm'), height = unit(heatmap_height, 'cm'),
  border = T,
  row_title = expression(paste(bold('\u0394'),
                               bold(G['fold,A']), 
                               sep = '')),
  column_title = expression(paste(bold('\u0394'),
                                  bold(G['fold,B']), 
                                  sep = '')),
  row_title_gp = gpar(fontsize=row_title_size, fontface = 'bold'),
  column_title_gp = gpar(fontsize=column_title_size, fontface = 'bold'),
  column_title_side = 'bottom',
  row_names_rot = 90, 
  row_names_centered = T,
  row_names_gp = gpar(fontsize=row_name_size, fontface = 'bold', align = 'left'),
  column_names_gp = gpar(fontsize=column_name_size, fontface='bold', align = 'left'),
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(
    # title = 'log2(Fitness)',
    title = '',
    title_gp = gpar(fontsize = legend_title_size),
    legend_height = unit(legend_height, "cm"),
    grid_width = unit(legend_width, "cm"),
    border='black',
    lwd=1.7,
    labels_gp = gpar(fontsize = legend_label_size),
    title_position = "leftcenter-rot"
  ), 
  row_labels = new_row_labels[seq(from = 4, to = 2000, by = 4)],
  column_labels = new_column_labels[seq(from = 1, to = 2000, by = 4)], 
  rect_gp = gpar(col = 'transparent',
                 lwd = 0
  ), use_raster = TRUE
)
p_stab_logw

#### Repeat the whole process when modifying binding affinity ####

#### Figure 1G ####

# Load the data on binding affinity
data_binding <- read_delim('Results_solution_space/solutionSpace_bindingEnergy.tsv', delim = '\t')

# Add columns for the percentages of complexes, without considering monomers
data_binding %<>% mutate(pct_AA = 100 * cAA / (cAA + cBB + cAB),
                         pct_AB = 100 * cAB / (cAA + cBB + cAB),
                         pct_BB = 100 * cBB / (cAA + cBB + cAB))

data_binding_wide <- data_binding %>% select(new_binding_AB, new_binding_BB, pct_AB) %>%
  arrange(desc(new_binding_AB)) %>%
  ## Limit to -15 to -5 to see the most interesting region
  filter(new_binding_AB <= -5, new_binding_AB >= -15,
         new_binding_BB <= -5, new_binding_BB >= -15) %>%
  pivot_wider(names_from = new_binding_BB, values_from = pct_AB)

# Move the first column to the index
first_column <- data_binding_wide$new_binding_AB
data_binding_wide <- as.matrix(data_binding_wide[, 2:ncol(data_binding_wide)])
rownames(data_binding_wide) <- first_column

#### Draw the heatmap ####
# Define column labels
column_labels = seq(from = -15, to = -5.01, by = 0.01)
new_column_labels = column_labels
for(i in 1:length(column_labels)){
  
  cur_val = column_labels[i]
  new_val = ifelse(mod(cur_val, 2) == 0 && cur_val != -15,
                   toString(cur_val), '')
  
  new_column_labels[i] = new_val
  
}

# Define row labels
row_labels = rev(seq(from = -15, to = -5.01, by = 0.01))
new_row_labels = row_labels
for(i in 1:length(row_labels)){
  
  cur_val = row_labels[i]
  new_val = ifelse(mod(cur_val, 2) == 0 && cur_val != -15,
                   toString(cur_val), '')
  
  new_row_labels[i] = new_val
  
}
p_binding_HET <- Heatmap(
  data_binding_wide,
  cluster_columns = F, cluster_rows = F, 
  col = colorRamp2(
    breaks = c(0, 20, 40, 60, 80, 100),
    colors = viridis(6)),
  show_column_names = T, row_names_side = 'left',
  width=unit(heatmap_width, 'cm'), height = unit(heatmap_height, 'cm'),
  border = T,
  row_title = expression(paste(bold('\u0394'),
                               bold(G['bind, HET AB']), 
                               sep = '')),
  column_title = expression(paste(bold('\u0394'),
                                  bold(G['bind,HM BB']), 
                                  sep = '')),
  row_title_gp = gpar(fontsize=row_title_size, fontface = 'bold'),
  column_title_gp = gpar(fontsize=column_title_size, fontface = 'bold'),
  column_title_side = 'bottom',
  row_names_rot = 90, 
  row_names_centered = T,
  row_names_gp = gpar(fontsize=row_name_size, fontface = 'bold', align = 'left'),
  column_names_gp = gpar(fontsize=column_name_size, fontface='bold', align = 'left'),
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(
    at = c(0, 20, 40, 60, 80, 100),
    # title = "Percentage of HET AB (%)", 
    title = '',
    title_gp = gpar(fontsize = legend_title_size),
    legend_height = unit(legend_height, "cm"),
    grid_width = unit(legend_width, "cm"),
    border='black',
    lwd=1.7,
    labels_gp = gpar(fontsize = legend_label_size),
    title_position = "leftcenter-rot"
  ), 
  row_labels = new_row_labels,
  column_labels = new_column_labels, 
  rect_gp = gpar(col = 'transparent',
                 lwd = 0
  ), use_raster = TRUE, 
  name = 'p_binding_HET'
)
p_binding_HET

#### Figure 1H ####
#### Percentage of HM BB ####

data_binding_wide <- data_binding %>% select(new_binding_AB, new_binding_BB, pct_BB) %>%
  arrange(desc(new_binding_AB)) %>%
  ## Limit to -15 to -5 to see the most interesting region
  filter(new_binding_AB <= -5, new_binding_AB >= -15,
         new_binding_BB <= -5, new_binding_BB >= -15) %>%
  pivot_wider(names_from = new_binding_BB, values_from = pct_BB)

# Move the first column to the index
first_column <- data_binding_wide$new_binding_AB
data_binding_wide <- as.matrix(data_binding_wide[, 2:ncol(data_binding_wide)])
rownames(data_binding_wide) <- first_column

#### Draw the heatmap ####
# Define column labels
column_labels = seq(from = -15, to = -5.01, by = 0.01)
new_column_labels = column_labels
for(i in 1:length(column_labels)){
  
  cur_val = column_labels[i]
  new_val = ifelse(mod(cur_val, 2) == 0 && cur_val != -15,
                   toString(cur_val), '')
  
  new_column_labels[i] = new_val
  
}

# Define row labels
row_labels = rev(seq(from = -15, to = -5.01, by = 0.01))
new_row_labels = row_labels
for(i in 1:length(row_labels)){
  
  cur_val = row_labels[i]
  new_val = ifelse(mod(cur_val, 2) == 0 && cur_val != -15,
                   toString(cur_val), '')
  
  new_row_labels[i] = new_val
  
}

p_binding_HM_BB <- Heatmap(
  data_binding_wide,
  cluster_columns = F, cluster_rows = F, 
  col = colorRamp2(
    breaks = c(0, 20, 40, 60, 80, 100),
    colors = viridis(6)),
  show_column_names = T, row_names_side = 'left',
  width=unit(heatmap_width, 'cm'), height = unit(heatmap_height, 'cm'),
  border = T,
  row_title = expression(paste(bold('\u0394'),
                               bold(G['bind,HET AB']), 
                               sep = '')),
  column_title = expression(paste(bold('\u0394'),
                                  bold(G['bind,HM BB']), 
                                  sep = '')),
  row_title_gp = gpar(fontsize=row_title_size, fontface = 'bold'),
  column_title_gp = gpar(fontsize=column_title_size, fontface = 'bold'),
  column_title_side = 'bottom',
  row_names_rot = 90, 
  row_names_centered = T,
  row_names_gp = gpar(fontsize=row_name_size, fontface = 'bold', align = 'left'),
  column_names_gp = gpar(fontsize=column_name_size, fontface='bold', align = 'left'),
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(
    at = c(0, 20, 40, 60, 80, 100),
    # title = "Percentage of HM BB (%)", 
    title = '',
    title_gp = gpar(fontsize = legend_title_size),
    legend_height = unit(legend_height, "cm"),
    grid_width = unit(legend_width, "cm"),
    border='black',
    lwd=1.7,
    labels_gp = gpar(fontsize = legend_label_size),
    title_position = "leftcenter-rot"
  ), 
  row_labels = new_row_labels,
  column_labels = new_column_labels, 
  rect_gp = gpar(col = 'transparent',
                 lwd = 0
  ), use_raster = TRUE
)
p_binding_HM_BB

#### Figure 1I ####
#### A figure for total concentration of complexes ####

data_binding_wide <- data_binding %>% 
  mutate(total_complexes = cAA + cAB + cBB, 
         total_activity = 0.1*cA + 0.1*cB + 1*cAA + 1*cAB + 1*cBB) %>%
  select(new_binding_AB, new_binding_BB, total_activity) %>%
  arrange(desc(new_binding_AB)) %>%
  ## Limit to -15 to -5 to see the most interesting region
  filter(new_binding_AB <= -5, new_binding_AB >= -15,
         new_binding_BB <= -5, new_binding_BB >= -15) %>%
  pivot_wider(names_from = new_binding_BB, values_from = total_activity)


# Move the first column to the index
first_column <- data_binding_wide$new_binding_AB
data_binding_wide <- as.matrix(data_binding_wide[, 2:ncol(data_binding_wide)])
rownames(data_binding_wide) <- first_column

#### Draw the heatmap ####
# Define column labels
column_labels = seq(from = -15, to = -5.01, by = 0.01)
new_column_labels = column_labels
for(i in 1:length(column_labels)){
  
  cur_val = column_labels[i]
  new_val = ifelse(mod(cur_val, 2) == 0 && cur_val != -15,
                   toString(cur_val), '')
  
  new_column_labels[i] = new_val
  
}

# Define row labels
row_labels = rev(seq(from = -15, to = -5.01, by = 0.01))
new_row_labels = row_labels
for(i in 1:length(row_labels)){
  
  cur_val = row_labels[i]
  new_val = ifelse(mod(cur_val, 2) == 0 && cur_val != -15,
                   toString(cur_val), '')
  
  new_row_labels[i] = new_val
  
}

p_binding_total_complexes <- Heatmap(
  data_binding_wide,
  cluster_columns = F, cluster_rows = F, 
  col = colorRamp2(
    breaks = c(0, 20, 40, 60, 80, 100),
    colors = magma(6)),
  show_column_names = T, row_names_side = 'left',
  width=unit(heatmap_width, 'cm'), height = unit(heatmap_height, 'cm'),
  border = T,
  row_title = expression(paste(bold('\u0394'),
                               bold(G['bind,HET AB']), 
                               sep = '')),
  column_title = expression(paste(bold('\u0394'),
                                  bold(G['bind,HM BB']), 
                                  sep = '')),
  row_title_gp = gpar(fontsize=row_title_size, fontface = 'bold'),
  column_title_gp = gpar(fontsize=column_title_size, fontface = 'bold'),
  column_title_side = 'bottom',
  row_names_rot = 90, 
  row_names_centered = T,
  row_names_gp = gpar(fontsize=row_name_size, fontface = 'bold', align = 'left'),
  column_names_gp = gpar(fontsize=column_name_size, fontface='bold', align = 'left'),
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(
    at = c(0, 20, 40, 60, 80, 100),
    # title = 'Total activity',
    title = '', 
    title_gp = gpar(fontsize = legend_title_size),
    legend_height = unit(legend_height, "cm"),
    grid_width = unit(legend_width, "cm"),
    border='black',
    lwd=1.7,
    labels_gp = gpar(fontsize = legend_label_size),
    title_position = "leftcenter-rot"
  ), 
  row_labels = new_row_labels,
  column_labels = new_column_labels, 
  rect_gp = gpar(col = 'transparent',
                 lwd = 0
  ), use_raster = TRUE
)
p_binding_total_complexes

#### Figure 1J ####
#### Add a panel for the fitness with alpha = 80 ####

alpha = 80
beta = 0.5

data_binding_wide <- data_binding %>% mutate(total_complexes = cAA + cAB + cBB) %>%
  mutate(total_activity = cAA + cAB + cBB + (0.1*cA) + (0.1*cB)) %>%
  mutate(logw = round(log2(beta^(log2(total_activity / alpha)^2)), 5)) %>%
  filter(new_binding_AB <= -5, new_binding_AB >= -15,
         new_binding_BB <= -5, new_binding_BB >= -15) %>%
  select(new_binding_AB, new_binding_BB, logw) %>%
  arrange(desc(new_binding_AB)) %>%
  pivot_wider(names_from = new_binding_BB, values_from = logw)

# Move the first column to the index
first_column <- data_binding_wide$new_binding_AB
data_binding_wide <- as.matrix(data_binding_wide[, 2:ncol(data_binding_wide)])
rownames(data_binding_wide) <- first_column

#### Draw the heatmap ####
# Define column labels
column_labels = seq(from = -15, to = -5.01, by = 0.01)
new_column_labels = column_labels
for(i in 1:length(column_labels)){
  
  cur_val = column_labels[i]
  new_val = ifelse(mod(cur_val, 2) == 0 && cur_val != -15,
                   toString(cur_val), '')
  
  new_column_labels[i] = new_val
  
}

# Define row labels
row_labels = rev(seq(from = -15, to = -5.01, by = 0.01))
new_row_labels = row_labels
for(i in 1:length(row_labels)){
  
  cur_val = row_labels[i]
  new_val = ifelse(mod(cur_val, 2) == 0 && cur_val != -15,
                   toString(cur_val), '')
  
  new_row_labels[i] = new_val
  
}

p_binding_logw <- Heatmap(
  data_binding_wide,
  cluster_columns = F, cluster_rows = F, 
  col = colorRamp2(
    breaks = c(-1, -0.75, -0.5, -0.25, 0),
    colors = cividis(5)
  ),
  show_column_names = T, row_names_side = 'left',
  width=unit(heatmap_width, 'cm'), height = unit(heatmap_height, 'cm'),
  border = T,
  row_title = expression(paste(bold('\u0394'),
                               bold(G['bind,HET AB']), 
                               sep = '')),
  column_title = expression(paste(bold('\u0394'),
                                  bold(G['bind,HM BB']), 
                                  sep = '')),
  row_title_gp = gpar(fontsize=row_title_size, fontface = 'bold'),
  column_title_gp = gpar(fontsize=column_title_size, fontface = 'bold'),
  column_title_side = 'bottom',
  row_names_rot = 90, 
  row_names_centered = T,
  row_names_gp = gpar(fontsize=row_name_size, fontface = 'bold', align = 'left'),
  column_names_gp = gpar(fontsize=column_name_size, fontface='bold', align = 'left'),
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(
    # title = 'log2(Fitness)',
    title = '',
    title_gp = gpar(fontsize = legend_title_size),
    legend_height = unit(legend_height, "cm"),
    grid_width = unit(legend_width, "cm"),
    border='black',
    lwd=1.7,
    labels_gp = gpar(fontsize = legend_label_size),
    title_position = "leftcenter-rot"
  ), 
  row_labels = new_row_labels,
  column_labels = new_column_labels,
  rect_gp = gpar(col = 'transparent',
                 lwd = 0
  ), use_raster = TRUE
)
p_binding_logw

#### Put panels together ####

p_fig1_top <- plot_grid(fig_1A, fig_1B, ncol = 2, labels = c('A', 'B'), 
                        label_fontface = 'bold', label_size = panel_label_size)

p_fig1_bottom <- p_all <- plot_grid(draw_CHeatmap(p_stab_HET), draw_CHeatmap(p_stab_HM_BB), draw_CHeatmap(p_stab_total_complexes), draw_CHeatmap(p_stab_logw),
                                    draw_CHeatmap(p_binding_HET), draw_CHeatmap(p_binding_HM_BB), draw_CHeatmap(p_binding_total_complexes),draw_CHeatmap(p_binding_logw),
                                    nrow = 2, labels = c('C', 'D', 'E','F',
                                                         'G', 'H', 'I', 'J'),
                                    label_size = panel_label_size, label_fontface = 'bold')

p_fig1 <- plot_grid(p_fig1_top, p_fig1_bottom, nrow = 2)

ggsave(p_fig1, device = cairo_pdf, width = 19, height = 19, dpi = 300, units = 'cm',
       filename = 'Figures/Main_figures/1.Fig1_all_values_totalComplexes_fitness.pdf')

#### Supp. figure 1 ####
#### Diagram of the unintuitive parts of the solution space ####

## New variables for heatmaps (adjusting for new panel size) ##
row_title_size = 12
column_title_size = 12
row_name_size = 10
column_name_size = 10
legend_title_size = 10
legend_label_size = 9
heatmap_width = 6
heatmap_height = 6
legend_height = 4
legend_width = 0.5
ht_opt$TITLE_PADDING = unit(c(8.5, 8.5), "points")


data_binding_wide <- data_binding %>% select(new_binding_AB, new_binding_BB, pct_AB) %>%
  arrange(desc(new_binding_AB)) %>%
  ## Limit to -15 to -5 to see the most interesting region
  filter(new_binding_AB <= -5, new_binding_AB >= -15,
         new_binding_BB <= -5, new_binding_BB >= -15) %>%
  pivot_wider(names_from = new_binding_BB, values_from = pct_AB)

# Move the first column to the index
first_column <- data_binding_wide$new_binding_AB
data_binding_wide <- as.matrix(data_binding_wide[, 2:ncol(data_binding_wide)])
rownames(data_binding_wide) <- first_column

#### Draw the heatmap ####
# Define column labels
column_labels = seq(from = -15, to = -5.01, by = 0.01)
new_column_labels = column_labels
for(i in 1:length(column_labels)){
  
  cur_val = column_labels[i]
  new_val = ifelse(mod(cur_val, 2) == 0 && cur_val != -15,
                   toString(cur_val), '')
  
  new_column_labels[i] = new_val
  
}

# Define row labels
row_labels = rev(seq(from = -15, to = -5.01, by = 0.01))
new_row_labels = row_labels
for(i in 1:length(row_labels)){
  
  cur_val = row_labels[i]
  new_val = ifelse(mod(cur_val, 2) == 0 && cur_val != -15,
                   toString(cur_val), '')
  
  new_row_labels[i] = new_val
  
}

p_binding_HET <- Heatmap(
  data_binding_wide,
  cluster_columns = F, cluster_rows = F, 
  col = colorRamp2(
    breaks = c(0, 20, 40, 60, 80, 100),
    colors = viridis(6)),
  show_column_names = T, row_names_side = 'left',
  width=unit(heatmap_width, 'cm'), height = unit(heatmap_height, 'cm'),
  border = T,
  row_title = expression(paste(bold('\u0394'),
                               bold(G['bind, HET AB']), 
                               bold(' (kcal/mol)'), sep = '')),
  column_title = expression(paste(bold('\u0394'),
                                  bold(G['bind,HM BB']), 
                                  bold(' (kcal/mol)'), sep = '')),
  row_title_gp = gpar(fontsize=row_title_size, fontface = 'bold'),
  column_title_gp = gpar(fontsize=column_title_size, fontface = 'bold'),
  column_title_side = 'bottom',
  row_names_rot = 90, 
  row_names_centered = T,
  row_names_gp = gpar(fontsize=row_name_size, fontface = 'bold', align = 'left'),
  column_names_gp = gpar(fontsize=column_name_size, fontface='bold', align = 'left'),
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(
    at = c(0, 20, 40, 60, 80, 100),
    title = "Percentage of HET AB (%)", 
    title_gp = gpar(fontsize = legend_title_size),
    legend_height = unit(legend_height, "cm"),
    grid_width = unit(legend_width, "cm"),
    border='black',
    lwd=1.7,
    labels_gp = gpar(fontsize = legend_label_size),
    title_position = "leftcenter-rot"
  ), 
  row_labels = new_row_labels,
  column_labels = new_column_labels, 
  rect_gp = gpar(col = 'transparent',
                 lwd = 0
  ), use_raster = TRUE, 
  name = 'p_binding_HET'
)
p_binding_HET

## Save an annotated version of panel A
cairo_pdf(filename = 'Figures/Supplementary_figures/Fig_S1A.pdf', width = 12, height = 12)

p_binding_HET

decorate_heatmap_body("p_binding_HET", {
  grid.lines(x = c(0, 1), y = c(0.5, 0.5), gp = gpar(lty = 2, lwd = 3, col = 'red'))
  
  ## Alternative option with the numbers
  grid.text(label = '1)', x = 0.5, y = 0.5, gp = gpar(fontface = 'bold', fontsize = 10))
  grid.text(label = '2)', x = 0.05, y = 0.95, gp = gpar(fontface = 'bold', fontsize = 10, col = 'white'))
  grid.text(label = '3)', x = 0.95, y = 0.05, gp = gpar(fontface = 'bold', fontsize = 10))
  grid.text(label = '4)', x = 0.05, y = 0.45, gp = gpar(fontface = 'bold', fontsize = 10, col = 'white'))
  grid.text(label = '5)', x = 0.95, y = 0.55, gp = gpar(fontface = 'bold', fontsize = 10))
}
)

dev.off()

#### Prepare axes for the rest of the panels ####
## Cartoons of homomers and heteromers were added manually

p_figS1_axes <- ggplot() +
  theme(axis.line.x = element_line(), 
        axis.line.y = element_line(), 
        axis.ticks.y = element_line(), 
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(), 
        panel.grid.major = element_blank()
  ) +
  ylim(-5, 105) +
  scale_x_continuous(breaks = c(1, 2, 3), limits = c(0.5, 3.5)) +
  xlab('Complex') +
  ylab('Percentage (%)')
p_figS1_axes

## Fig. S1B: Starting point
p_figS1B <- p_figS1_axes + 
  theme(panel.grid.major = element_blank(), 
        plot.title = element_text(hjust = 0.5, size = 6)) +
  ggtitle(expression(atop(bold('Scenario 1): '), paste(
                           plain('\u0394'), plain(G['bind,HM AA']), 
                           plain(' = '), bold('\u0394'), bold(G['bind,HET AB']), 
                           plain(' = \u0394'), plain(G['bind,HM BB']), 
                           sep = ''))
  )
  )
p_figS1B

## Fig. S1C: HET has the weakest binding affinity
p_figS1C <- p_figS1_axes + 
  theme(panel.grid.major = element_blank(), 
        plot.title = element_text(hjust = 0.5, size = 7)) +
  ggtitle(expression(atop(bold('Scenario 2): '), paste(
                           bold('\u0394'), bold(G['bind,HET AB']), 
                           plain(' >> \u0394'), plain(G['bind,HM AA']), 
                           plain(' >> \u0394'), plain(G['bind,HM BB']),
                           sep = ''))
  )
  )
p_figS1C

## Fig. S1D: HET has the strongest binding affinity
p_figS1D <- p_figS1_axes + 
  theme(panel.grid.major = element_blank(), 
        plot.title = element_text(hjust = 0.5, size = 7)) +
  ggtitle(expression(atop(bold('Scenario 3): '), paste(
                           bold('\u0394'), bold(G['bind,HET AB']), 
                           plain(' << \u0394'), plain(G['bind,HM AA']), 
                           plain(' << \u0394'), plain(G['bind,HM BB']), 
                           sep = ''))
  )
  )
p_figS1D

## Fig. S1E
p_figS1E <- p_figS1_axes + 
  theme(panel.grid.major = element_blank(), 
        plot.title = element_text(hjust = 0.5, size = 7)) +
  ggtitle(expression(atop(bold('Scenario 4): '), paste(
    plain('\u0394'),
                           plain(G['bind,HM AA']),
                           plain(' > '), bold('\u0394'), bold(G['bind,HET AB']),
                           plain(' >> \u0394'), plain(G['bind,HM BB']), 
                           sep = ''))
  )
  )
p_figS1E

## Fig. S1F
p_figS1F <- p_figS1_axes + 
  theme(panel.grid.major = element_blank(), 
        plot.title = element_text(hjust = 0.5, size = 7)) +
  ggtitle(expression(atop(bold('Scenario 5): '), paste(
    plain('\u0394'),
                           plain(G['bind,HM AA']), 
                           plain(' < '), bold('\u0394'), bold(G['bind,HET AB']), 
                           plain(' << \u0394'), plain(G['bind,HM BB']), 
                           sep = ''))
  )
  )
p_figS1F


p_figS1 <- plot_grid(NULL, p_figS1B, p_figS1C,
                     p_figS1D, p_figS1E, p_figS1F, 
                     nrow = 2, ncol = 3, labels = c('A', 'B', 'C','D', 'E', 'F'), 
                     label_size = panel_label_size, label_fontface = 'bold'
)
ggsave(p_figS1, device = cairo_pdf, width = 20, height = 12, units = 'cm',  dpi = 300, 
       filename = 'Figures/Supplementary_figures/Appendix_Fig_S1_noPanelA_Empty.pdf')

#### Figure S2 ####

# Load the data on binding affinity
data_binding <- read_delim('Data/Results_solution_space/results_bindingEnergy_completeRange.tsv', delim = '\t')

# Add columns for the percentages of complexes, without considering monomers
data_binding %<>% mutate(total_complexes = (cAA + cBB + cAB)) %>%
  mutate(pct_AA = 100 * cAA / (cAA + cBB + cAB),
         pct_AB = 100 * cAB / (cAA + cBB + cAB),
         pct_BB = 100 * cBB / (cAA + cBB + cAB)) %>%
  mutate(pct_AA = ifelse(total_complexes > 80, NA, pct_AA), 
         pct_AB = ifelse(total_complexes > 80, NA, pct_AB), 
         pct_BB = ifelse(total_complexes > 80, NA, pct_BB)
  )

data_binding_wide <- data_binding %>% 
  mutate(total_complexes = cAA + cBB + cAB) %>%
  select(new_binding_AB, new_binding_AA, pct_AB) %>%
  arrange(desc(new_binding_AB)) %>%
  pivot_wider(names_from = new_binding_AA, values_from = pct_AB)

# Move the first column to the index
first_column <- data_binding_wide$new_binding_AB
data_binding_wide <- as.matrix(data_binding_wide[, 2:ncol(data_binding_wide)])
rownames(data_binding_wide) <- first_column

#### Draw the heatmap ####
# Define column labels
column_labels = seq(from = -15, to = 14.9, by = 0.1)
new_column_labels = column_labels
for(i in 1:length(column_labels)){
  
  cur_val = column_labels[i]
  new_val = ifelse(mod(cur_val, 4) == 0 && cur_val != -15,
                   toString(cur_val), '')
  
  new_column_labels[i] = new_val
  
}

# Define row labels
row_labels = rev(seq(from = -15, to = 14.9, by = 0.1))

new_row_labels = row_labels
for(i in 1:length(row_labels)){
  
  cur_val = row_labels[i]
  new_val = ifelse(mod(cur_val, 4) == 0 && cur_val != -15,
                   toString(cur_val), '')
  
  new_row_labels[i] = new_val
  
}

p_binding_HET <- Heatmap(
  data_binding_wide,
  cluster_columns = F, cluster_rows = F, 
  col = colorRamp2(
    breaks = c(0, 20, 40, 60, 80, 100),
    colors = viridis(6)),
  show_column_names = T, row_names_side = 'left',
  width=unit(heatmap_width, 'cm'), height = unit(heatmap_height, 'cm'),
  border = T,
  row_title = expression(paste(bold('\u0394'),
                               bold(G['bind, HET AB']), 
                               bold(' (kcal/mol)'), sep = '')),
  column_title = expression(paste(bold('\u0394'),
                                  bold(G['bind,HM AA']), 
                                  bold(' (kcal/mol)'), sep = '')),
  row_title_gp = gpar(fontsize=row_title_size, fontface = 'bold'),
  column_title_gp = gpar(fontsize=column_title_size, fontface = 'bold'),
  column_title_side = 'bottom',
  row_names_rot = 90, 
  row_names_centered = T,
  row_names_gp = gpar(fontsize=row_name_size, fontface = 'bold', align = 'left'),
  column_names_gp = gpar(fontsize=column_name_size, fontface='bold', align = 'left'),
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(
    at = c(0, 20, 40, 60, 80, 100),
    title = "Percentage of HET AB (%)", 
    title_gp = gpar(fontsize = legend_title_size),
    legend_height = unit(legend_height, "cm"),
    grid_width = unit(legend_width, "cm"),
    border='black',
    lwd=1.7,
    labels_gp = gpar(fontsize = legend_label_size),
    title_position = "leftcenter-rot"
  ), 
  row_labels = new_row_labels,
  column_labels = new_column_labels, 
  rect_gp = gpar(col = 'transparent',
                 lwd = 0
  ), use_raster = TRUE, 
  name = 'p_binding_HET'
)
p_binding_HET

#### Percentage of HM BB ####

data_binding_wide <- data_binding %>% select(new_binding_AB, new_binding_AA, pct_BB) %>%
  arrange(desc(new_binding_AB)) %>%
  pivot_wider(names_from = new_binding_AA, values_from = pct_BB)

# Move the first column to the index
first_column <- data_binding_wide$new_binding_AB
data_binding_wide <- as.matrix(data_binding_wide[, 2:ncol(data_binding_wide)])
rownames(data_binding_wide) <- first_column

#### Draw the heatmap ####
# Define column labels
column_labels = seq(from = -15, to = 14.9, by = 0.1)

new_column_labels = column_labels
for(i in 1:length(column_labels)){
  
  cur_val = column_labels[i]
  new_val = ifelse(mod(cur_val, 4) == 0 && cur_val != -15,
                   toString(cur_val), '')
  
  new_column_labels[i] = new_val
  
}

# Define row labels
row_labels = rev(seq(from = -15, to = 14.9, by = 0.1))

new_row_labels = row_labels
for(i in 1:length(row_labels)){
  
  cur_val = row_labels[i]
  new_val = ifelse(mod(cur_val, 4) == 0 && cur_val != -15,
                   toString(cur_val), '')
  
  new_row_labels[i] = new_val
  
}

p_binding_HM_BB <- Heatmap(
  data_binding_wide,
  cluster_columns = F, cluster_rows = F, 
  col = colorRamp2(
    breaks = c(0, 20, 40, 60, 80, 100),
    colors = viridis(6)),
  show_column_names = T, row_names_side = 'left',
  width=unit(heatmap_width, 'cm'), height = unit(heatmap_height, 'cm'),
  border = T,
  row_title = expression(paste(bold('\u0394'),
                               bold(G['bind,HET AB']), 
                               bold(' (kcal/mol)'), sep = '')),
  column_title = expression(paste(bold('\u0394'),
                                  bold(G['bind,HM AA']), 
                                  bold(' (kcal/mol)'), sep = '')),
  row_title_gp = gpar(fontsize=row_title_size, fontface = 'bold'),
  column_title_gp = gpar(fontsize=column_title_size, fontface = 'bold'),
  column_title_side = 'bottom',
  row_names_rot = 90, 
  row_names_centered = T,
  row_names_gp = gpar(fontsize=row_name_size, fontface = 'bold', align = 'left'),
  column_names_gp = gpar(fontsize=column_name_size, fontface='bold', align = 'left'),
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(
    at = c(0, 20, 40, 60, 80, 100),
    title = "Percentage of HM BB (%)", 
    title_gp = gpar(fontsize = legend_title_size),
    legend_height = unit(legend_height, "cm"),
    grid_width = unit(legend_width, "cm"),
    border='black',
    lwd=1.7,
    labels_gp = gpar(fontsize = legend_label_size),
    title_position = "leftcenter-rot"
  ), 
  row_labels = new_row_labels,
  column_labels = new_column_labels, 
  rect_gp = gpar(col = 'transparent',
                 lwd = 0
  ), use_raster = TRUE
)
p_binding_HM_BB

#### A figure for total concentration of complexes ####

data_binding_wide <- data_binding %>% 
  mutate(total_activity = 0.1*cA + 0.1*cB + 1*cAA + 1*cAB + 1*cBB) %>%
  mutate(total_activity = ifelse(or(total_complexes > 80, total_activity == 0), NA, total_activity)
  ) %>%
  select(new_binding_AB, new_binding_AA, total_activity) %>%
  arrange(desc(new_binding_AB)) %>%
  pivot_wider(names_from = new_binding_AA, values_from = total_activity)


# Move the first column to the index
first_column <- data_binding_wide$new_binding_AB
data_binding_wide <- as.matrix(data_binding_wide[, 2:ncol(data_binding_wide)])
rownames(data_binding_wide) <- first_column

#### Draw the heatmap ####
# Define column labels
column_labels = seq(from = -15, to = 14.9, by = 0.1)

new_column_labels = column_labels
for(i in 1:length(column_labels)){
  
  cur_val = column_labels[i]
  new_val = ifelse(mod(cur_val, 4) == 0 && cur_val != -15,
                   toString(cur_val), '')
  
  new_column_labels[i] = new_val
  
}

# Define row labels
row_labels = rev(seq(from = -15, to = 14.9, by = 0.1))

new_row_labels = row_labels
for(i in 1:length(row_labels)){
  
  cur_val = row_labels[i]
  new_val = ifelse(mod(cur_val, 4) == 0 && cur_val != -15,
                   toString(cur_val), '')
  
  new_row_labels[i] = new_val
  
}

p_binding_total_complexes <- Heatmap(
  data_binding_wide,
  cluster_columns = F, cluster_rows = F, 
  col = colorRamp2(
    breaks = c(0, 20, 40, 60, 80, 100),
    colors = magma(6)),
  show_column_names = T, row_names_side = 'left',
  width=unit(heatmap_width, 'cm'), height = unit(heatmap_height, 'cm'),
  border = T,
  row_title = expression(paste(bold('\u0394'),
                               bold(G['bind,HET AB']), 
                               bold(' (kcal/mol)'), sep = '')),
  column_title = expression(paste(bold('\u0394'),
                                  bold(G['bind,HM AA']), 
                                  bold(' (kcal/mol)'), sep = '')),
  row_title_gp = gpar(fontsize=row_title_size, fontface = 'bold'),
  column_title_gp = gpar(fontsize=column_title_size, fontface = 'bold'),
  column_title_side = 'bottom',
  row_names_rot = 90, 
  row_names_centered = T,
  row_names_gp = gpar(fontsize=row_name_size, fontface = 'bold', align = 'left'),
  column_names_gp = gpar(fontsize=column_name_size, fontface='bold', align = 'left'),
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(
    at = c(0, 20, 40, 60, 80, 100),
    title = 'Total activity',
    title_gp = gpar(fontsize = legend_title_size),
    legend_height = unit(legend_height, "cm"),
    grid_width = unit(legend_width, "cm"),
    border='black',
    lwd=1.7,
    labels_gp = gpar(fontsize = legend_label_size),
    title_position = "leftcenter-rot"
  ), 
  row_labels = new_row_labels,
  column_labels = new_column_labels, 
  rect_gp = gpar(col = 'transparent',
                 lwd = 0
  ), use_raster = TRUE
)
p_binding_total_complexes

#### Add a panel for the fitness with alpha = 80 ####

alpha = 80
beta = 0.5

data_binding_wide <- data_binding %>% mutate(total_complexes = cAA + cAB + cBB) %>%
  mutate(total_activity = cAA + cAB + cBB + (0.1*cA) + (0.1*cB)) %>%
  mutate(logw = round(log2(beta^(log2(total_activity / alpha)^2)), 5)) %>%
  mutate(logw = ifelse(or(total_complexes > 80, total_activity == 0), NA, logw)) %>%
  select(new_binding_AB, new_binding_AA, logw) %>%
  arrange(desc(new_binding_AB)) %>%
  pivot_wider(names_from = new_binding_AA, values_from = logw)

# Move the first column to the index
first_column <- data_binding_wide$new_binding_AB
data_binding_wide <- as.matrix(data_binding_wide[, 2:ncol(data_binding_wide)])
rownames(data_binding_wide) <- first_column

#### Draw the heatmap ####
# Define column labels
column_labels = seq(from = -15, to = 14.9, by = 0.1)

new_column_labels = column_labels
for(i in 1:length(column_labels)){
  
  cur_val = column_labels[i]
  new_val = ifelse(mod(cur_val, 4) == 0 && cur_val != -15,
                   toString(cur_val), '')
  
  new_column_labels[i] = new_val
  
}

# Define row labels
row_labels = rev(seq(from = -15, to = 14.9, by = 0.1))

new_row_labels = row_labels
for(i in 1:length(row_labels)){
  
  cur_val = row_labels[i]
  new_val = ifelse(mod(cur_val, 4) == 0 && cur_val != -15,
                   toString(cur_val), '')
  
  new_row_labels[i] = new_val
  
}

p_binding_logw <- Heatmap(
  data_binding_wide,
  cluster_columns = F, cluster_rows = F, 
  col = colorRamp2(
    breaks = c(-1, -0.75, -0.5, -0.25, 0),
    colors = cividis(5)
  ),
  show_column_names = T, row_names_side = 'left',
  width=unit(heatmap_width, 'cm'), height = unit(heatmap_height, 'cm'),
  border = T,
  row_title = expression(paste(bold('\u0394'),
                               bold(G['bind,HET AB']), 
                               bold(' (kcal/mol)'), sep = '')),
  column_title = expression(paste(bold('\u0394'),
                                  bold(G['bind,HM AA']), 
                                  bold(' (kcal/mol)'), sep = '')),
  row_title_gp = gpar(fontsize=row_title_size, fontface = 'bold'),
  column_title_gp = gpar(fontsize=column_title_size, fontface = 'bold'),
  column_title_side = 'bottom',
  row_names_rot = 90, 
  row_names_centered = T,
  row_names_gp = gpar(fontsize=row_name_size, fontface = 'bold', align = 'left'),
  column_names_gp = gpar(fontsize=column_name_size, fontface='bold', align = 'left'),
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(
    title = 'log2(Fitness)',
    title_gp = gpar(fontsize = legend_title_size),
    legend_height = unit(legend_height, "cm"),
    grid_width = unit(legend_width, "cm"),
    border='black',
    lwd=1.7,
    labels_gp = gpar(fontsize = legend_label_size),
    title_position = "leftcenter-rot"
  ), 
  row_labels = new_row_labels,
  column_labels = new_column_labels,
  rect_gp = gpar(col = 'transparent',
                 lwd = 0
  ), use_raster = TRUE
)
p_binding_logw

p_heatmaps_full_range <- plot_grid(draw_CHeatmap(p_binding_HET), draw_CHeatmap(p_binding_HM_BB),
                                   draw_CHeatmap(p_binding_total_complexes), draw_CHeatmap(p_binding_logw), 
                                   ncol = 2, nrow = 2, labels = c('A', 'B', 'C', 'D'), label_size = panel_label_size, 
                                   label_fontface = 'bold'
)

ggsave(p_heatmaps_full_range, device = cairo_pdf, width = 12, height = 9, units = 'cm', dpi = 300, 
       filename = 'Figures/Supplementary_figures/Appendix_FigS02.solution_space_allBindingEnergies_complete.pdf')


#### Figure 2: Parametric simulations ####

## New data with pdup = 1, alpha = 80
all_results_parametric_sd <- read_delim('Results_simulations/Results_pdup1/Parametric_simulations/008_simulations_parametric_sd_final_80opt/all_results_all_sims.tsv', 
                                        delim = '\t')

## Remove some unneeded columns
all_results_parametric_sd %<>% 
  select(-sA, -dA, -aA, -dAA, -aAA, -sB, -dB, -aB, -dBB, -aBB, -dAB)

all_results_parametric_sd %<>% 
  mutate(fixed_mut = rep(1:200, nrow(all_results_parametric_sd) / 200))

# Calculate the percentage of HET
all_results_parametric_sd %<>% filter(!(is.na(aAB))) %>%
  mutate(Total_complexes = cAA + cAB + cBB + cA + cB) %>% # Add the total concentration of complexes
  mutate(pct_HET = round((cAB *100 / Total_complexes), 2), 
         pct_HMs = round((cAA + cBB) * 100 / Total_complexes, 2), 
         pct_mono = round((cA + cB) * 100 / Total_complexes, 2)
  )

# Label the different replicates depending on where they finish
all_results_parametric_final_points_sd <- all_results_parametric_sd %>% filter(fixed_mut == 200) %>%
  mutate(Outcome = ifelse(pct_HET >= 70, 'HET dominant', 
                          ifelse(pct_HMs >= 70, 'HM dominant', 
                                 ifelse(pct_mono >= 70, 'Monomers', 
                                        ifelse((pct_HET + pct_HMs) >= 70, 'Both HM and HET', 
                                               'Ambiguous')
                                 )
                          )
  )
  )

# Add the info on the final points
all_results_parametric_outcome_sd <- left_join(x = all_results_parametric_sd,
                                               y = all_results_parametric_final_points_sd %>% select(Replicate, intHM_param, sdHM, Outcome), 
                                               by = c('Replicate' = 'Replicate', 'intHM_param' = 'intHM_param', 'sdHM' = 'sdHM'))


#### Figure S3 ####
p_figS3 <- all_results_parametric_outcome_sd %>%
  mutate(Outcome = factor(Outcome, levels = c('HET dominant', 'Both HM and HET', 'HM dominant', 'Monomers', 'Ambiguous'))) %>%
  ggplot(aes(x = fixed_mut, y = pct_HET, group = Replicate, colour = Outcome)) +
  geom_line() +
  scale_colour_manual(values = c('#4575B4', '#008000', '#FC8D59', 'black', '#B300B3')) +
  facet_grid(sdHM~intHM_param) +
  xlab('Fixed mutations') + ylab('HET percentage (%)') +
  theme(legend.position = 'top', legend.justification = 'center', 
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(angle = 90, size = 8, hjust = 1, vjust = 0.5), 
        axis.text.y = element_text(size = 8), 
        strip.text = element_text(size = 10)) +
  labs(colour = '')
p_figS3
ggsave(p_figS3, device = cairo_pdf, width = 20, height = 20, units = 'cm', dpi = 300,
       filename = 'Figures/Supplementary_figures/Appendix_FigS03.pdf')

## Heatmap for figure 2B ##
median_outcome_parametric_sd <- all_results_parametric_final_points_sd %>% ungroup() %>%
  group_by(intHM_param, intHET_param, sdHM, sdHET) %>%
  summarise(median_HET = median(pct_HET), 
            mean_HET = mean(pct_HET), 
            median_HMs = median(pct_HMs), 
            mean_HMs = mean(pct_HMs), 
            median_mono = median(pct_mono),
            mean_mono = mean(pct_mono)
  )
write.table(median_outcome_parametric_sd, append = F, quote = F, sep = '\t', row.names = F, col.names = T,
            file = 'Figures/Main_figures/data_fig_2B.tsv'
)

## Convert to a matrix
matrix_median_outcome_param_df <- median_outcome_parametric_sd %>% ungroup() %>%
  mutate(intHM_param = format(intHM_param, nsmall = 2)) %>%
  select(intHM_param, sdHM, mean_HET) %>%
  pivot_wider(names_from = intHM_param, values_from = mean_HET)

needed_rownames <- matrix_median_outcome_param_df$sdHM
matrix_median_outcome_param_df %<>% select(-sdHM) %>% as.matrix()
rownames(matrix_median_outcome_param_df) <- needed_rownames

## New parameters
heatmap_width_tmp = 5
heatmap_height_tmp = 5

row_title_size = 12
column_title_size = 12
row_name_size = 10
column_name_size = 10

legend_title_size = 10
legend_label_size = 8
legend_height = 3
ht_opt$TITLE_PADDING = unit(c(0.5, 0.5), "points")
legend_width = 0.15

p_fig2B <- Heatmap(matrix_median_outcome_param_df,
                   cluster_columns = F, cluster_rows = F, 
                   col = colorRamp2(
                     breaks = c(0, 20, 40, 60, 80, 100),
                     colors = viridis(6)
                   ),
                   show_column_names = T, row_names_side = 'left',
                   width=unit(heatmap_width_tmp, 'cm'), height = unit(heatmap_height_tmp, 'cm'),
                   border = T,
                   column_title = expression(paste(bold('Mean \u0394\u0394'), 
                                                   bold(G['bind,HM']), ' ',
                                                   bold('(kcal/mol)'), sep = '')
                   ),
                   row_title = expression(paste(bold('Std. dev. \u0394\u0394'), 
                                                bold(G['bind,HM']), ' ',
                                                bold('(kcal/mol)'), sep = '')
                   ),
                   row_title_gp = gpar(fontsize=row_title_size, fontface = 'bold'),
                   column_title_gp = gpar(fontsize=column_title_size, fontface = 'bold'),
                   column_title_side = 'bottom',
                   row_names_rot = 0, 
                   row_names_centered = T,
                   row_names_gp = gpar(fontsize=row_name_size, fontface = 'bold', align = 'left'),
                   column_names_gp = gpar(fontsize=column_name_size, fontface='bold', align = 'left'),
                   show_heatmap_legend = TRUE,
                   heatmap_legend_param = list(
                     at = c(0, 20, 40, 60, 80, 100),
                     title = "Percentage of HET AB (%)", 
                     title_gp = gpar(fontsize = legend_title_size),
                     legend_height = unit(legend_height, "cm"),
                     grid_width = unit(legend_width, "cm"),
                     border='black',
                     lwd=1.7,
                     labels_gp = gpar(fontsize = legend_label_size),
                     title_position = "leftcenter-rot"
                   ), 
                   ## Draw dots to indicate paramaters being the same as HET and the observed ones
                   cell_fun = function(i, j, x, y, w, h, col) { # add text to each grid
                     if (i == 3 && j == 3){ ## HM would have the same parameters as the HET
                       grid.points(x, y, pch = 19, 
                                   size = unit(1, 'char'),
                                   gp = gpar(col = 'white'))
                     }else if(i == 7 && j == 7){ ## Observed parameters for the HM
                       grid.points(x, y, pch = 17, 
                                   size = unit(1, 'char'),
                                   gp = gpar(col = 'white'))
                     }
                   },  rect_gp = gpar(col = 'transparent',
                                      lwd = 0
                   ), use_raster = TRUE
)
p_fig2B

## Load the mutational effects
mut_effects <- read_delim('Data/all_mut_matrices.tsv', '\t')

## New data with pdup = 1, alpha = 80, regular simulations
all_results <- read_delim('Results_simulations/Results_pdup1/Simulations_with_structures/008_simulations_80opt/all_results_all_sims.tsv', delim = '\t')

## Remove some unneeded columns
all_results %<>% 
  select(-sA, -dA, -aA, -dAA, -aAA, -sB, -dB, -aB, -dBB, -aBB, -dAB, -aAB, -t)

## New data with pdup = 1, alpha = 80, biased simulations
all_results_biased <- read_delim('Results_simulations/Results_pdup1/Simulations_with_structures/008_biased_simulations_80opt/all_results_all_sims.tsv', delim = '\t')

all_results_biased %<>% 
  select(-sA, -dA, -aA, -dAA, -sB, -dB, -aB, -dBB, -dAB, -t)

structures_normal <- unique(all_results$Complex)
structures_biased <- unique(all_results_biased$Complex)

## Check which simulations to remove to have the same data everywhere
structures_normal[which(!(structures_normal %in% structures_biased))]
structures_biased[which(!(structures_biased %in% structures_normal))]

## The following structures should be removed because of redundancy / homology
# Homologous to 3gi4: 1bdr, 2z54, 1hps
# Homologous to 3w5z: 6h5v
# Homologous to 4nak: 5mrn
# Homologous but disagreement in the structures: 1lhv, 1f5f
structures_remove <- structures_normal[which(!(structures_normal %in% structures_biased))]

## Include these structures in the list of structures to remove,
## also include 2c9s, 3f0h, and 3f1l because of the FoldX bad distributions of
## mutational effects (simulated HM exactly like HET)
## 3qfl is the coiled coil that has too many destabilizing mutations
## 3cyl is homologous to 4rfp (48.8% sequence identity)
## 4jdr, 6jbj are structures added for some tests but do not clear our thresholds
structures_remove <- c(structures_remove, '1bdr', '2z54', '1hps', '6h5v', '5mrn',
                       '1lhv', '1f5f', '2c9s', '3f1l', '3f0h', '3qfl', 
                       '3cyl', '4jdr', '6jbj')

# Add a column indicating how many mutations have fixed up to that point
all_results %<>% group_by(Complex, Replicate, Start_stab, Start_aff) %>%
  mutate(fixed_mut = row_number()) %>%
  filter(!(Complex %in% structures_remove))

## Save the list of structures I am considering for the simulations
structures_keep <- as.data.frame(unique(all_results$Complex))
write.table(structures_keep, append = F, quote = F, sep = '\t', row.names = F, col.names = F, 
            file = 'Data/final_104_structures.txt')

## Keep only the structures that are still in the results matrix
mut_effects %<>% filter(Complex %in% all_results$Complex)


### Check the parameters of the distribution of mutational effects
## All structures together
mean(mut_effects$Mean_ddG_stab_HET)
# 2.596058

sd(mut_effects$Mean_ddG_stab_HET)
# 4.633035

mean(mut_effects$Mean_ddG_int_HM)
# 0.3506799

sd(mut_effects$Mean_ddG_int_HM)
# 2.287689

mean(mut_effects$Mean_ddG_int_HET)
# 0.1682575

sd(mut_effects$Mean_ddG_int_HET)
# 1.118576

## Check the average of the correlations when separating individual proteins
avg_corr <- mut_effects %>% group_by(Complex) %>%
  summarise(cor_intHM_intHET = cor(Mean_ddG_int_HM, Mean_ddG_int_HET), 
            cor_intHM_fold = cor(Mean_ddG_int_HM, Mean_ddG_stab_HET), 
            cor_intHET_fold = cor(Mean_ddG_int_HET, Mean_ddG_stab_HET))

mean(avg_corr$cor_intHM_intHET)
# 0.9196066

mean(avg_corr$cor_intHM_fold)
# 0.2604355

mean(avg_corr$cor_intHET_fold)
# 0.2677022

## Get the ratios of protein lengths to 200 mutations
prot_lengths <- read_delim('Data/Structures/final_dataset_lengths.tsv', delim = '\t')

which(!(tolower(prot_lengths$PDB) %in% structures_keep$`unique(all_results$Complex)`))

summary(prot_lengths$Ratio_muts_length)
summary(prot_lengths$Length)

#### Figure EV1 ####

## Show the correlation between effects on binding HM and binding HET
mut_effects_nonzero <- mut_effects %>% 
  filter(!(Mean_ddG_int_HET == 0), !(Mean_ddG_int_HM == 0))

p_binding_corr_center <- mut_effects_nonzero %>%
  ggplot(aes(x = Mean_ddG_int_HET, y = Mean_ddG_int_HM)) +
  geom_hex(aes(fill = log10(..count..))) +
  scale_fill_gradientn(limits = c(0, 6), 
                       colors = c('#000000', '#F6F1F1', '#19A7CE', '#146C94')) +
  stat_cor(method = 'pearson', p.accuracy = 0.001, r.accuracy = 0.01, cor.coef.name = 'r',
           label.y.npc = 0.005, label.x.npc = 0.42, size = 3) +
  geom_smooth(method = 'lm', alpha = 0.5) +
  theme(
    panel.grid.major = element_blank(), 
    axis.text = element_text(size = 10), 
    axis.title = element_text(size = 12), 
    legend.position = 'none', legend.justification = 0.5) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  xlim(-5, 10) + ylim(-5, 10) +
  labs(x = expression(paste( 
    bold('\u0394\u0394'), 
    bold(G['bind,HET']), sep = '')),
    y = expression(paste( 
      bold('\u0394\u0394'), 
      bold(G['bind,HM']), sep = '')), 
    fill = 'log10(Mutation count)'
  )
p_binding_corr_center  

p_density_right <- mut_effects_nonzero %>%
  ggplot(aes(y = Mean_ddG_int_HM)) +
  geom_density(alpha = 0.5
  ) +
  ylim(-5, 10) +
  theme(legend.position = 'none', 
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = margin(l = -0.7, r = 0, b = 0.85, t = 0.2, unit = 'cm')
  ) +
  xlab('') + ylab('')
p_density_right

p_density_top<- mut_effects_nonzero %>%
  ggplot(aes(x = Mean_ddG_int_HET)) +
  geom_density(alpha = 0.5
  ) +
  xlim(-5, 10) +
  theme(legend.position = 'none', 
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin = margin(t = 0, b = -0.5, l = 0.85, r = 0.15, unit = 'cm')
  ) +
  xlab('') + ylab('')
p_density_top

p_binding_corr <- plot_grid(p_density_top, NULL, 
                            p_binding_corr_center, p_density_right, ncol = 2, nrow = 2, 
                            rel_widths = c(1, 0.15), rel_heights = c(0.15, 1))
p_binding_corr

## Binding HM versus folding energy
corr_pearson_intHM_fold <- cor(mut_effects_nonzero$Mean_ddG_int_HM, mut_effects_nonzero$Mean_ddG_stab_HET)

p_bindingHM_folding_center <- mut_effects_nonzero %>%
  ggplot(aes(x = Mean_ddG_stab_HET, y = Mean_ddG_int_HM)) +
  geom_hex(aes(fill = log10(..count..))) +
  scale_fill_gradientn(limits = c(0, 6), 
                       colors = c('#000000', '#F6F1F1', '#19A7CE', '#146C94')) +
  stat_cor(method = 'pearson', p.accuracy = 0.001, r.accuracy = 0.01, cor.coef.name = 'r',
           label.y.npc = 0.005, label.x.npc = 0.42, size = 3) +
  geom_smooth(method = 'lm', alpha = 0.5) +
  theme(
    panel.grid.major = element_blank(),
    axis.text = element_text(size = 10), 
    axis.title = element_text(size = 12), 
    legend.position = 'none', legend.justification = 0.5, 
    legend.title = element_text(size = 9), 
    legend.text = element_text(size = 8)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  xlim(-5, 10) + ylim(-5, 10) +
  labs(x = expression(paste(
    bold('\u0394\u0394'), 
    bold(G['fold']), sep = '')),
    y = expression(paste( 
      bold('\u0394\u0394'), 
      bold(G['bind,HM']), sep = '')), 
    fill = 'log10(Mutation count)'
  )
p_bindingHM_folding_center  

p_density_right <- mut_effects_nonzero %>%
  ggplot(aes(y = Mean_ddG_int_HM)) +
  geom_density(alpha = 0.5
  ) +
  ylim(-5, 10) +
  theme(legend.position = 'none', 
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = margin(l = -0.7, r = 0, b = 0.85, t = 0.2, unit = 'cm')
  ) +
  xlab('') + ylab('')
p_density_right

p_density_top<- mut_effects_nonzero %>%
  ggplot(aes(x = Mean_ddG_stab_HET)) +
  geom_density(alpha = 0.5
  ) +
  xlim(-5, 10) +
  theme(legend.position = 'none', 
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(),
        axis.text.y = element_blank(),
        plot.margin = margin(t = 0, b = -0.5, l = 0.85, r = 0.15, unit = 'cm')
  ) +
  xlab('') + ylab('')
p_density_top

p_bindingHM_folding_corr <- plot_grid(p_density_top, NULL, 
                                      p_bindingHM_folding_center, p_density_right, ncol = 2, nrow = 2, 
                                      rel_widths = c(1, 0.15), rel_heights = c(0.15, 1))
p_bindingHM_folding_corr

## Binding HET vs folding effects
corr_pearson_intHET_fold <- cor(mut_effects_nonzero$Mean_ddG_int_HET, mut_effects_nonzero$Mean_ddG_stab_HET)

p_bindingHET_folding_center <- mut_effects_nonzero %>%
  ggplot(aes(x = Mean_ddG_stab_HET, y = Mean_ddG_int_HET)) +
  geom_hex(aes(fill = log10(..count..))) +
  scale_fill_gradientn(limits = c(0, 6), 
                       colors = c('#000000', '#F6F1F1', '#19A7CE', '#146C94')) +
  stat_cor(method = 'pearson', p.accuracy = 0.001, r.accuracy = 0.01, cor.coef.name = 'r',
           label.y.npc = 0.005, label.x.npc = 0.42, size = 3) +
  geom_smooth(method = 'lm', alpha = 0.5) +
  theme(
    panel.grid.major = element_blank(),
    axis.text = element_text(size = 10), 
    axis.title = element_text(size = 12), 
    legend.position = 'none', legend.justification = 0.5) +
  xlim(-5, 10) + ylim(-5, 10) +
  labs(x = expression(paste( 
    bold('\u0394\u0394'), 
    bold(G['fold']), sep = '')),
    y = expression(paste(
      bold('\u0394\u0394'), 
      bold(G['bind,HET']), sep = '')), 
    fill = 'log10(Mutation count)'
  ) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed')
p_bindingHET_folding_center  

p_density_right <- mut_effects_nonzero %>%
  ggplot(aes(y = Mean_ddG_int_HET)) +
  geom_density(alpha = 0.5
  ) +
  ylim(-5, 10) +
  theme(legend.position = 'none', 
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = margin(l = -0.7, r = 0, b = 0.85, t = 0.2, unit = 'cm')
  ) +
  xlab('') + ylab('')
p_density_right

p_density_top<- mut_effects_nonzero %>%
  ggplot(aes(x = Mean_ddG_stab_HET)) +
  geom_density(alpha = 0.5
  ) +
  xlim(-5, 10) +
  theme(legend.position = 'none', 
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(),
        axis.text.y = element_blank(),
        plot.margin = margin(t = 0, b = -0.5, l = 0.85, r = 0.15, unit = 'cm')
  ) +
  xlab('') + ylab('')
p_density_top

p_bindingHET_folding_corr <- plot_grid(p_density_top, NULL, 
                                       p_bindingHET_folding_center, p_density_right, ncol = 2, nrow = 2, 
                                       rel_widths = c(1, 0.15), rel_heights = c(0.15, 1))
p_bindingHET_folding_corr

p_figEV1_nolegend <- plot_grid(p_binding_corr + theme(aspect.ratio = 1, legend.position = 'none'),
                              p_bindingHM_folding + theme(aspect.ratio = 1, legend.position = 'none'),
                              p_bindingHET_folding + theme(aspect.ratio = 1, legend.position = 'none'),
                              labels = c('A', 'B', 'C'), label_size = panel_label_size, ncol = 3,
                              label_fontface = 'bold')

legend_hexbins = get_legend(p_bindingHM_folding + theme(legend.justification = 0.5, 
                                                        legend.key.height = unit(0.15, 'cm')
)
)

p_figEV1 <- plot_grid(legend_hexbins, p_figEV1_nolegend, nrow = 2,
                    rel_heights = c(0.15, 1))

ggsave(p_figEV1, device = cairo_pdf, width = 20, height = 8, dpi = 300, units = 'cm',
       filename = 'Figures/ExpandedView_figures/FigEV1_corr_mut_effects.pdf')


## Look at the final points to check the outcomes of simulations
all_data_final_points <- all_results %>% ungroup() %>%
  mutate(Total_complexes = cAA + cAB + cBB + cA + cB) %>% # Add the total concentration of complexes
  filter(fixed_mut == 200, Start_stab == -5, Start_aff == -10) %>%
  mutate(outcome = ifelse(cAB / Total_complexes >= 0.7, 'HET dominant', 
                          ifelse( (cAA + cBB) / Total_complexes >= 0.7, 'HM dominant', 
                                  ifelse( (cA + cB) / Total_complexes >= 0.7, 'Monomers', 
                                          ifelse( (cAA + cAB + cBB) / Total_complexes >= 0.7, 'Both HM and HET', 
                                                  'Ambiguous')
                                  )
                          )
  )
  )

all_data_summary <- all_results %>% 
  ungroup() %>%
  mutate(Total_complexes = cAA + cAB + cBB + cA + cB) %>% # Add the total concentration of complexes
  mutate(pct_HET = round((cAB *100 / Total_complexes), 2), 
         pct_HMs = round((cAA + cBB) * 100 / Total_complexes, 2), 
         pct_mono = round((cA + cB) * 100 / Total_complexes, 2)) %>%
  group_by(Start_stab, Start_aff, Complex, fixed_mut) %>%
  summarise(
    med_pct_HET = median(pct_HET), 
    mean_pct_HET = mean(pct_HET), 
    sem_pct_HET = sd(pct_HET) / sqrt(n()), 
    
    mean_pct_HMs = mean(pct_HMs), 
    sem_pct_HMs = sd(pct_HMs) / sqrt(n()),
    
    mean_pct_mono = mean(pct_mono), 
    sem_pct_mono = sd(pct_mono) / sqrt(n())
  )

# Get the order of the figures
all_data_summary_final_points <- all_data_summary %>% 
  filter(fixed_mut == 200, Start_stab == -5, Start_aff == -10) %>%
  arrange(desc(mean_pct_HET))

# Relevel according to the order from the other plot
all_data_summary %<>% 
  mutate(Complex = factor(toupper(Complex), levels = toupper(all_data_summary_final_points$Complex)))

#### Figure S6: Different sets of starting conditions ####
p_figS6 <- all_data_summary %>%
  mutate(Conditions = str_c('Start dG_stab = ', Start_stab, ', Start dG aff = ', Start_aff, sep = '')) %>%
  ggplot(aes(x = fixed_mut, y = med_pct_HET, colour = Conditions)) +
  xlab('Fixed mutations') + ylab('HET percentage (%)') +
  geom_line(linewidth = 2) +
  theme(axis.text.y = element_text(size = 24),
        axis.text.x = element_text(size = 24, angle = 90, vjust = 0.5, hjust = 1),
        axis.title = element_text(size = 26),
        legend.title = element_text(size = 26),
        legend.text = element_text(size = 24, hjust = 0.5),
        legend.position = 'top',legend.justification = 'center',
        plot.title = element_text(hjust = 0.5),
        strip.text = element_text(size = 28, face = 'bold'), 
        panel.grid.major = element_blank()) +
  facet_wrap(~Complex, nrow = 8) +
  guides(colour = guide_legend(nrow = 2, byrow = T)) +
  geom_hline(yintercept = 50, linewidth = 1, linetype = 'dashed')
p_figS6
ggsave(p_figS6, device = cairo_pdf, width = 45, height = 28, dpi = 500,
       filename = 'Figures/Supplementary_figures/Appendix_FigS06.pdf')

### Present results as stacked barplots
## Arrange and relevel
all_data_summary_final_points %<>% arrange(desc(mean_pct_HET))

all_data_summary_final_points %<>%
  mutate(Complex = factor(toupper(Complex), levels = toupper(all_data_summary_final_points$Complex)))

all_data_summary_final_points_long <- all_data_summary_final_points %>% ungroup() %>%
  pivot_longer(cols = c(mean_pct_HET, mean_pct_HMs, mean_pct_mono), 
               names_to = 'assembly', values_to = 'mean_value')

# Stacked barplot of outcomes per structure
p_fig2C <- all_data_summary_final_points_long %>% ungroup() %>%
  mutate(assembly_label = ifelse(assembly == 'mean_pct_HET', 'HET', 
                                 ifelse(assembly == 'mean_pct_HMs', 'HMs (AA+BB)',
                                        'Monomers (A+B)')), 
         Complex = toupper(Complex)) %>%
  mutate(assembly_label = factor(assembly_label, levels = c('HET', 'Monomers (A+B)', 'HMs (AA+BB)')),
         Complex = factor(Complex, levels = toupper(all_data_summary_final_points$Complex))
  ) %>%
  ggplot(aes(x = Complex, y = mean_value, fill = assembly_label)) +
  geom_bar(stat = 'identity') +
  xlab('PDB ID') + ylab('Percentage (%)') +
  labs(fill = '') +
  geom_hline(yintercept = 50, linetype = 'dashed') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 9),
        axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 9),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.position = 'top', 
        legend.justification = 'center') +
  guides(fill = guide_legend(override.aes = list(size = 0.25))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  scale_fill_manual(values = c('#4575B4', 'grey', '#FC8D59' ))
p_fig2C

#### Put all parts of the new figure 2 together ####

## Load figure 2A
p_fig2A <- ggdraw() + draw_image(
  'Figures/Diagrams/Diagram2A.png'
)
p_fig2A


p_fig2_top <- plot_grid(p_fig2A, draw_CHeatmap(p_fig2B), 
                        labels = c('A', 'B'), label_size = panel_label_size, ncol = 2,
                        label_fontface = 'bold', rel_widths = c(1, 0.8))

p_fig2 <- plot_grid(p_fig2_top, p_fig2C, 
                    labels = c('', 'C'), label_size = panel_label_size, nrow = 2,
                    label_fontface = 'bold')

ggsave(p_fig2, device = cairo_pdf, width = 28, height = 20, dpi = 300, units = 'cm',
       filename = 'Figures/Main_figures/2.Fig2.pdf')
## Cartoons were added manually

#### Figure S4: Use the data from the previous figure to classify structures by their outcome ####
all_data_summary_final_points %<>% ungroup() %>%
  mutate(Outcome = case_when(
    mean_pct_HET >= 70 ~ 'HET dominant', 
    mean_pct_HMs >= 70 ~ 'HM dominant',
    mean_pct_HET + mean_pct_HMs >= 70 ~ 'Both HM and HET',
    mean_pct_mono >= 70 ~ 'Monomers', 
    TRUE ~ 'Ambiguous'
  ))

all_results %<>% 
  ungroup() %>%
  mutate(Total_complexes = cAA + cAB + cBB + cA + cB) %>% # Add the total concentration of complexes
  mutate(pct_HET = round((cAB *100 / Total_complexes), 2))

# Add the outcome of the replicate
# Classify the final data points accordingly
all_data_final_points <- all_results %>%
  filter(fixed_mut == 200)

all_data_final_points %<>% mutate(HET_pct_total = 100*cAB / (cAA + cAB + cBB + cA + cB)) %>%
  mutate(outcome = ifelse(cAB / Total_complexes >= 0.7, 'HET dominant', 
                          ifelse( (cAA + cBB) / Total_complexes >= 0.7, 'HM dominant', 
                                  ifelse( (cA + cB) / Total_complexes >= 0.7, 'Monomers', 
                                          ifelse( (cAA + cAB + cBB) / Total_complexes >= 0.7, 'Both HM and HET', 
                                                  'Ambiguous')
                                  )
                          )
  )
  ) %>%
  mutate(Complex = factor(toupper(Complex), levels = toupper(all_data_summary_final_points$Complex)),
         outcome = factor(outcome, levels = c('HET dominant', 'Both HM and HET', 'HM dominant', 'Monomers', 'Ambiguous'))
  )

all_results_outcomes <- left_join(x = all_results %>%
                                    mutate(Complex = toupper(Complex)), 
                                  y = all_data_final_points %>% select(Complex, Start_stab, Start_aff, Replicate, outcome), 
                                  by = c('Complex' = 'Complex', 
                                         'Replicate' = 'Replicate',
                                         'Start_stab' = 'Start_stab', 
                                         'Start_aff' = 'Start_aff'))

# Relevel
all_results_outcomes %<>% 
  mutate(Complex = factor(toupper(Complex), levels = toupper(all_data_summary_final_points$Complex)))

p_all_trajectories <- all_results_outcomes %>% 
  filter(Start_stab == -5, Start_aff == -10) %>%
  ggplot(aes(x = fixed_mut, y = pct_HET, group = Replicate, colour = outcome)) +
  xlab('Fixed mutations') + ylab('HET percentage (%)') +
  geom_line() +
  theme(axis.text.y = element_text(size = 24),
        axis.text.x = element_text(size = 24, angle = 90, vjust = 0.5, hjust = 1),
        axis.title = element_text(size = 26),
        legend.title = element_text(size = 26),
        legend.text = element_text(size = 24, hjust = 0.5),
        legend.position = 'top',
        legend.justification = 'center',
        plot.title = element_text(hjust = 0.5),
        strip.text = element_text(size = 28, face = 'bold'), 
        panel.grid.major = element_blank()) +
  facet_wrap(~Complex, nrow = 8) +
  scale_colour_manual(values = c('#4575B4', '#008000', '#FC8D59', 'black', '#b300b3')) +
  labs(colour = 'Outcome') +
  geom_hline(yintercept = 50, linewidth = 1, linetype = 'dashed') +
  guides(colour = guide_legend(override.aes = list(size = c(3, 3, 3, 3, 3))))
p_all_trajectories
ggsave(p_all_trajectories, device = cairo_pdf, width = 45, height = 28, dpi = 500, 
       filename = 'Figures/Supplementary_figures/Appendix_FigS04.pdf')

## Figure S5: No clear association of structural parameters with the outcome ## 

# Load data about interfaces
interface_table = read_delim('Data/Structures/interface_summary.tsv', delim = '\t')

interface_table %<>% mutate(fraction_core = Interface_residues / Total_residues, 
                            fraction_rim = Rim_residues / Total_residues, 
                            fraction_core_rim = (Interface_residues + Rim_residues) / Total_residues)

all_data_summary_final_points_outcome <- all_data_summary_final_points %>% ungroup() %>%
  rowwise() %>%
  mutate(Outcome = ifelse(mean_pct_HET >= 70, 'HET dominant', 
                          ifelse(mean_pct_HMs >= 70, 'HM dominant', 
                                 ifelse(mean_pct_mono >= 70, 'Monomers', 
                                        ifelse((mean_pct_HET + mean_pct_HMs) >= 70, 'Both HM and HET', 
                                               'Ambiguous')
                                 )
                          )
  )
  )

data_final_points_outcome <- left_join(x = all_data_summary_final_points_outcome, 
                                       y = interface_table %>%
                                         mutate(Structure = toupper(Structure)), 
                                       by = c('Complex' = 'Structure'))

# Let's see the distribution of interface sizes in each category
p_figS5A <- data_final_points_outcome %>% 
  ## To remove ambiguous
  filter(Outcome != 'Ambiguous') %>%
  mutate(Outcome = factor(Outcome, levels = c('HET dominant', 'Both HM and HET', 'HM dominant'))) %>%
  ggplot(aes(x = Outcome, y = fraction_core, fill = Outcome)) +
  geom_point(position = position_jitter(width = 0.2), aes(colour = Outcome)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  geom_violin(alpha = 0.3, scale = 'width') +
  stat_compare_means(method = 'wilcox.test', paired = F, size = 3,
                     comparisons = list(c('Both HM and HET', 'HET dominant'), 
                                        c('Both HM and HET', 'HM dominant'), 
                                        c('HET dominant', 'HM dominant'))) +
  scale_fill_manual(values = c('#4575B4', '#008000', '#FC8D59', '#b300b3')) +
  scale_colour_manual(values = c('#4575B4', '#008000', '#FC8D59', '#b300b3')) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8),
        axis.title = element_text(size = 8),
        legend.position =  'none', 
        plot.margin = margin(t = 0.7, l = 0.5, unit = 'cm'), 
        panel.grid.major = element_blank()) +
  xlab('Outcome') + ylab('Fraction of interface core residues') +
  ylim(0, 0.6)
p_figS5A

#### Check if structures with different outcomes have different distributions of stickiness ####

stickiness_interfaces <- read_delim('Data/Structures/interface_stickiness.tsv', delim = '\t')

stickiness_summary <- stickiness_interfaces %>% group_by(PDB, Region) %>%
  summarise(mean_stickiness = mean(Levy_propensity), 
            sem_stickiness = sd(Levy_propensity) / sqrt(n()))

data_final_points_outcome_stickiness <- inner_join(x = all_data_summary_final_points_outcome, 
                                                   y = stickiness_summary %>%
                                                     mutate(PDB = toupper(PDB)), 
                                                   by = c('Complex' = 'PDB'))


p_figS5B <- data_final_points_outcome_stickiness %>% 
  ## To remove ambiguous
  filter(Outcome != 'Ambiguous') %>%
  mutate(Outcome = factor(Outcome, levels = c('HET dominant', 'Both HM and HET', 'HM dominant'))) %>%
  ggplot(aes(x = Outcome, y = mean_stickiness, fill = Outcome)) +
  geom_point(position = position_jitter(width = 0.2), aes(colour = Outcome)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  geom_violin(alpha = 0.3, scale = 'width') +
  scale_colour_manual(values = c('#4575B4', '#008000', '#FC8D59', '#b300b3')) +
  facet_wrap(~Region) +
  stat_compare_means(method = 'wilcox.test', paired = F, size = 3,
                     comparisons = list(c('Both HM and HET', 'HET dominant'), 
                                        c('Both HM and HET', 'HM dominant'), 
                                        c('HET dominant', 'HM dominant'))) +
  scale_fill_manual(values = c('#4575B4', '#008000', '#FC8D59', '#b300b3')) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8), 
        legend.position =  'none', 
        panel.grid.major = element_blank(), 
        axis.title = element_text(size = 8)
        ) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  xlab('Outcome') + ylab('Stickiness') +
  ylim(-0.8, 1.3)
p_figS5B

#### Add a panel for the average distance to self in each structure ####

all_data_summary_outcome <- all_data_summary_final_points %>% ungroup() %>%
  rowwise() %>%
  mutate(Outcome = ifelse(mean_pct_HET >= 70, 'HET dominant',
                          ifelse(mean_pct_HMs >= 70, 'HM dominant',
                                 ifelse(mean_pct_mono >= 70, 'Monomers',
                                        ifelse((mean_pct_HET + mean_pct_HMs) >= 70, 'Both HM and HET',
                                               'Ambiguous')
                                 )
                          )
  )
  )

## Load data on distance to self ##
self_distance <- read_delim('Data/Structures/interface_distance_self.tsv', delim = '\t')

self_distance_summ <- self_distance %>% group_by(PDB, Region) %>%
  mutate(self_contact = (Min_dist <= 4)) %>%
  summarise(mean_distance = mean(Min_dist), 
            sem_distance = sd(Min_dist) / sqrt(n()), 
            prop_self_contact = mean(self_contact)) %>%
  arrange(PDB, Region)

## Join the data
distance_means <- left_join(x = self_distance_summ %>%
                              mutate(PDB = toupper(PDB)),
                            y = all_data_summary_final_points %>% select(Complex, med_pct_HET, mean_pct_HET, sem_pct_HET), 
                            by = c('PDB' = 'Complex'))

distance_means_final <- inner_join(x = distance_means, 
                                   y = all_data_summary_outcome %>% select(Complex, Outcome), 
                                   by = c('PDB' = 'Complex'))

## Try the proportion of residues in self contact
p_figS5C <- distance_means_final %>%
  mutate(Region = case_when(
    Region == 0.75 ~ 'Rim', 
    Region == 1 ~ 'Core'
  ), 
  Outcome = factor(Outcome, levels = c('HET dominant', 'Both HM and HET', 'HM dominant'))) %>%
  filter(Region == 'Core') %>%
  ggplot(aes(x = Outcome, y = prop_self_contact, fill = Outcome)) +
  geom_jitter(width = 0.2, aes(colour = Outcome)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  geom_violin(alpha = 0.3, scale = 'width') +
  # facet_wrap(~Region) +
  stat_compare_means(method = 'wilcox.test', paired = F, size = 3,
                     comparisons = list(c('Both HM and HET', 'HET dominant'), 
                                        c('Both HM and HET', 'HM dominant'), 
                                        c('HET dominant', 'HM dominant'))) +
  scale_fill_manual(values = c('#4575B4', '#008000', '#FC8D59', '#b300b3')) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8), 
        legend.position =  'none', 
        plot.margin = margin(t = 0.7, l = 0.5, unit = 'cm'), 
        panel.grid.major = element_blank(), 
        axis.title = element_text(size = 8)
  ) +
  scale_colour_manual(values = c('#4575B4', '#008000', '#FC8D59', '#b300b3')) +
  xlab('Outcome') + ylab('Proportion of self-contacting residues') +
  ylim(0, 0.6)
p_figS5C

#### New panel on secondary structure content ####
## Load the DSSP data
dssp_data <- read_delim('Data/Structures/summary_table_DSSP_all_structures.txt', delim = '\t')
dssp_data %<>% mutate(total_residues = none + Missing + `Alpha helix` + `Beta bridge` + `Beta ladder` + `3/10 helix` + `5-helix` + `H-bonded turn` + Bend) %>%
  mutate(prop_alpha = `Alpha helix` / total_residues, 
         prop_beta = (`Beta bridge` + `Beta ladder`) / total_residues)

dssp_data_long <- dssp_data %>% select(Structure, prop_alpha, prop_beta) %>%
  pivot_longer(cols = c('prop_alpha', 'prop_beta'), names_to = 'sec_struc', values_to = 'proportion')

## Add the data for the outcome
dssp_data_final <- inner_join(x = dssp_data_long %>%
                                mutate(Structure = toupper(Structure)), 
                              y = all_data_summary_outcome %>% select(Complex, Outcome), 
                              by = c('Structure' = 'Complex'))

p_figS5D <- dssp_data_final %>%
  mutate(Outcome = factor(Outcome, levels = c('HET dominant', 'Both HM and HET', 'HM dominant')), 
         sec_struc = case_when(
           sec_struc == 'prop_beta' ~ 'Beta strands',
           sec_struc == 'prop_alpha' ~ 'Alpha helices'
         )) %>%
  ggplot(aes(x = Outcome, y = proportion, fill = Outcome)) +
  geom_jitter(width = 0.2, aes(colour = Outcome)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  geom_violin(alpha = 0.3, scale = 'width') +
  facet_wrap(~sec_struc) +
  stat_compare_means(method = 'wilcox.test', paired = F, size = 3,
                     comparisons = list(c('Both HM and HET', 'HET dominant'), 
                                        c('Both HM and HET', 'HM dominant'), 
                                        c('HET dominant', 'HM dominant'))) +
  scale_fill_manual(values = c('#4575B4', '#008000', '#FC8D59', '#b300b3')) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8), 
        legend.position =  'none', 
        panel.grid.major = element_blank(), 
        axis.title = element_text(size = 8)
        ) +
  scale_colour_manual(values = c('#4575B4', '#008000', '#FC8D59', '#b300b3')) +
  xlab('Outcome') + ylab('Secondary structure fraction') +
  # ylim(0, 125)
  ylim(0, 1.25)
p_figS5D

## Read the ECOD annotations
## Keep in mind that ECOD data are organized by hierarchical classifications (top to bottom):
# Architecture
# Possible homology (X)
# Homology (H) 
# Topology (T)
# Family (F).

# Read the data from ecod
ecod_data <- read_delim('Data/Structures/ecod.latest.domains.txt', delim = '\t', skip = 4)

# Use an inner join with the data from the simulations
all_data_summary_ecod <- inner_join(x = all_data_summary_final_points, 
                                    y = ecod_data %>% mutate(pdb = toupper(pdb)),
                                    by = c('Complex' = 'pdb'))

## Get the median for each architecture type to show them by decreasing HET dominance
all_data_summary_ecod_medians <- all_data_summary_ecod %>% filter(chain == 'A') %>%
  group_by(arch_name) %>% summarise(median_HET_pct = median(mean_pct_HET)) %>%
  arrange(desc(median_HET_pct))

## Check how the HET_pct varies for different ECOD groups
p_figS5E <- all_data_summary_ecod %>% 
  filter(chain == 'A') %>%
  mutate(arch_name = factor(arch_name, levels = all_data_summary_ecod_medians$arch_name)) %>% rowwise() %>%
  mutate(arch_type = ifelse(str_detect(pattern = 'alpha', string = arch_name), 'Alpha', 
                            ifelse(str_detect(pattern = 'beta', string = arch_name), 'Beta', 
                                   ifelse(or(str_detect(pattern = 'a/b', string = arch_name), str_detect(pattern = 'a\\+b', string = arch_name)), 'Mixed',
                                          'Other')))) %>%
  mutate(arch_type = factor(arch_type, levels = c('Alpha', 'Beta', 'Mixed', 'Other'))) %>%
  ggplot(aes(y = arch_name, x = mean_pct_HET, fill = arch_type)) +
  scale_fill_manual(values = c('#2c8221', '#0a2f35', '#f7a325', 'gray')) +
  geom_rect(inherit.aes = F, fill = 'grey', colour = NA, xmin = -Inf, xmax = Inf, ymin = 30, ymax = 70, alpha = 0.3) +
  geom_jitter(height = 0.2, aes(colour  = arch_type)) +
  scale_colour_manual(values = c('#2c8221', '#0a2f35', '#f7a325', 'gray')) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_vline(xintercept = 50, linetype = 'dashed') +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.position = 'top', legend.justification = 'center', 
        axis.title = element_text(size = 8)
        # plot.margin = margin(l = 1, r = 0.5, t = 0.5, b = 0.5, unit = 'cm')
        ) +
  ylab('ECOD architecture') + xlab('HET percentage (%)') +
  labs(fill = '', colour = '')
p_figS5E
p_figS5_top <- plot_grid(p_figS5A, p_figS5B, p_figS5C, p_figS5D, 
                         nrow = 2, labels = c('A', 'B', 'C', 'D'), 
                         label_fontface = 'bold', label_size = panel_label_size)

p_figS5 <- plot_grid(p_figS5_top, p_figS5E, nrow = 2, labels = c('', 'E'), 
                     label_fontface = 'bold', label_size = panel_label_size, 
                     rel_heights = c(2, 1))

ggsave(p_figS5, device = cairo_pdf, width = 24, height = 24, units = 'cm', dpi = 300, 
       filename = 'Figures/Supplementary_figures/Appendix_FigS05.structural_properties.pdf')

#### Save the data for the structures (dataset EV1) ####
## Need to put together the data for the different items
ecod_interface_core <- inner_join(
  x = all_data_summary_ecod %>% ungroup() %>%
    mutate(Complex = as.character(Complex)) %>%
    select(Complex, mean_pct_HET, mean_pct_HMs, ecod_domain_id,
           pdb_range, seqid_range, arch_name, x_name, h_name, t_name, f_name), 
  y = data_final_points_outcome %>% ungroup() %>%
    mutate(Total_residue_count = Total_residues, 
           Core_interface_residue_count = Interface_residues, 
           Rim_residue_count = Rim_residues, 
           Fraction_core = fraction_core, 
           Fraction_rim = fraction_rim, 
           Complex = as.character(Complex)
           ) %>%
    select(Complex, Total_residue_count, Core_interface_residue_count, Rim_residue_count, 
           Fraction_core, Fraction_rim), 
  by = c('Complex' = 'Complex'))

## Add the stickiness of interface cores and rims
ecod_interface_stickiness <- left_join(x = ecod_interface_core, 
                                       y = data_final_points_outcome_stickiness %>%
                                         select(Complex, Region, mean_stickiness) %>%
                                         pivot_wider(names_from = Region, values_from = mean_stickiness, 
                                                     names_prefix = 'mean_stickiness_'),
                                       by = c('Complex' = 'Complex')
                                       )

ecod_interface_stickiness_self_contacts <- left_join(x = ecod_interface_stickiness,
                                                     y = distance_means_final %>% filter(Region == 1) %>%
                                                       select(PDB, prop_self_contact), 
                                                     by = c('Complex' = 'PDB')
                                                     )

dataset_EV1 <- left_join(x = ecod_interface_stickiness_self_contacts, 
                      y = dssp_data_final %>% select(Structure, sec_struc, proportion) %>%
                        pivot_wider(names_from = sec_struc, values_from = proportion, 
                                    names_prefix = 'proportion_'), 
                      by = c('Complex' = 'Structure'))

write.table(dataset_EV1, append = F, quote = F, sep = '\t', 
            row.names = F, col.names = T, 
            file = 'dataset_EV1.tsv')


#### Figure EV2: Similar to Figure 2C but with alpha = 60, correlation with alpha =80 ####

## New data with pdup = 1, regular simulations
all_results_60opt <- read_delim('Results_simulations/Results_pdup1/Simulations_with_structures/008_simulations_60opt/all_results_all_sims.tsv', delim = '\t')

## Remove some unneeded columns
all_results_60opt %<>% 
  select(-sA, -dA, -aA, -dAA, -aAA, -sB, -dB, -aB, -dBB, -aBB, -dAB, -aAB, -t)

# Add a column indicating how many mutations have fixed up to that point
all_results_60opt %<>% group_by(Complex, Replicate, Start_stab, Start_aff) %>%
  mutate(fixed_mut = row_number()) %>%
  filter(!(Complex %in% structures_remove))

## Look at the final points to check the outcomes of simulations
all_data_final_points_60opt <- all_results_60opt %>% ungroup() %>%
  mutate(Total_complexes = cAA + cAB + cBB + cA + cB) %>% # Add the total concentration of complexes
  filter(fixed_mut == 200, Start_stab == -5, Start_aff == -10) %>%
  mutate(outcome = ifelse(cAB / Total_complexes >= 0.7, 'HET dominant', 
                          ifelse( (cAA + cBB) / Total_complexes >= 0.7, 'HM dominant', 
                                  ifelse( (cA + cB) / Total_complexes >= 0.7, 'Monomers', 
                                          ifelse( (cAA + cAB + cBB) / Total_complexes >= 0.7, 'Both HM and HET', 
                                                  'Ambiguous')
                                  )
                          )
  )
  )

all_data_summary_60opt <- all_results_60opt %>% 
  ungroup() %>%
  mutate(Total_complexes = cAA + cAB + cBB + cA + cB) %>% # Add the total concentration of complexes
  mutate(pct_HET = round((cAB *100 / Total_complexes), 2), 
         pct_HMs = round((cAA + cBB) * 100 / Total_complexes, 2), 
         pct_mono = round((cA + cB) * 100 / Total_complexes, 2)) %>%
  group_by(Start_stab, Start_aff, Complex, fixed_mut) %>%
  summarise(
    med_pct_HET = median(pct_HET), 
    mean_pct_HET = mean(pct_HET), 
    sem_pct_HET = sd(pct_HET) / sqrt(n()), 
    
    mean_pct_HMs = mean(pct_HMs), 
    sem_pct_HMs = sd(pct_HMs) / sqrt(n()),
    
    mean_pct_mono = mean(pct_mono), 
    sem_pct_mono = sd(pct_mono) / sqrt(n())
  )

# Relevel according to the order from the other plot
all_data_summary_60opt %<>% 
  mutate(Complex = factor(toupper(Complex), levels = toupper(all_data_summary_final_points$Complex)))

# Get the order of the figures
all_data_summary_final_points_60opt <- all_data_summary_60opt %>% 
  filter(fixed_mut == 200, Start_stab == -5, Start_aff == -10) %>%
  arrange(desc(mean_pct_HET))

all_data_summary_final_points_60opt %<>% arrange(desc(mean_pct_HET))

all_data_summary_final_points_60opt %<>%
  mutate(Complex = factor(toupper(Complex), levels = toupper(all_data_summary_final_points$Complex)))

all_data_summary_final_points_long_60opt <- all_data_summary_final_points_60opt %>% ungroup() %>%
  pivot_longer(cols = c(mean_pct_HET, mean_pct_HMs, mean_pct_mono), 
               names_to = 'assembly', values_to = 'mean_value')

# Stacked barplot of outcomes per structure
p_figEV2A <- all_data_summary_final_points_long_60opt %>% ungroup() %>%
  mutate(assembly_label = ifelse(assembly == 'mean_pct_HET', 'HET', 
                                 ifelse(assembly == 'mean_pct_HMs', 'HMs (AA+BB)',
                                        'Monomers (A+B)')), 
         Complex = toupper(Complex)) %>%
  mutate(assembly_label = factor(assembly_label, levels = c('HET', 'Monomers (A+B)', 'HMs (AA+BB)')),
         Complex = factor(Complex, levels = toupper(all_data_summary_final_points$Complex))
  ) %>%
  ggplot(aes(x = Complex, y = mean_value, fill = assembly_label)) +
  geom_bar(stat = 'identity') +
  xlab('PDB ID') + ylab('Percentage (%)') +
  labs(fill = '') +
  geom_hline(yintercept = 50, linetype = 'dashed') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.position = 'top', 
        legend.justification = 'center') +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  scale_fill_manual(values = c('#4575B4', 'grey', '#FC8D59' ))
p_figEV2A

#### Figure EV2B: The correlation between the results with alpha 60 and alpha 80 ####

corr_check_diff_opt <- inner_join(x = all_data_summary_final_points_60opt %>% ungroup() %>%
                                    mutate(Complex = toupper(Complex)) %>%
                                    select('Complex', 'mean_pct_HET', 'mean_pct_HMs', 'mean_pct_mono'), 
                                  y = all_data_summary_final_points %>% ungroup() %>%
                                    select('Complex', 'mean_pct_HET', 'mean_pct_HMs', 'mean_pct_mono'), 
                                  by = c('Complex' = 'Complex')
)

p_figEV2B <- corr_check_diff_opt %>%
  ggplot(aes(x = mean_pct_HET.x, y = mean_pct_HET.y)) +
  geom_point() +
  xlab('HET percentage (pdup = 1, alpha = 60)') +
  ylab('HET percentage (pdup = 1, alpha = 80)') +
  theme(panel.grid.major = element_blank(), 
        axis.title = element_text(size = 9)) +
  geom_hline(yintercept =50, linetype = 'dashed') +
  geom_vline(xintercept = 50, linetype = 'dashed') +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', col = 'red') +
  xlim(0, 100) + ylim(0, 100) +
  stat_cor(method = 'spearman', p.accuracy = 0.001, r.accuracy = 0.01, cor.coef.name = 'rho', size = 3,
           label.y.npc = 0.7, label.x.npc = 0.05) +
  stat_cor(method = 'pearson', p.accuracy = 0.001, r.accuracy = 0.01, cor.coef.name = 'r', size = 3,
           label.y.npc = 0.8, label.x.npc = 0.05) +
  geom_smooth(method = 'lm', alpha = 0.5)
p_figEV2B

#### Figure EV2C: Decrease in folding energy in simulations with alpha = 60 ####

all_results_stab_long_60opt <- all_results_60opt %>% ungroup() %>%
  filter(Start_stab == -5, Start_aff == -10) %>% rowwise() %>%
  select(Replicate, Complex, fixed_mut, curr_stab_A, curr_stab_B) %>%
  mutate(less_stable_subunit = max(curr_stab_A, curr_stab_B), 
         more_stable_subunit = min(curr_stab_A, curr_stab_B)) %>%
  pivot_longer(cols = c('less_stable_subunit', 'more_stable_subunit'), names_to = 'chain',
               values_to = 'curr_stab') %>%
  mutate(chain = case_when(
    chain == 'less_stable_subunit' ~ 'Less stable subunit',
    chain == 'more_stable_subunit' ~ 'More stable subunit'
  ))

## Draw the figure
p_figEV2C <- all_results_stab_long_60opt %>% 
  filter(Complex == '1gpr') %>%
  ggplot(aes(x = fixed_mut, y = curr_stab, colour = chain)) +
  geom_point(alpha = 0.1) +
  facet_wrap(~chain) +
  geom_line(aes(group = interaction(Replicate, Complex, chain))) +
  theme(legend.position = 'none', panel.grid.major = element_blank(), 
        plot.title = element_text(hjust = 0.5, face = 'bold', size = 14),
        strip.text = element_text(size = 12)
        ) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  ylim(-8, 2) +
  xlab('Fixed mutation') +
  labs(y = expression(
    paste(
      bold('\u0394'), 
      bold(G['fold']), 
      sep = ''
    )
  )) +
  ggtitle(expression(paste(bold('\u0394'),
                           bold(G[fold]), 
                           bold(' in simulations with alpha = 60'),
                           sep = ''
  )
  )
  )
p_figEV2C

#### Figure EV2D: No decrease in folding energy in simulations with alpha = 80 ####

## Distinguish the least stable of the two chains from the more stable
all_results_stab_long <- all_results %>% ungroup() %>%
  filter(Start_stab == -5, Start_aff == -10) %>% rowwise() %>%
  select(Replicate, Complex, fixed_mut, curr_stab_A, curr_stab_B) %>%
  mutate(less_stable_subunit = max(curr_stab_A, curr_stab_B), 
         more_stable_subunit = min(curr_stab_A, curr_stab_B)) %>%
  pivot_longer(cols = c('less_stable_subunit', 'more_stable_subunit'), names_to = 'chain',
               values_to = 'curr_stab') %>%
  mutate(chain = case_when(
    chain == 'less_stable_subunit' ~ 'Less stable subunit',
    chain == 'more_stable_subunit' ~ 'More stable subunit'
  ))

## Draw the figure
p_figEV2D <- all_results_stab_long %>% 
  filter(Complex == '1gpr') %>%
  ggplot(aes(x = fixed_mut, y = curr_stab, colour = chain)) +
  geom_point(alpha = 0.1) +
  facet_wrap(~chain) +
  geom_line(aes(group = interaction(Replicate, Complex, chain))) +
  theme(legend.position = 'none', panel.grid.major = element_blank(), 
        plot.title = element_text(hjust = 0.5, size = 14, face = 'bold'), 
        strip.text = element_text(size = 12)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  ylim(-8, 2) +
  xlab('Fixed mutation') +
  labs(y = expression(
    paste(
      bold('\u0394'), 
      bold(G['fold']), 
      sep = ''
    )
  )) +
  ggtitle(expression(paste(bold('\u0394'),
                           bold(G[fold]), 
                           bold(' in simulations with alpha = 80'),
                           sep = ''
  )
  )
  )

p_figEV2D

#### Figure EV2E: Transient decrease in HET percentage in the simulations with alpha = 60 ####

## Get the list of fixed mutational effects for each structure
all_results_fixed_effects_60opt <- all_results_60opt %>% ungroup() %>%
  mutate(pct_HET = cAB * 100 / (cA + cB + cAA + cBB + cAB)) %>%
  mutate(mut_eff_stab_A = curr_stab_A - lag(curr_stab_A, default = NA),
         mut_eff_stab_B = curr_stab_B - lag(curr_stab_B, default = NA),
         mut_eff_binding_energy_AA = curr_binding_energy_AA - lag(curr_binding_energy_AA, default = NA),
         mut_eff_binding_energy_AB = curr_binding_energy_AB - lag(curr_binding_energy_AB, default = NA),
         mut_eff_binding_energy_BB = curr_binding_energy_BB - lag(curr_binding_energy_BB, default = NA),
         delta_cAB = cAB - lag(cAB, default = NA),
         delta_pct_HET = pct_HET - lag(pct_HET, default = NA)
  ) %>%
  group_by(Replicate, Start_stab, Start_aff, Complex) %>%
  filter(fixed_mut >= 3) # The first line is the starting point, the second one is the duplication 

values_check <- c(5, 10, 25, 50, 100, 198)
all_plots_60opt <- list()
all_het_pct_60opt <- c()

for(i in 1:length(values_check)){
  n_mut_check <- values_check[i]
  
  # Get the net effect on HM and HET binding
  all_results_fixed_effects_tmp <- all_results_fixed_effects_60opt %>% ungroup() %>%
    filter(Start_stab == -5, Start_aff == -10) %>%
    mutate(netEffect_HETbind = mut_eff_binding_energy_AB, 
           # Mutations can only affect one of the two HMs at a time, so the registered change
           # has to be zero for one of them
           netEffect_HMbind = mut_eff_binding_energy_AA + mut_eff_binding_energy_BB, 
           total_complexes = cA + cB + cAA + cAB + cBB) %>%
    mutate(pct_HET = (cAB * 100) / total_complexes) %>%
    select(Replicate, Complex,
           netEffect_HETbind, netEffect_HMbind, pct_HET) %>% 
    group_by(Replicate, Complex) %>%
    # Calculate residuals
    mutate(linear_residuals = netEffect_HETbind - (0.5 * netEffect_HMbind), 
           mut_num = row_number())
  
  # Get the percentage of HET after the mutations that accumulated
  current_pct_HET <- all_results_fixed_effects_tmp %>%
    filter(mut_num == n_mut_check)
  
  # Get the first n mutations
  all_results_fixed_effects_summary <- all_results_fixed_effects_tmp %>%
    group_by(Replicate, Complex) %>%
    filter(mut_num <= n_mut_check) %>%
    summarise(fixed_residual_sum = sum(linear_residuals))
  
  
  all_results_fixed_effects_summary <- inner_join(x = all_results_fixed_effects_summary,
                                                  y = current_pct_HET %>%
                                                    select(Replicate, Complex, pct_HET),
                                                  by = c('Replicate' = 'Replicate',
                                                         'Complex' = 'Complex'))
  
  ## Show in a figure the relation between the first fixed residuals and the outcome
  p <- all_results_fixed_effects_summary %>%
    ggplot(aes(x = fixed_residual_sum, y = pct_HET)) +
    geom_point(alpha = 0.2) +
    xlim(-7, 7) +
    ggtitle(str_c('First ', toString(n_mut_check), ' mutations', sep = '')) +
    xlab('Cumulative residual sum') + ylab('HET percentage (%)')
  
  ## Save the figure to the list
  all_plots_60opt[[i]] <- p
  
  ## Accumulate all the values for pct_HET to show the distributions
  all_het_pct_60opt <- bind_rows(all_het_pct_60opt, 
                                 current_pct_HET %>% select(Replicate, Complex, pct_HET, mut_num))
  
}

## Draw boxplots of the distributions of percentages of HET
p_figEV2E <- all_het_pct_60opt %>% 
  ggplot(aes(x = as.factor(mut_num), y = pct_HET)) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  geom_violin(alpha = 0.3, scale = 'width') +
  theme(panel.grid.major = element_blank(), 
        plot.title = element_text(face = 'bold', size = 14, hjust = 0.5)) +
  xlab('Fixed mutations') + ylab('HET percentage (%)') +
  geom_hline(yintercept = 50, linetype = 'dashed') +
  ggtitle('HET percentage in simulations with alpha = 60')
p_figEV2E


#### Figure EV2F: No transient decrease in HET percentage in simulations with alpha = 80 ####

## Get the list of fixed mutational effects for each structure
all_results_fixed_effects <- all_results %>% ungroup() %>%
  mutate(pct_HET = cAB * 100 / (cA + cB + cAA + cBB + cAB)) %>%
  mutate(mut_eff_stab_A = curr_stab_A - lag(curr_stab_A, default = NA),
         mut_eff_stab_B = curr_stab_B - lag(curr_stab_B, default = NA),
         mut_eff_binding_energy_AA = curr_binding_energy_AA - lag(curr_binding_energy_AA, default = NA),
         mut_eff_binding_energy_AB = curr_binding_energy_AB - lag(curr_binding_energy_AB, default = NA),
         mut_eff_binding_energy_BB = curr_binding_energy_BB - lag(curr_binding_energy_BB, default = NA),
         delta_cAB = cAB - lag(cAB, default = NA),
         delta_pct_HET = pct_HET - lag(pct_HET, default = NA)
  ) %>%
  group_by(Replicate, Start_stab, Start_aff, Complex) %>%
  filter(fixed_mut >= 3) # The first line is the starting point, the second one is the duplication 

values_check <- c(5, 10, 25, 50, 100, 198)
all_plots <- list()
all_het_pct <- c()

for(i in 1:length(values_check)){
  n_mut_check <- values_check[i]
  
  # Get the net effect on HM and HET binding
  all_results_fixed_effects_tmp <- all_results_fixed_effects %>% ungroup() %>%
    filter(Start_stab == -5, Start_aff == -10) %>%
    mutate(netEffect_HETbind = mut_eff_binding_energy_AB, 
           # Mutations can only affect one of the two HMs at a time, so the registered change
           # has to be zero for one of them
           netEffect_HMbind = mut_eff_binding_energy_AA + mut_eff_binding_energy_BB, 
           total_complexes = cA + cB + cAA + cAB + cBB) %>%
    mutate(pct_HET = (cAB * 100) / total_complexes) %>%
    select(Replicate, Complex,
           netEffect_HETbind, netEffect_HMbind, pct_HET) %>% 
    group_by(Replicate, Complex) %>%
    # Calculate residuals
    mutate(linear_residuals = netEffect_HETbind - (0.5 * netEffect_HMbind), 
           mut_num = row_number())
  
  # Get the percentage of HET after the mutations that accumulated
  current_pct_HET <- all_results_fixed_effects_tmp %>%
    filter(mut_num == n_mut_check)
  
  # Get the first n mutations
  all_results_fixed_effects_summary <- all_results_fixed_effects_tmp %>%
    group_by(Replicate, Complex) %>%
    filter(mut_num <= n_mut_check) %>%
    summarise(fixed_residual_sum = sum(linear_residuals))
  
  
  all_results_fixed_effects_summary <- inner_join(x = all_results_fixed_effects_summary,
                                                  y = current_pct_HET %>%
                                                    select(Replicate, Complex, pct_HET),
                                                  by = c('Replicate' = 'Replicate',
                                                         'Complex' = 'Complex'))
  
  ## Show in a figure the relation between the first fixed residuals and the outcome
  p <- all_results_fixed_effects_summary %>%
    ggplot(aes(x = fixed_residual_sum, y = pct_HET)) +
    geom_point(alpha = 0.2) +
    xlim(-7, 7) +
    ggtitle(str_c('First ', toString(n_mut_check), ' mutations', sep = '')) +
    xlab('Cumulative residual sum') + ylab('HET percentage (%)')
  
  ## Save the figure to the list
  all_plots[[i]] <- p
  
  ## Accumulate all the values for pct_HET to show the distributions
  all_het_pct <- bind_rows(all_het_pct, 
                           current_pct_HET %>% select(Replicate, Complex, pct_HET, mut_num))
  
}

## Draw boxplots of the distributions of percentages of HET
p_figEV2F <- all_het_pct %>% 
  ggplot(aes(x = as.factor(mut_num), y = pct_HET)) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  geom_violin(alpha = 0.3, scale = 'width') +
  theme(panel.grid.major = element_blank(), 
        plot.title = element_text(face = 'bold', size = 14, hjust = 0.5)) +
  xlab('Fixed mutations') + ylab('HET percentage (%)') +
  geom_hline(yintercept = 50, linetype = 'dashed') +
  ggtitle('HET percentage in simulations with alpha = 80')
p_figEV2F

## Draw figure S5 ##

p_figEV2_top <- plot_grid(p_figEV2A, p_figEV2B, rel_widths = c(3, 1.2), ncol = 2, 
                         labels = c('A', 'B'), label_fontface = 'bold', label_size = panel_label_size)

p_figEV2_rest <- plot_grid(p_figEV2C, p_figEV2D, p_figEV2E, p_figEV2F, 
                          nrow = 2, labels = c('C', 'D', 'E', 'F'), label_fontface = 'bold',
                          label_size = panel_label_size)

p_figEV2 <- plot_grid(p_figEV2_top, p_figEV2_rest, nrow = 2, rel_heights = c(1, 2))

ggsave(p_figEV2, device = cairo_pdf, width = 38, height = 24, units = 'cm', dpi = 300, 
       filename = 'Figures/ExpandedView_figures/FigEV2_new.pdf')

#### Figure 3 ####

### Factors that influence the outcome of simulations for each structure
## Load the mutational effects
mut_effects <- read_delim('Data/all_mut_matrices.tsv', '\t')

## Keep only the structures that are still in the results matrix
mut_effects %<>% filter(Complex %in% all_results$Complex)

## Let's group by structure
means_effects <- mut_effects %>% ungroup() %>%
  group_by(Complex) %>%
  summarise(Mean_ddG_stab_HET = mean(Mean_ddG_stab_HET), 
            Mean_ddG_int_HM = mean(Mean_ddG_int_HM), 
            Mean_ddG_int_HET = mean(Mean_ddG_int_HET), 
            mean_dev_diagonal = mean(Mean_ddG_int_HET - 0.5*Mean_ddG_int_HM), 
            median_dev_diagonal = median(Mean_ddG_int_HET - 0.5*Mean_ddG_int_HM))

#### Check if the means of the distributions of deviations from the diagonal relate to the outcome ####

## Join the previous table with the outcomes from the previous figures
mut_effects_means <- inner_join(x = means_effects,
                                y = all_data_summary_final_points %>% select(Complex, med_pct_HET, mean_pct_HET, sem_pct_HET), 
                                by = c('Complex' = 'Complex'))

## Figure 3A: Distribution of available effects for the two most extreme cases ##
#### Work with the residuals ####
compare_outliers <- function(pdb1, pdb2){
  pdb1 <- toupper(pdb1)
  pdb2 <- toupper(pdb2)
  
  mut_effects_extreme <- mut_effects %>%
    mutate(Complex = toupper(Complex)) %>%
    filter(Complex %in% c(pdb1, pdb2), Mean_ddG_int_HET != 0, Mean_ddG_int_HM != 0) %>%
    mutate(Complex = factor(Complex, levels = c(pdb1, pdb2)))
  
  p_center <- mut_effects_extreme %>% 
    mutate(Complex = ifelse(Complex == pdb1, str_c(pdb1, ' (HET dominant)', sep = ''),
                            str_c(pdb2, ' (HM dominant)', sep = ''))) %>%
    mutate(Complex = factor(Complex, levels = c(
      str_c(pdb1, ' (HET dominant)', sep = ''), 
      str_c(pdb2, ' (HM dominant)', sep = '')
    ))) %>%
    ggplot(aes(x = Mean_ddG_int_HM, y = Mean_ddG_int_HET, colour = Complex)) +
    geom_point(alpha = 0.5) +
    scale_colour_manual(values = c('#4575B4', '#FC8D59')) +
    xlim(-3, 5) + ylim(-3, 5) +
    geom_vline(xintercept = 0, linetype = 'dashed') +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    geom_abline(slope = 0.5, intercept = 0, linetype = 'dashed') +
    theme(legend.position = 'none', 
          panel.grid.major = element_blank(), 
          plot.title = element_text(hjust = 0.5), 
          axis.title = element_text(size = 9)) +
    labs(x  = expression(paste(bold('\u0394\u0394'), 
                               bold(G['bind,HM']), 
                               bold(' (kcal/mol)'))), 
         y = expression(paste(bold('\u0394\u0394'), 
                              bold(G['bind,HET']), 
                              bold(' (kcal/mol)'))))

  p_center_new <- p_center + 
    theme(legend.position = 'top',
          legend.justification = 'center', 
          legend.spacing.x = unit(0.12, 'cm')) + 
    labs(colour = 'Complex', fill = 'Complex')
  
  p_legend <- get_legend(p_center_new)
  
  p_density_right <- mut_effects_extreme %>%
    ggplot(aes(y = Mean_ddG_int_HET, fill = Complex)) +
    scale_fill_manual(values = c('#4575B4', '#FC8D59')) +
    geom_density(alpha = 0.5
    ) +
    ylim(-3, 5) +
    theme(legend.position = 'none', 
          panel.grid.major = element_blank(),
          panel.border = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          plot.margin = margin(l = -0.7, r = 0, b = 0.7, unit = 'cm')
    ) +
    xlab('') + ylab('')
  
  
  p_density_top<- mut_effects_extreme %>%
    ggplot(aes(x = Mean_ddG_int_HM, fill = Complex)) +
    geom_density(alpha = 0.5
    ) +
    scale_fill_manual(values = c('#4575B4', '#FC8D59')) +
    xlim(-3, 5) +
    theme(legend.position = 'none', 
          panel.grid.major = element_blank(),
          panel.border = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          plot.margin = margin(t = 0, b = -0.5, l = 0.5, unit = 'cm')
    ) +
    xlab('') + ylab('')
  
  p_distributions_comp <- plot_grid(p_density_top, NULL, 
                                    p_center, p_density_right, ncol = 2, nrow = 2, 
                                    rel_widths = c(1, 0.15), rel_heights = c(0.15, 1))

  
  return(p_distributions_comp)
}

## Most extreme comparison with alpha = 80
p_fig3A <- compare_outliers('2b18', '3ulh')
p_fig3A

## Figure 3B: Distribution of fixed effects for the two most extreme cases ##
## Get the list of fixed mutational effects for each structure
all_results_fixed_effects <- all_results %>% ungroup() %>%
  mutate(pct_HET = cAB * 100 / (cA + cB + cAA + cBB + cAB)) %>%
  mutate(mut_eff_stab_A = curr_stab_A - lag(curr_stab_A, default = NA),
         mut_eff_stab_B = curr_stab_B - lag(curr_stab_B, default = NA),
         mut_eff_binding_energy_AA = curr_binding_energy_AA - lag(curr_binding_energy_AA, default = NA),
         mut_eff_binding_energy_AB = curr_binding_energy_AB - lag(curr_binding_energy_AB, default = NA),
         mut_eff_binding_energy_BB = curr_binding_energy_BB - lag(curr_binding_energy_BB, default = NA),
         delta_cAB = cAB - lag(cAB, default = NA),
         delta_pct_HET = pct_HET - lag(pct_HET, default = NA)
  ) %>%
  group_by(Replicate, Start_stab, Start_aff, Complex) %>%
  filter(fixed_mut >= 3) # The first line is the starting point, the second one is the duplication 

all_results_fixed_effects %<>% ungroup() %>%
  filter(mut_eff_binding_energy_AB != 0) %>%
  mutate(mut_eff_HM = case_when(
    abs(mut_eff_binding_energy_AA) >= abs(mut_eff_binding_energy_BB) ~ mut_eff_binding_energy_AA, 
    abs(mut_eff_binding_energy_BB) > abs(mut_eff_binding_energy_AA) ~ mut_eff_binding_energy_BB
  )) %>%
  mutate(diff_eff_HET_HM = mut_eff_binding_energy_AB - mut_eff_HM)

p_center <- all_results_fixed_effects %>%
  filter(Complex %in% c('2b18', '3ulh')) %>%
  mutate(Complex = case_when(
    Complex == '2b18' ~ '2B18 (HET dominant)', 
    Complex == '3ulh' ~ '3ULH (HM dominant)'
  )) %>%
  ggplot(aes(x = mut_eff_HM, y = mut_eff_binding_energy_AB, colour = Complex)) +
  geom_point(alpha = 0.3) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_abline(slope = 0.5, intercept = 0, linetype = 'dashed') +
  scale_colour_manual(values = c('#4575B4', '#FC8D59')) +
  theme(panel.grid.major = element_blank(), 
        legend.position = 'none', legend.justification = 'center', 
        plot.title = element_text(hjust = 0.5, size = 14), 
        axis.title = element_text(size = 9)
  ) +
  labs(
    x = expression(paste(bold('\u0394\u0394'), 
                         bold(G['bind,HM']), 
                         bold(' (kcal/mol)'), sep = '')), 
    y = expression(paste(bold('\u0394\u0394'), 
                         bold(G['bind,HET']), 
                         bold(' (kcal/mol)'), sep = ''))
  ) +
  xlim(-3, 5) + ylim(-3, 5)
p_center

p_center_new <- p_center + 
  theme(legend.position = 'top',
        legend.justification = 'center', 
        legend.spacing.x = unit(0.12, 'cm')) + 
  labs(colour = 'Complex', fill = 'Complex') +
  guides(colour = guide_legend(override.aes = list(alpha = 1)))


p_legend <- get_legend(p_center_new)

p_density_right <- all_results_fixed_effects %>%
  filter(Complex %in% c('2b18', '3ulh')) %>%
  mutate(Complex = case_when(
    Complex == '2b18' ~ '2B18 (HET dominant)', 
    Complex == '3ulh' ~ '3ULH (HM dominant)'
  )) %>%
  ggplot(aes(y = mut_eff_binding_energy_AB, fill = Complex)) +
  scale_fill_manual(values = c('#4575B4', '#FC8D59')) +
  geom_density(alpha = 0.5# , scale = 'width'
  ) +
  ylim(-3, 5) +
  theme(legend.position = 'none', 
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = margin(l = -0.7, r = 0, b = 0.7, unit = 'cm')
  ) +
  xlab('') + ylab('')


p_density_top<- all_results_fixed_effects %>%
  filter(Complex %in% c('2b18', '3ulh')) %>%
  mutate(Complex = case_when(
    Complex == '2b18' ~ '2B18 (HET dominant)', 
    Complex == '3ulh' ~ '3ULH (HM dominant)'
  )) %>%
  ggplot(aes(x = mut_eff_HM, fill = Complex)) +
  geom_density(alpha = 0.5 # , scale = 'width'
  ) +
  scale_fill_manual(values = c('#4575B4', '#FC8D59')) +
  xlim(-3, 5) +
  theme(legend.position = 'none', 
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        # axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(),
        # axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        plot.margin = margin(t = 0, b = -0.5, l = 0.5, unit = 'cm')
  ) +
  xlab('') + ylab('')

p_fig3B <- plot_grid(p_density_top, NULL, 
                     p_center, p_density_right, ncol = 2, nrow = 2, 
                     rel_widths = c(1, 0.15), rel_heights = c(0.15, 1))


p_fig3B

## Put panels A and B together with the appropriate legend
p_fig3AB_nolegend <- plot_grid(p_fig3A, p_fig3B, ncol = 2, labels = c('A', 'B'), 
                               label_size = panel_label_size, label_fontface = 'bold')

p_fig3AB_legend <- plot_grid(p_legend, p_fig3AB_nolegend, nrow = 2, rel_heights = c(0.085, 1))

p_fig3AB_legend

## Figure 3C: Cumulative distribution of residuals for available mutations ##
mut_effects_new <- left_join(x = mut_effects %>%
                               mutate(Complex = toupper(Complex)), 
                             y = all_data_summary_final_points %>% ungroup() %>%
                               mutate(mean_HET_HM = mean_pct_HET + mean_pct_HMs) %>%
                               mutate(Outcome = case_when(
                                 mean_pct_HET >= 70 ~ 'HET dominant', 
                                 mean_pct_HMs >= 70 ~ 'HM dominant', 
                                 mean_HET_HM >= 70 ~ 'Both HM and HET',
                                 mean_pct_mono >= 70 ~ 'Monomer',
                                 TRUE ~ 'Ambiguous')
                               ) %>%
                               select(Complex, Outcome),
                             by = c('Complex' = 'Complex'))

## Add another boolean variable to distinguish alpha
mut_effects_new %<>% mutate(alpha_check = case_when(
  Complex %in% c('2B18', '3ULH') ~ TRUE,
  TRUE ~ FALSE
))

mut_effects_new %<>% mutate(linear_residuals = Mean_ddG_int_HET - (0.5 * Mean_ddG_int_HM)) %>% 
  filter(Mean_ddG_int_HET != 0, Mean_ddG_int_HM != 0) %>%
  mutate(Outcome = factor(Outcome, levels = c('HET dominant', 'Both HM and HET', 'HM dominant')))

p_fig3C <- mut_effects_new %>% ggplot(aes(x = linear_residuals, group = Complex, colour = Outcome,
                                          alpha = alpha_check, linewidth = alpha_check)) +
  geom_vline(xintercept = 0.2, linetype = 'dashed') +
  geom_vline(xintercept = -0.2, linetype = 'dashed') +
  annotate(geom = 'rect', xmin = 0.2, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = 'grey', alpha = 0.3, colour = NA) +
  annotate(geom = 'rect', xmin = -Inf, xmax = -0.2, ymin = -Inf, ymax = Inf,
           fill = 'grey', alpha = 0.3) +
  stat_ecdf() +
  stat_ecdf(data = mut_effects_new %>% filter(alpha_check == T)) +
  scale_linewidth_manual(values = c(0.5, 0.75)) +
  scale_colour_manual(values = c('#4575B4', '#008000', '#FC8D59')) +
  scale_alpha_manual(values = c(0.05, 1)) +
  xlim(-0.5, 0.5) +
  annotate(geom = 'text', x = -0.375, y = 0.5, label =  'Favoring HET\n(negative residual)', size = 3) +
  annotate(geom = 'text', x = 0.375, y = 0.5, label = 'Favoring HM\n(positive residual)', size = 3) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  theme(legend.position = 'none',
        legend.justification = 'center', 
        panel.grid.major = element_blank(), 
        plot.title = element_text(hjust = 0.5, size = 14), 
        axis.title = element_text(size = 9)) +
  xlab('Residuals') + ylab('Cumulative density') +
  guides(alpha = 'none', linewidth = 'none', fill = 'none')
p_fig3C

## Figure 3D: Cumulative distribution of residuals for fixed mutations ##

## Add the outcomes for each structure
all_results_fixed_effects_outcome <- left_join(x = all_results_fixed_effects %>% ungroup() %>%
                                                 filter(Start_stab == -5, Start_aff == -10) %>%
                                                 mutate(netEffect_HETbind = mut_eff_binding_energy_AB, 
                                                        # Mutations can only affect one of the two HMs at a time, so the registered change
                                                        # has to be zero for one of them
                                                        netEffect_HMbind = mut_eff_binding_energy_AA + mut_eff_binding_energy_BB, 
                                                        Complex = toupper(Complex)
                                                 ) %>%
                                                 select(Replicate, Complex,
                                                        netEffect_HETbind, netEffect_HMbind) %>% 
                                                 # Calculate residuals
                                                 mutate(linear_residuals = netEffect_HETbind - (0.5 * netEffect_HMbind)), 
                                               y = all_data_summary_final_points %>% ungroup() %>%
                                                 mutate(mean_HET_HM = mean_pct_HET + mean_pct_HMs) %>%
                                                 mutate(Outcome = case_when(
                                                   mean_pct_HET >= 70 ~ 'HET dominant', 
                                                   mean_pct_HMs >= 70 ~ 'HM dominant', 
                                                   mean_HET_HM >= 70 ~ 'Both HM and HET',
                                                   mean_pct_mono >= 70 ~ 'Monomer',
                                                   TRUE ~ 'Ambiguous')
                                                 ) %>%
                                                 select(Complex, Outcome),
                                               by = c('Complex' = 'Complex'))

## Add another boolean variable to distinguish alpha
all_results_fixed_effects_outcome %<>% mutate(alpha_check = case_when(
  Complex %in% c('2b18', '3ulh') ~ TRUE,
  TRUE ~ FALSE
))

p_fig3D <- all_results_fixed_effects_outcome %>%
  filter(netEffect_HETbind != 0, netEffect_HMbind != 0) %>%
  mutate(Outcome = factor(Outcome, levels = c('HET dominant', 'Both HM and HET', 'HM dominant'))) %>%
  ggplot(aes(x = linear_residuals, group = Complex, colour = Outcome,)) +
  geom_vline(xintercept = 0.2, linetype = 'dashed') +
  annotate(geom = 'rect', xmin = 0.2, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = 'grey', alpha = 0.3, colour = NA) +
  geom_vline(xintercept = -0.2, linetype = 'dashed') +
  annotate(geom = 'rect', xmin = -Inf, xmax = -0.2, ymin = -Inf, ymax = Inf,
           fill = 'grey', alpha = 0.3) +
  stat_ecdf(linewidth = 0.5, alpha = 0.05) +
  stat_ecdf(data = all_results_fixed_effects_outcome %>%
              filter(netEffect_HETbind != 0, netEffect_HMbind != 0, 
                     Complex %in% c('2B18', '3ULH')
              ) %>%
              mutate(Outcome = factor(Outcome, levels = c('HET dominant', 'HM dominant'))), 
            alpha = 1, linewidth = 0.75
  ) +
  scale_colour_manual(values = c('#4575B4', '#008000', '#FC8D59')) +
  xlim(-0.5, 0.5) +
  ylim(0, 1) +
  annotate(geom = 'text', x = -0.375, y = 0.5, label =  'Favoring HET\n(negative residual)', size = 3) +
  annotate(geom = 'text', x = 0.375, y = 0.5, label = 'Favoring HM\n(positive residual)', size = 3) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  theme(legend.position = 'none',
        legend.justification = 'center', 
        panel.grid.major = element_blank(), 
        plot.title = element_text(hjust = 0.5, size = 14), 
        axis.title = element_text(size = 9)) +
  xlab('Residuals') + ylab('Cumulative density') +
  guides(alpha = 'none', linewidth = 'none', fill = 'none')
p_fig3D

p_legend_fig3CD <- get_legend(p_fig3D + theme(legend.position = 'top'))

p_fig3CD_nolegend <- plot_grid(p_fig3C, p_fig3D, ncol = 2, labels = c('C', 'D'), 
                               label_size = panel_label_size, label_fontface = 'bold')

p_fig3CD <- plot_grid(p_legend_fig3CD, p_fig3CD_nolegend, nrow = 2, rel_heights = c(0.085, 1))
p_fig3CD

## Figure 3E: Correlation between enrichment of available HET favoring residuals and HET percentage ##
fractions_mut_effects <- mut_effects_new %>% ungroup() %>%
  group_by(Complex, Outcome) %>%
  summarise(fraction_favor_HET = mean(linear_residuals < -0.2), 
            fraction_favor_HM = mean(linear_residuals > 0.2)) %>%
  mutate(delta_favor_HET = fraction_favor_HET - fraction_favor_HM)

fractions_mut_effects_final_HET <- left_join(x = fractions_mut_effects %>%
                                               mutate(Complex = toupper(Complex)), 
                                             y = all_data_summary_final_points %>% 
                                               select(Complex, mean_pct_HET), 
                                             by = c('Complex' = 'Complex'))

p_fig3E <- fractions_mut_effects_final_HET %>%
  ggplot(aes(x = delta_favor_HET, y = mean_pct_HET)) +
  geom_point() +
  xlab('Enrichment of available mutations with HET-favoring residuals') + ylab('Mean HET percentage (%)') +
  theme(panel.grid.major = element_blank(), 
        plot.title = element_text(hjust = 0.5), 
        axis.title = element_text(size = 9)) +
  stat_cor(method = 'spearman', p.accuracy = 0.001, r.accuracy = 0.01, cor.coef.name = 'rho',
           size = 3, label.x.npc = 0.5, label.y.npc = 0.1) +
  geom_smooth(method = 'lm', alpha = 0.5) +
  ylim(0, 100)
p_fig3E

## Figure 3F: Correlation between enrichment of fixed HET favoring residuals and HET percentage ##
fractions_fixed_mut_effects <- 
  all_results_fixed_effects_outcome %>%
  ungroup() %>%
  group_by(Complex, Outcome) %>%
  summarise(fraction_favor_HET = mean(linear_residuals < -0.2), 
            fraction_favor_HM = mean(linear_residuals > 0.2)) %>%
  mutate(delta_favor_HET = fraction_favor_HET - fraction_favor_HM)

fractions_fixed_mut_effects_final_HET <- left_join(x = fractions_fixed_mut_effects %>%
                                                     mutate(Complex = toupper(Complex)), 
                                                   y = all_data_summary_final_points %>% 
                                                     select(Complex, mean_pct_HET), 
                                                   by = c('Complex' = 'Complex'))

p_fig3F <- fractions_fixed_mut_effects_final_HET %>%
  ggplot(aes(x = delta_favor_HET, y = mean_pct_HET)) +
  geom_point() +
  xlab('Enrichment of fixed mutations with HET-favoring residuals') + ylab('Mean HET percentage (%)') +
  theme(panel.grid.major = element_blank(), 
        plot.title = element_text(hjust = 0.5), 
        axis.title = element_text(size = 9)) +
  stat_cor(method = 'spearman', p.accuracy = 0.001, r.accuracy = 0.01, cor.coef.name = 'rho',
           size = 3, label.x.npc = 0.5, label.y.npc = 0.1) +
  geom_smooth(method = 'lm', alpha = 0.5) +
  ylim(0, 100)
p_fig3F

## Put all the panels together
p_fig3EF <- plot_grid(p_fig3E, p_fig3F, ncol = 2, labels = c('E', 'F'), 
                      label_size = panel_label_size, label_fontface = 'bold')

p_fig3 <-plot_grid(p_fig3AB_legend, p_fig3CD, p_fig3EF, nrow = 3)

ggsave(p_fig3, device = cairo_pdf, width = 24, height = 24, units = 'cm', dpi = 300, 
       filename = 'Figures/Main_figures/3.Fig3_available_fixed_effects.pdf')
## Labels of fixed and available mutations were added manually

#### Figure EV3: Like figure 3 but only with the very high quality structures ####

axis_title_size_fig3 = 10

# Load QSbio data
qsbio <- read_delim(file = 'Data/QSbio/QSbio.csv',
                    delim = '\t', locale = locale(decimal_mark = ','))

qsbio_hm <- qsbio %>%
  filter(pdb_nsub == 2) %>% 
  mutate(code_sub = toupper(code_sub))

all_data_summary_final_points_outcome_qsbio <- left_join(x = all_data_summary_final_points_outcome,
                                                         y = qsbio_hm %>% 
                                                           separate(code, into = c('code_sub_new', 'bio_assembly_num'), sep = '_') %>%
                                                           filter(bio_assembly_num == 1) %>%
                                                           select(code_sub, error_estimated, QSbio.confidence), 
                                                         by = c("Complex" = "code_sub"))
all_data_summary_final_points_long_qsbio <- all_data_summary_final_points_outcome_qsbio %>% ungroup() %>%
  pivot_longer(cols = c(mean_pct_HET, mean_pct_HMs, mean_pct_mono), 
               names_to = 'assembly', values_to = 'mean_value')

## Fig. EV3A: Like fig. 2C but only with the very high quality structures
p_figEV3A <- all_data_summary_final_points_long_qsbio %>% ungroup() %>%
  filter(QSbio.confidence == 'Very high') %>%
  mutate(assembly_label = ifelse(assembly == 'mean_pct_HET', 'HET', 
                                 ifelse(assembly == 'mean_pct_HMs', 'HMs (AA+BB)',
                                        'Monomers (A+B)')), 
         Complex = toupper(Complex)) %>%
  mutate(assembly_label = factor(assembly_label, levels = c('HET', 'Monomers (A+B)', 'HMs (AA+BB)')),
         Complex = factor(Complex, levels = toupper(all_data_summary_final_points$Complex))
  ) %>%
  ggplot(aes(x = Complex, y = mean_value, fill = assembly_label)) +
  geom_bar(stat = 'identity') +
  xlab('PDB ID') + ylab('Percentage (%)') +
  labs(fill = '') +
  geom_hline(yintercept = 50, linetype = 'dashed') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 9),
        axis.title = element_text(size = axis_title_size_fig3),
        axis.text.y = element_text(size = 9),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.position = 'top', 
        legend.justification = 'center') +
  guides(fill = guide_legend(override.aes = list(size = 0.25))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  scale_fill_manual(values = c('#4575B4', 'grey', '#FC8D59' ))
p_figEV3A

## Fig. EV3B ##
mut_effects_new <- left_join(x = mut_effects %>%
                               mutate(Complex = toupper(Complex)), 
                             ## When filtering for only QSBio very high confidence structures
                             y = all_data_summary_final_points_outcome_qsbio %>% ungroup() %>%
                               filter(QSbio.confidence == 'Very high') %>%
                               mutate(mean_HET_HM = mean_pct_HET + mean_pct_HMs) %>%
                               mutate(Outcome = case_when(
                                 mean_pct_HET >= 70 ~ 'HET dominant', 
                                 mean_pct_HMs >= 70 ~ 'HM dominant', 
                                 mean_HET_HM >= 70 ~ 'Both HM and HET',
                                 mean_pct_mono >= 70 ~ 'Monomer',
                                 TRUE ~ 'Ambiguous')
                               ) %>%
                               select(Complex, Outcome, QSbio.confidence),
                             by = c('Complex' = 'Complex'))

## Add another boolean variable to distinguish alpha
mut_effects_new %<>% mutate(alpha_check = case_when(
  Complex %in% c('2B18', '3ULH') ~ TRUE,
  TRUE ~ FALSE
))

mut_effects_new %<>% mutate(linear_residuals = Mean_ddG_int_HET - (0.5 * Mean_ddG_int_HM)) %>% 
  filter(Mean_ddG_int_HET != 0, Mean_ddG_int_HM != 0) %>%
  mutate(Outcome = factor(Outcome, levels = c('HET dominant', 'Both HM and HET', 'HM dominant'))) %>%
  filter(QSbio.confidence == 'Very high')

p_figEV3B <- mut_effects_new %>% ggplot(aes(x = linear_residuals, group = Complex, colour = Outcome,
                                           alpha = alpha_check, linewidth = alpha_check)) +
  geom_vline(xintercept = 0.2, linetype = 'dashed') +
  geom_vline(xintercept = -0.2, linetype = 'dashed') +
  annotate(geom = 'rect', xmin = 0.2, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = 'grey', alpha = 0.3, colour = NA) +
  annotate(geom = 'rect', xmin = -Inf, xmax = -0.2, ymin = -Inf, ymax = Inf,
           fill = 'grey', alpha = 0.3) +
  stat_ecdf() +
  stat_ecdf(data = mut_effects_new %>% filter(alpha_check == T)) +
  scale_linewidth_manual(values = c(0.5, 0.75)) +
  scale_colour_manual(values = c('#4575B4', '#008000', '#FC8D59')) +
  scale_alpha_manual(values = c(0.2, 1)) +
  xlim(-0.5, 0.5) +
  annotate(geom = 'text', x = -0.375, y = 0.5, label =  'Favoring HET\n(negative residual)', size = 3.5) +
  annotate(geom = 'text', x = 0.375, y = 0.5, label = 'Favoring HM\n(positive residual)', size = 3.5) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  theme(legend.position = 'none',
        legend.justification = 'center', 
        panel.grid.major = element_blank(), 
        plot.title = element_text(hjust = 0.5, size = 14), 
        axis.title = element_text(size = axis_title_size_fig3)) +
  xlab('Residuals') + ylab('Cumulative density') +
  guides(alpha = 'none', linewidth = 'none', fill = 'none')
p_figEV3B

## Fig. EV3C ##
## Add the outcomes for each structure
all_results_fixed_effects_outcome <- left_join(x = all_results_fixed_effects %>% ungroup() %>%
                                                 filter(Start_stab == -5, Start_aff == -10) %>%
                                                 mutate(netEffect_HETbind = mut_eff_binding_energy_AB, 
                                                        # Mutations can only affect one of the two HMs at a time, so the registered change
                                                        # has to be zero for one of them
                                                        netEffect_HMbind = mut_eff_binding_energy_AA + mut_eff_binding_energy_BB, 
                                                        Complex = toupper(Complex)
                                                 ) %>%
                                                 select(Replicate, Complex,
                                                        netEffect_HETbind, netEffect_HMbind) %>% 
                                                 # Calculate residuals
                                                 mutate(linear_residuals = netEffect_HETbind - (0.5 * netEffect_HMbind)), 
                                               y = all_data_summary_final_points_outcome_qsbio %>% ungroup() %>%
                                                 filter(QSbio.confidence == 'Very high') %>%
                                                 mutate(mean_HET_HM = mean_pct_HET + mean_pct_HMs) %>%
                                                 mutate(Outcome = case_when(
                                                   mean_pct_HET >= 70 ~ 'HET dominant', 
                                                   mean_pct_HMs >= 70 ~ 'HM dominant', 
                                                   mean_HET_HM >= 70 ~ 'Both HM and HET',
                                                   mean_pct_mono >= 70 ~ 'Monomer',
                                                   TRUE ~ 'Ambiguous')
                                                 ) %>%
                                                 select(Complex, Outcome, QSbio.confidence),
                                               by = c('Complex' = 'Complex'))

## Add another boolean variable to distinguish alpha
all_results_fixed_effects_outcome %<>% mutate(alpha_check = case_when(
  Complex %in% c('2b18', '3ulh') ~ TRUE,
  TRUE ~ FALSE
)) %>%
  filter(QSbio.confidence == 'Very high')

p_figEV3C <- all_results_fixed_effects_outcome %>%
  filter(netEffect_HETbind != 0, netEffect_HMbind != 0) %>%
  mutate(Outcome = factor(Outcome, levels = c('HET dominant', 'Both HM and HET', 'HM dominant'))) %>%
  ggplot(aes(x = linear_residuals, group = Complex, colour = Outcome,)) +
  geom_vline(xintercept = 0.2, linetype = 'dashed') +
  annotate(geom = 'rect', xmin = 0.2, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = 'grey', alpha = 0.3, colour = NA) +
  geom_vline(xintercept = -0.2, linetype = 'dashed') +
  annotate(geom = 'rect', xmin = -Inf, xmax = -0.2, ymin = -Inf, ymax = Inf,
           fill = 'grey', alpha = 0.3) +
  stat_ecdf(linewidth = 0.5, alpha = 0.2) +
  stat_ecdf(data = all_results_fixed_effects_outcome %>%
              filter(netEffect_HETbind != 0, netEffect_HMbind != 0, 
                     Complex %in% c('2B18', '3ULH')
              ) %>%
              mutate(Outcome = factor(Outcome, levels = c('HET dominant', 'HM dominant'))), 
            alpha = 1, linewidth = 0.75
  ) +
  scale_colour_manual(values = c('#4575B4', '#008000', '#FC8D59')) +
  xlim(-0.5, 0.5) +
  ylim(0, 1) +
  annotate(geom = 'text', x = -0.375, y = 0.5, label =  'Favoring HET\n(negative residual)', size = 3.5) +
  annotate(geom = 'text', x = 0.375, y = 0.5, label = 'Favoring HM\n(positive residual)', size = 3.5) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  theme(legend.position = 'none',
        legend.justification = 'center', 
        panel.grid.major = element_blank(), 
        plot.title = element_text(hjust = 0.5, size = 14), 
        axis.title = element_text(size = axis_title_size_fig3)) +
  xlab('Residuals') + ylab('Cumulative density') +
  guides(alpha = 'none', linewidth = 'none', fill = 'none')
p_figEV3C

p_legend_figEV3BC <- get_legend(p_figEV3C + theme(legend.position = 'top'))

p_figEV3BC_nolegend <- plot_grid(p_figEV3B, p_figEV3C, ncol = 2, labels = c('B', 'C'), 
                                label_size = panel_label_size, label_fontface = 'bold')

p_figEV3BC <- plot_grid(p_legend_figEV3BC, p_figEV3BC_nolegend, nrow = 2, rel_heights = c(0.085, 1))
p_figEV3BC

## Fig. EV3D ##
fractions_mut_effects <- mut_effects_new %>% ungroup() %>%
  group_by(Complex, Outcome) %>%
  summarise(fraction_favor_HET = mean(linear_residuals < -0.2), 
            fraction_favor_HM = mean(linear_residuals > 0.2)) %>%
  mutate(delta_favor_HET = fraction_favor_HET - fraction_favor_HM)

fractions_mut_effects_final_HET <- left_join(x = fractions_mut_effects %>%
                                               mutate(Complex = toupper(Complex)), 
                                             y = all_data_summary_final_points %>% 
                                               select(Complex, mean_pct_HET), 
                                             by = c('Complex' = 'Complex'))

p_figEV3D <- fractions_mut_effects_final_HET %>%
  ggplot(aes(x = delta_favor_HET, y = mean_pct_HET)) +
  geom_point() +
  xlab('Enrichment of available mutations with HET-favoring residuals') + ylab('Mean HET percentage (%)') +
  theme(panel.grid.major = element_blank(), 
        plot.title = element_text(hjust = 0.5), 
        axis.title = element_text(size = axis_title_size_fig3)) +
  stat_cor(method = 'spearman', p.accuracy = 0.001, r.accuracy = 0.01, cor.coef.name = 'rho',
           size = 3.5, label.x.npc = 0.5, label.y.npc = 0.1) +
  geom_smooth(method = 'lm', alpha = 0.5) +
  ylim(0, 100)
p_figEV3D

## Fig. EV3E ##
fractions_fixed_mut_effects <- 
  all_results_fixed_effects_outcome %>%
  ungroup() %>%
  group_by(Complex, Outcome) %>%
  summarise(fraction_favor_HET = mean(linear_residuals < -0.2), 
            fraction_favor_HM = mean(linear_residuals > 0.2)) %>%
  mutate(delta_favor_HET = fraction_favor_HET - fraction_favor_HM)

fractions_fixed_mut_effects_final_HET <- left_join(x = fractions_fixed_mut_effects %>%
                                                     mutate(Complex = toupper(Complex)), 
                                                   y = all_data_summary_final_points %>% 
                                                     select(Complex, mean_pct_HET), 
                                                   by = c('Complex' = 'Complex'))

p_figEV3E <- fractions_fixed_mut_effects_final_HET %>%
  ggplot(aes(x = delta_favor_HET, y = mean_pct_HET)) +
  geom_point() +
  xlab('Enrichment of fixed mutations with HET-favoring residuals') + ylab('Mean HET percentage (%)') +
  theme(panel.grid.major = element_blank(), 
        plot.title = element_text(hjust = 0.5), 
        axis.title = element_text(size = axis_title_size_fig3)) +
  stat_cor(method = 'spearman', p.accuracy = 0.001, r.accuracy = 0.01, cor.coef.name = 'rho',
           size = 3.5, label.x.npc = 0.5, label.y.npc = 0.1) +
  geom_smooth(method = 'lm', alpha = 0.5) +
  ylim(0, 100)
p_figEV3E

## Put all the panels together
p_figEV3DE <- plot_grid(p_figEV3D, p_figEV3E, ncol = 2, labels = c('D', 'E'), 
                       label_size = panel_label_size, label_fontface = 'bold')

p_figEV3 <-plot_grid(p_figEV3A, p_figEV3BC, p_figEV3DE, nrow = 3, labels = c('A', '', ''), 
                    label_size = panel_label_size, label_fontface = 'bold')

ggsave(p_figEV3, device = cairo_pdf, width = 26, height = 24, units = 'cm', dpi = 300, 
       filename = 'Figures/ExpandedView_figures/FigEV3.Available_fixed_effects_QSBio_veryHigh.pdf')

#### Figure S07: Relationship between cumulative residuals and concentration of HET at the time ####

values_check <- c(5, 10, 25, 50, 100, 200)
all_plots <- list()
all_het_pct <- c()

for(i in 1:length(values_check)){
  n_mut_check <- values_check[i]
  
  # Get the net effect on HM and HET binding
  all_results_fixed_effects_tmp <- all_results_fixed_effects %>% ungroup() %>%
    filter(Start_stab == -5, Start_aff == -10) %>%
    mutate(netEffect_HETbind = mut_eff_binding_energy_AB, 
           # Mutations can only affect one of the two HMs at a time, so the registered change
           # has to be zero for one of them
           netEffect_HMbind = mut_eff_binding_energy_AA + mut_eff_binding_energy_BB, 
           total_complexes = cA + cB + cAA + cAB + cBB) %>%
    mutate(pct_HET = (cAB * 100) / total_complexes) %>%
    select(Replicate, Complex, fixed_mut,
           netEffect_HETbind, netEffect_HMbind, pct_HET) %>% 
    group_by(Replicate, Complex) %>%
    # Calculate residuals
    mutate(linear_residuals = netEffect_HETbind - (0.5 * netEffect_HMbind), 
           mut_num = row_number())
  
  
  # Get the percentage of HET after the mutations that accumulated
  current_pct_HET <- all_results_fixed_effects_tmp %>%
    filter(fixed_mut == n_mut_check)
  
  # Get the first n mutations
  all_results_fixed_effects_summary <- all_results_fixed_effects_tmp %>%
    group_by(Replicate, Complex) %>%
    filter(fixed_mut <= n_mut_check) %>%
    summarise(fixed_residual_sum = sum(linear_residuals))
  
  
  all_results_fixed_effects_summary <- inner_join(x = all_results_fixed_effects_summary,
                                                  y = current_pct_HET %>%
                                                    select(Replicate, Complex, pct_HET),
                                                  by = c('Replicate' = 'Replicate',
                                                         'Complex' = 'Complex'))
  
  ## Show in a figure the relation between the first fixed residuals and the outcome
  p <- all_results_fixed_effects_summary %>%
    ggplot(aes(x = fixed_residual_sum, y = pct_HET)) +
    geom_point(alpha = 0.2) +
    xlim(-7, 7) +
    ggtitle(str_c('First ', toString(n_mut_check), ' mutations', sep = '')) +
    theme(plot.title = element_text(size = 14, hjust = 0.5), 
          panel.grid.major = element_blank()) +
    geom_hline(yintercept = 50, linetype = 'dashed') +
    xlab('Cumulative residual sum') + ylab('HET percentage (%)')
  
  ## Save the figure to the list
  all_plots[[i]] <- p
  
  ## Accumulate all the values for pct_HET to show the distributions
  all_het_pct <- bind_rows(all_het_pct, 
                           current_pct_HET %>% select(Replicate, Complex, pct_HET, mut_num))
  
}

p_cumul_res <- plot_grid(all_plots[[1]], all_plots[[2]], all_plots[[3]],
                         all_plots[[4]], all_plots[[5]], all_plots[[6]], 
                         nrow = 2)
ggsave(p_cumul_res, device = cairo_pdf, width = 20, height = 14, units = 'cm', dpi = 300, 
       filename = 'Figures/Supplementary_figures/Appendix_FigS07_residuals_HET_at_the_time.pdf')

#### Figure 4 ####

#### Figure 4A ####

#### Repeat but using log2 of the synthesis ratio on the x-axis ####

data_syn_ratio_bind <- read_delim('Results_solution_space/solutionSpace_syn_ratio_diff_binding_log2.tsv',
                                  delim = '\t')

# Organize as a matrix for a heatmap
data_syn_ratio_bind_wide <- data_syn_ratio_bind %>% select(log2_val, binding_diff, pct_AB) %>%
  arrange(log2_val, desc(binding_diff)) %>%
  pivot_wider(names_from = log2_val, values_from = pct_AB)

# Move the first column to the index
first_column <- data_syn_ratio_bind_wide$binding_diff
data_syn_ratio_bind_wide <- as.matrix(data_syn_ratio_bind_wide[, 2:ncol(data_syn_ratio_bind_wide)])
rownames(data_syn_ratio_bind_wide) <- first_column

## Draw heatmaps for effect of synthesis rate

# Define column labels
column_labels = seq(from = -2, to = 1.95, by = 0.05)
new_column_labels = column_labels
for(i in 1:length(column_labels)){
  
  cur_val = column_labels[i]
  new_val = ifelse(mod(round(cur_val*10, 2), 5) == 0, toString(cur_val), '')
  
  new_column_labels[i] = new_val
  
}

# Define row labels
row_labels = rev(seq(from = -3, to = 2.9, by = 0.1))
new_row_labels = row_labels
for(i in 1:length(row_labels)){
  
  cur_val = row_labels[i]
  new_val = ifelse(mod(cur_val*10, 10) == 0, toString(cur_val), '')
  
  new_row_labels[i] = new_val
  
}

## New variables for heatmaps (adjusting for new panel size) ##
row_title_size = 9
column_title_size = 9
row_name_size = 8
column_name_size = 8
legend_title_size = 9
legend_label_size = 8
heatmap_width = 5
heatmap_height = 5
legend_height = 4
legend_width = 0.3
ht_opt$TITLE_PADDING = unit(c(0.5, 0.5), "points")

p_fig4A <- Heatmap(
  data_syn_ratio_bind_wide,
  cluster_columns = F, cluster_rows = F,
  col = colorRamp2(
    breaks = c(0, 20, 40, 60, 80, 100),
    colors = viridis(6)),
  show_column_names = T, column_names_side = 'bottom',
  column_names_rot = 90, column_names_centered = T,
  column_names_gp = gpar(fontsize=column_name_size, fontface='bold', align = 'bottom', colour = 'black'),
  column_labels = new_column_labels,

  show_row_names = T, row_names_side = 'left',
  row_names_rot = 90,  row_names_centered = T,
  row_names_gp = gpar(fontsize=row_name_size, fontface = 'bold', align = 'left', colour = 'black'),
  row_labels = new_row_labels,

  width=unit(heatmap_width, 'cm'), height = unit(heatmap_height, 'cm'),
  border = T,
  row_title = expression(paste(bold('\u0394'),
                               bold(G['bind,HET AB']), 
                               bold(' - \u0394'),
                               bold(G['bind,HM AA']), 
                               bold(' (kcal/mol)'), sep = '')),
  column_title = expression(paste(bold('log2('), 
                                  bold(s['B']),
                                  bold(' / '),
                                  bold(s['A']), 
                                  bold(')'))),
  row_title_gp = gpar(fontsize=row_title_size, fontface = 'bold'),
  column_title_gp = gpar(fontsize=column_title_size, fontface = 'bold'),
  column_title_side = 'bottom',
  
  
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(
    at = c(0, 20, 40, 60, 80, 100),
    title = "Percentage of HET AB (%)",
    title_gp = gpar(fontsize = legend_title_size),
    legend_height = unit(legend_height, "cm"),
    grid_width = unit(legend_width, 'cm'),
    border='black',
    lwd=1.7,
    labels_gp = gpar(fontsize = legend_label_size),
    title_position = "leftcenter-rot"
  ),
  rect_gp = gpar(col = 'transparent',
                 lwd = 0
  ), use_raster = T
)
p_fig4A

#### Figure 4B ####
## Repeat for the HM BB

# Organize as a matrix for a heatmap
data_syn_ratio_bind_wide <- data_syn_ratio_bind %>%
  mutate(pct_BB = 100 * cBB / (cA + cB + cAA + cAB + cBB)) %>%
  select(log2_val, binding_diff, pct_BB) %>%
  arrange(log2_val, desc(binding_diff)) %>%
  pivot_wider(names_from = log2_val, values_from = pct_BB)

# Move the first column to the index
first_column <- data_syn_ratio_bind_wide$binding_diff
data_syn_ratio_bind_wide <- as.matrix(data_syn_ratio_bind_wide[, 2:ncol(data_syn_ratio_bind_wide)])
rownames(data_syn_ratio_bind_wide) <- first_column

## Draw heatmaps for effect of synthesis rate

# Define column labels
column_labels = seq(from = -2, to = 1.95, by = 0.05)
new_column_labels = column_labels
for(i in 1:length(column_labels)){
  
  cur_val = column_labels[i]
  new_val = ifelse(mod(round(cur_val*10, 2), 5) == 0, toString(cur_val), '')
  
  new_column_labels[i] = new_val
  
}

# Define row labels
row_labels = rev(seq(from = -3, to = 2.9, by = 0.1))
new_row_labels = row_labels
for(i in 1:length(row_labels)){
  
  cur_val = row_labels[i]
  new_val = ifelse(mod(cur_val*10, 10) == 0, toString(cur_val), '')
  
  new_row_labels[i] = new_val
  
}

p_fig4B <- Heatmap(
  data_syn_ratio_bind_wide,
  cluster_columns = F, cluster_rows = F,
  col = colorRamp2(
    breaks = c(0, 20, 40, 60, 80, 100),
    colors = viridis(6)),
  show_column_names = T, column_names_side = 'bottom',
  column_names_rot = 90, column_names_centered = T,
  column_names_gp = gpar(fontsize=column_name_size, fontface='bold', align = 'bottom', colour = 'black'),
  column_labels = new_column_labels,

  show_row_names = T, row_names_side = 'left',
  row_names_rot = 90,  row_names_centered = T,
  row_names_gp = gpar(fontsize=row_name_size, fontface = 'bold', align = 'left', colour = 'black'),
  row_labels = new_row_labels,

  width=unit(heatmap_width, 'cm'), height = unit(heatmap_height, 'cm'),
  border = T,
  row_title = expression(paste(bold('\u0394'),
                               bold(G['bind,HET AB']), 
                               bold(' - \u0394'),
                               bold(G['bind,HM AA']), 
                               bold(' (kcal/mol)'), sep = '')),
  column_title = expression(paste(bold('log2('),
                                  bold(s['B']),
                                  bold(' / '),
                                  bold(s['A']), 
                                  bold(')'))),
  row_title_gp = gpar(fontsize=row_title_size, fontface = 'bold'),
  column_title_gp = gpar(fontsize=column_title_size, fontface = 'bold'),
  column_title_side = 'bottom',
  
  
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(
    at = c(0, 20, 40, 60, 80, 100),
    title = "Percentage of HM BB (%)",
    title_gp = gpar(fontsize = legend_title_size),
    legend_height = unit(legend_height, "cm"),
    grid_width = unit(legend_width, "cm"),
    border='black',
    lwd=1.7,
    labels_gp = gpar(fontsize = legend_label_size),
    title_position = "leftcenter-rot"
  ),
  rect_gp = gpar(col = 'transparent',
                 lwd = 0
  ), use_raster = T
)
p_fig4B

## Figure 4C: A heatmap of the total activity in the solution space ##

## Look at the total activity
# Organize as a matrix for a heatmap
data_syn_ratio_bind_wide <- data_syn_ratio_bind %>%
  mutate(total_activity = 0.1*cA + 0.1*cB + 1*cAA + 1*cBB + 1*cAB) %>% 
  select(log2_val, binding_diff, total_activity) %>%
  arrange(log2_val, desc(binding_diff)) %>%
  pivot_wider(names_from = log2_val, values_from = total_activity)

# Move the first column to the index
first_column <- data_syn_ratio_bind_wide$binding_diff
data_syn_ratio_bind_wide <- as.matrix(data_syn_ratio_bind_wide[, 2:ncol(data_syn_ratio_bind_wide)])
rownames(data_syn_ratio_bind_wide) <- first_column

## Draw heatmaps for effect of synthesis rate

# Define column labels
column_labels = seq(from = -2, to = 1.95, by = 0.05)
new_column_labels = column_labels
for(i in 1:length(column_labels)){
  
  cur_val = column_labels[i]
  new_val = ifelse(mod(round(cur_val*10, 2), 5) == 0, toString(cur_val), '')
  
  new_column_labels[i] = new_val
  
}

# Define row labels
row_labels = rev(seq(from = -3, to = 2.9, by = 0.1))
new_row_labels = row_labels
for(i in 1:length(row_labels)){
  
  cur_val = row_labels[i]
  new_val = ifelse(mod(cur_val*10, 10) == 0, toString(cur_val), '')
  
  new_row_labels[i] = new_val
  
}

p_fig4C_activity <- Heatmap(
  data_syn_ratio_bind_wide,
  cluster_columns = F, cluster_rows = F,
  col = colorRamp2(
    breaks = c(0, 40, 80, 120, 160),
    colors = magma(5)
  ),
  show_column_names = T, column_names_side = 'bottom',
  column_names_rot = 90, column_names_centered = T,
  column_names_gp = gpar(fontsize=column_name_size, fontface='bold', align = 'bottom', colour = 'black'),
  column_labels = new_column_labels,

  show_row_names = T, row_names_side = 'left',
  row_names_rot = 90,  row_names_centered = T,
  row_names_gp = gpar(fontsize=row_name_size, fontface = 'bold', align = 'left', colour = 'black'),
  row_labels = new_row_labels,

  width=unit(heatmap_width, 'cm'), height = unit(heatmap_height, 'cm'),
  border = T,
  row_title = expression(paste(
    bold('\u0394'),
    bold(G['bind,HET AB']), 
    bold(' - \u0394'),
    bold(G['bind,HM AA']), 
    bold(' (kcal/mol)'),
    sep = '')),
  column_title = expression(paste(bold('log2('),
                                  bold(s['B']),
                                  bold(' / '),
                                  bold(s['A']), 
                                  bold(')'))),
  row_title_gp = gpar(fontsize=row_title_size, fontface = 'bold'),
  column_title_gp = gpar(fontsize=column_title_size, fontface = 'bold'),
  column_title_side = 'bottom',
  
  
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(
    title = "Total activity",
    title_gp = gpar(fontsize = legend_title_size),
    legend_height = unit(legend_height, "cm"),
    grid_width = unit(legend_width, "cm"),
    border='black',
    lwd=1.7,
    labels_gp = gpar(fontsize = legend_label_size),
    title_position = "leftcenter-rot"
  ),
  rect_gp = gpar(col = 'transparent',
                 lwd = 0
  ), use_raster = T
)
p_fig4C_activity

#### Figure 4D: A similar heatmap but with log fitness ####

alpha <- 80
beta <- 0.5

## Look at the logw
# Organize as a matrix for a heatmap
data_syn_ratio_bind_wide <- data_syn_ratio_bind %>%
  mutate(total_activity = 0.1*cA + 0.1*cB + 1*cAA + 1*cBB + 1*cAB) %>% 
  mutate(logw = round(log2(beta^(log2(total_activity / alpha)^2)), 5)) %>%
  select(log2_val, binding_diff, logw) %>%
  arrange(log2_val, desc(binding_diff)) %>%
  pivot_wider(names_from = log2_val, values_from = logw)

# Move the first column to the index
first_column <- data_syn_ratio_bind_wide$binding_diff
data_syn_ratio_bind_wide <- as.matrix(data_syn_ratio_bind_wide[, 2:ncol(data_syn_ratio_bind_wide)])
rownames(data_syn_ratio_bind_wide) <- first_column

## Draw heatmaps for effect of synthesis rate

# Define column labels
column_labels = seq(from = -2, to = 1.95, by = 0.05)
new_column_labels = column_labels
for(i in 1:length(column_labels)){
  
  cur_val = column_labels[i]
  new_val = ifelse(mod(round(cur_val*10, 2), 5) == 0, toString(cur_val), '')
  
  new_column_labels[i] = new_val
  
}

# Define row labels
row_labels = rev(seq(from = -3, to = 2.9, by = 0.1))
new_row_labels = row_labels
for(i in 1:length(row_labels)){
  
  cur_val = row_labels[i]
  new_val = ifelse(mod(cur_val*10, 10) == 0, toString(cur_val), '')
  
  new_row_labels[i] = new_val
  
}


p_fig4D_logw <- Heatmap(
  data_syn_ratio_bind_wide,
  cluster_columns = F, cluster_rows = F,
  col = colorRamp2(
    breaks = c(-1, -0.75, -0.5, -0.25, 0),
    colors = cividis(5)
  ),
  show_column_names = T, column_names_side = 'bottom',
  column_names_rot = 90, column_names_centered = T,
  column_names_gp = gpar(fontsize=column_name_size, fontface='bold', align = 'bottom', colour = 'black'),
  column_labels = new_column_labels,
  show_row_names = T, row_names_side = 'left',
  row_names_rot = 90,  row_names_centered = T,
  row_names_gp = gpar(fontsize=row_name_size, fontface = 'bold', align = 'left', colour = 'black'),
  row_labels = new_row_labels,

  width=unit(heatmap_width, 'cm'), height = unit(heatmap_height, 'cm'),
  border = T,
  row_title = expression(paste(
    bold('\u0394'),
    bold(G['bind,HET AB']), 
    bold(' - \u0394'),
    bold(G['bind,HM AA']), 
    bold(' (kcal/mol)'),
    sep = '')),
  column_title = expression(paste(bold('log2('),
                                  bold(s['B']),
                                  bold(' / '),
                                  bold(s['A']), 
                                  bold(')'))),
  row_title_gp = gpar(fontsize=row_title_size, fontface = 'bold'),
  column_title_gp = gpar(fontsize=column_title_size, fontface = 'bold'),
  column_title_side = 'bottom',
  
  
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(
    title = "log2(Fitness)",
    title_gp = gpar(fontsize = legend_title_size),
    legend_height = unit(legend_height, "cm"),
    grid_width = unit(legend_width, "cm"),
    border='black',
    lwd=1.7,
    labels_gp = gpar(fontsize = legend_label_size),
    title_position = "leftcenter-rot"
  ),
  rect_gp = gpar(col = 'transparent',
                 lwd = 0
  ), use_raster = T
)
p_fig4D_logw

#### Figure 4E ####

## Load the data for parametric simulations with gene expression
all_results_parametric_gene_exp_dup <- read_delim('Results_simulations/Results_pdup1/Parametric_simulations/008_simulations_parametric_gene_exp_pdup1_80opt/all_results_all_sims.tsv', delim = '\t')

all_results_parametric_gene_exp_dup %<>% 
  select(-dA, -aA, -dAA, -aAA, -dB, -aB, -dBB, -aBB, -dAB, -t, -Start_aff, -Start_stab)

## Check for duplications
all_results_parametric_dup_check <- all_results_parametric_gene_exp_dup %>% ungroup() %>%
  mutate(notNA = !(is.na(sB))) %>%
  group_by(Replicate, Expression_mut_probs, intHM_param) %>%
  summarise(dup_mut_count = sum(notNA))

all_results_param_dup_check_final <- all_results_parametric_dup_check %>% ungroup() %>%
  mutate(dup_bool = (dup_mut_count > 0)) %>%
  group_by(Expression_mut_probs, intHM_param) %>%
  summarise(rep_dup_count = sum(dup_bool))

all_results_parametric_gene_exp_dup %<>% 
  mutate(fixed_mut = rep(1:200, nrow(all_results_parametric_gene_exp_dup) / 200))

# Calculate the percentage of HET
all_results_parametric_gene_exp_dup %<>% filter(!(is.na(aAB))) %>%
  mutate(Total_complexes = cAA + cAB + cBB + cA + cB) %>% # Add the total concentration of complexes
  mutate(pct_HET = round((cAB *100 / Total_complexes), 2), 
         pct_HMs = round((cAA + cBB) * 100 / Total_complexes, 2), 
         pct_mono = round((cA + cB) * 100 / Total_complexes, 2)
  )

# Label the different replicates depending on where they finish
all_results_parametric_gene_exp_dup_final_points <- all_results_parametric_gene_exp_dup %>% filter(fixed_mut == 200) %>%
  mutate(Outcome = ifelse(pct_HET >= 70, 'HET dominant', 
                          ifelse(pct_HMs >= 70, 'HM dominant', 
                                 ifelse(pct_mono >= 70, 'Monomers', 
                                        ifelse((pct_HET + pct_HMs) >= 70, 'Both HM and HET', 
                                               'Ambiguous')
                                 )
                          )
  )
  )

# Add the info on the final points
all_results_parametric_gene_exp_dup_outcome <- left_join(x = all_results_parametric_gene_exp_dup,
                                                         y = all_results_parametric_gene_exp_dup_final_points %>% 
                                                           select(Replicate, intHM_param, Expression_mut_probs, Outcome), 
                                                         by = c('Replicate' = 'Replicate', 'intHM_param' = 'intHM_param', 
                                                                'Expression_mut_probs' = 'Expression_mut_probs'))

## Figure S08: Trajectories for parametric simulations allowing changes in synthesis rates ##
p_figS08 <- all_results_parametric_gene_exp_dup_outcome %>%
  mutate(intHM_param = str_c('Mean ddG HM = ', intHM_param, sep = ''), 
         Expression_mut_probs = str_c('p_exp = ', Expression_mut_probs, sep = '')) %>%
  mutate(Outcome = factor(Outcome, levels = c('HET dominant', 'Both HM and HET', 'HM dominant', 'Monomers', 'Ambiguous'))) %>%
  ggplot(aes(x = fixed_mut, y = pct_HET, group = Replicate, colour = Outcome)) +
  geom_line() +
  scale_colour_manual(values = c('#4575B4', '#008000', '#FC8D59', 'black', '#B300B3')) +
  facet_grid(intHM_param~Expression_mut_probs) +
  geom_hline(yintercept = 50, linetype ='dashed') +
  xlab('Fixed mutations') + ylab('HET percentage (%)') +
  theme(legend.position = 'top', legend.justification = 'center', 
        panel.grid.major = element_blank(), 
        strip.text = element_text(size = 18), 
        legend.title = element_text(size = 22), 
        legend.text = element_text(size = 20), 
        axis.title = element_text(size = 24), 
        axis.text = element_text(size = 16))
p_figS08
ggsave(p_figS08, device = cairo_pdf, width = 35, height = 28, dpi = 500,
       filename = 'Figures/Supplementary_figures/Appendix_FigS08.Trajectories_parametric_simulations_new.pdf')

#### Convert the results of the simulations above to a heatmap ####
#### x-axis can be probability of changes to expression, y axis can be effect on HM ####

data_fig_4E <- all_results_parametric_gene_exp_dup_outcome %>% 
  filter(fixed_mut == 200) %>% ungroup() %>%
  group_by(Expression_mut_probs, intHM_param) %>%
  summarise(median_HET = median(pct_HET), 
            mean_HET = mean(pct_HET)) %>%
  select(-median_HET)

param_gene_exp_final_points <- data_fig_4E %>%
  pivot_wider(names_from = intHM_param, values_from = mean_HET)

needed_rownames <- param_gene_exp_final_points$Expression_mut_probs
param_gene_exp_final_points <- as.matrix(param_gene_exp_final_points %>% ungroup() %>% select(-Expression_mut_probs))
rownames(param_gene_exp_final_points) <- needed_rownames

p_fig4E <- Heatmap(t(param_gene_exp_final_points),
                   cluster_columns = F, cluster_rows = F, 
                   col = colorRamp2(
                     breaks = c(0, 20, 40, 60, 80, 100),
                     colors = viridis(6)),
                   show_column_names = T, row_names_side = 'left',
                   width=unit(heatmap_width, 'cm'), height = unit(heatmap_height, 'cm'),
                   border = T,
                   column_title = "Probability of mutations affecting gene expression",
                   row_title = expression(paste(
                     bold('Mean \u0394\u0394'), 
                     bold(G['bind,HM']), 
                     sep = ''
                   )),
                   row_title_gp = gpar(fontsize=row_title_size, fontface = 'bold'),
                   column_title_gp = gpar(fontsize=column_title_size, fontface = 'bold'),
                   column_title_side = 'bottom',
                   row_names_rot = 0, 
                   row_names_centered = T,
                   row_names_gp = gpar(fontsize=row_name_size, fontface = 'bold', align = 'left'),
                   column_names_gp = gpar(fontsize=column_name_size, fontface='bold', align = 'left'),
                   show_heatmap_legend = TRUE,
                   heatmap_legend_param = list(
                     at = c(0, 20, 40, 60, 80, 100),
                     title = "Percentage of HET AB (%)", 
                     title_gp = gpar(fontsize = legend_title_size),
                     legend_height = unit(legend_height, "cm"),
                     grid_width = unit(legend_width, "cm"),
                     border='black',
                     lwd=1.7,
                     labels_gp = gpar(fontsize = legend_label_size),
                     title_position = "leftcenter-rot"
                   )
)
p_fig4E


## Save the data
write.table(data_fig_4E, append = F, quote = F, sep = '\t', row.names = F, col.names = T, 
            file = 'Figures/Main_figures/data_fig_4E.tsv')

#### Figure 4F ####

## Load data for structures when allowing changes in gene expression ##
all_results_gene_exp <- read_delim('Results_simulations/Results_pdup1/Simulations_with_structures/008_simulations_gene_expression_pdup1_80opt/all_results_all_sims.tsv', delim = '\t')

## Remove structures
all_results_gene_exp %<>% filter(!(Complex %in% structures_remove))

all_results_gene_exp %<>% 
  select(-dA, -aA, -dAA, -aAA, -dB, -aB, -dBB, -aBB, -dAB, -t)

# Add a column indicating how many mutations have fixed up to that point
all_results_gene_exp %<>% group_by(Complex, Replicate, Start_stab, Start_aff, Expression_mut_probs) %>%
  mutate(fixed_mut = row_number())

## Look at the final points to check the outcomes of simulations
all_data_final_points_gene_exp<- all_results_gene_exp %>% ungroup() %>%
  mutate(Complex = toupper(Complex)) %>%
  filter(!(is.na(sB))) %>% # Seems like in some replicates the duplication did not fix
  mutate(Total_complexes = cAA + cAB + cBB + cA + cB) %>% # Add the total concentration of complexes
  filter(fixed_mut == 200, Start_stab == -5, Start_aff == -10) %>%
  mutate(outcome = ifelse(cAB / Total_complexes >= 0.7, 'HET dominant', 
                          ifelse( (cAA + cBB) / Total_complexes >= 0.7, 'HM dominant', 
                                  ifelse( (cA + cB) / Total_complexes >= 0.7, 'Monomers', 
                                          ifelse( (cAA + cAB + cBB) / Total_complexes >= 0.7, 'Both HM and HET', 
                                                  'Ambiguous')
                                  )
                          )
  )
  )

#### Figure S09: Check how much the two paralogs diverged in expression ####

## Look at how much the synthesis rates of the paralogs diverged

gene_exp_synth_rates_param <- all_data_final_points_gene_exp %>%
  ## When classifying by the subunit with higher synthesis rate
  mutate(s_abundant_HM = ifelse(sA >= sB, sA, sB), 
         s_low_HM = ifelse(sA > sB, sB, sA))

## Get an idea of how different synthesis rates are throughout
test_exp_diff <- gene_exp_synth_rates_param %>% filter(Expression_mut_probs == 0.9) %>%
  mutate(diff_synthesis = s_abundant_HM - s_low_HM)
summary(test_exp_diff$diff_synthesis)

## Look at it from the perspective of a log2-fold difference
test_exp_diff %<>% mutate(log2FC = log2(s_abundant_HM / s_low_HM))
summary(test_exp_diff$log2FC)

p_figS09 <- gene_exp_synth_rates_param %>% ungroup() %>% 
  pivot_longer(cols = c('s_abundant_HM', 's_low_HM'), names_to = 'Checked_paralog', values_to = 'Synthesis_rate') %>%
  mutate(Checked_paralog = case_when(
    Checked_paralog == 's_abundant_HM' ~ 'More abundant subunit',
    Checked_paralog == 's_low_HM' ~ 'Less abundant subunit'
  )) %>%
  ggplot(aes(x = as.factor(Expression_mut_probs), y = Synthesis_rate, fill = Checked_paralog)) +
  geom_boxplot() +
  geom_hline(yintercept = 100, linetype = 'dashed') +
  theme(legend.position = 'top', legend.justification = 'center', 
        panel.grid.major = element_blank()
  ) +
  xlab('Probability of mutations affecting synthesis rate') + ylab('Synthesis rate') +
  labs(fill = '') +
  ylim(50, 200)
p_figS09
ggsave(p_figS09, device = cairo_pdf, width = 18, height = 10, units = 'cm', dpi = 300, 
       filename = 'Figures/Supplementary_figures/Appendix_FigS09.Divergence_paralogs.pdf')

## Continue working with the simulations with data on gene expression
all_data_summary_gene_exp <- all_results_gene_exp %>% 
  ungroup() %>%
  filter(!(is.na(sB))) %>% # Seems like in some replicates the duplication did not fix
  mutate(Total_complexes = cAA + cAB + cBB + cA + cB) %>% # Add the total concentration of complexes
  mutate(pct_HET = round((cAB *100 / Total_complexes), 2), 
         pct_HMs = round((cAA + cBB) * 100 / Total_complexes, 2), 
         pct_mono = round((cA + cB) * 100 / Total_complexes, 2)) %>%
  group_by(Start_stab, Start_aff, Expression_mut_probs, Complex, fixed_mut) %>%
  summarise(
    med_pct_HET = median(pct_HET), 
    mean_pct_HET = mean(pct_HET), 
    sem_pct_HET = sd(pct_HET) / sqrt(n()), 
    
    mean_pct_HMs = mean(pct_HMs), 
    sem_pct_HMs = sd(pct_HMs) / sqrt(n()),
    
    mean_pct_mono = mean(pct_mono), 
    sem_pct_mono = sd(pct_mono) / sqrt(n())
  )

# Get the order of the figures
all_data_summary_final_points_gene_exp <- all_data_summary_gene_exp %>% 
  filter(fixed_mut == 200, Start_stab == -5, Start_aff == -10) %>%
  arrange(desc(mean_pct_HET))

all_data_summary_final_points_gene_exp2 <- all_data_summary_gene_exp %>% 
  filter(fixed_mut == 200, Start_stab == -5, Start_aff == -10, Expression_mut_probs == 0.1) %>%
  arrange(desc(mean_pct_HET))

# Relevel according to the order from the other plot
all_data_summary_gene_exp %<>% 
  mutate(Complex = toupper(Complex)) %>%
  mutate(Complex =  factor(Complex, levels = unique(all_data_summary_final_points$Complex)))

## Use the stacked barplot visualization
all_data_summary_final_points_gene_exp %<>% arrange(desc(mean_pct_HET))

all_data_summary_final_points_gene_exp %<>%
  mutate(Complex = toupper(Complex)) %>%
  mutate(Complex = factor(Complex, levels = unique(all_data_summary_final_points_gene_exp$Complex)))

all_data_summary_final_points_long_gene_exp <- all_data_summary_final_points_gene_exp %>% ungroup() %>%
  pivot_longer(cols = c(mean_pct_HET, mean_pct_HMs, mean_pct_mono), names_to = 'assembly', values_to = 'mean_value')

#### Figure showing the percentages of each molecular species ####

gene_exp_dup_check <- all_results_gene_exp %>%
  filter(fixed_mut == 200) %>%
  mutate(dup_check = !(is.na(sB))) %>%
  select(Complex, Replicate, Expression_mut_probs, dup_check)

gene_exp_all_final_points <- all_results_gene_exp %>%
  filter(fixed_mut == 200) %>% 
  mutate(Total_complexes = cAA + cAB + cBB + cA + cB) %>% # Add the total concentration of complexes
  mutate(pct_HET = round((cAB *100 / Total_complexes), 2), 
         pct_HMs = round((cAA + cBB) * 100 / Total_complexes, 2), 
         pct_mono = round((cA + cB) * 100 / Total_complexes, 2))

gene_exp_long_filtered_final_points <- inner_join(x = gene_exp_all_final_points, 
                                                  y = gene_exp_dup_check, 
                                                  by = c('Complex' = 'Complex', 'Replicate' = 'Replicate', 
                                                         'Expression_mut_probs' = 'Expression_mut_probs'))

avg_outcome_reps_gene_exp <- gene_exp_long_filtered_final_points %>% ungroup() %>%
  filter(dup_check) %>%
  mutate(pct_abundant_HM = ifelse(cAA > cBB, cAA * 100 / Total_complexes, cBB * 100 / Total_complexes), 
         pct_low_HM = ifelse(cAA > cBB, cBB * 100 / Total_complexes, cAA * 100 / Total_complexes), 
         pct_abundant_mono = ifelse(cAA > cBB, cA * 100 / Total_complexes, cB * 100 / Total_complexes), 
         pct_low_mono = ifelse(cAA > cBB, cB * 100 / Total_complexes, cA * 100 / Total_complexes)) %>%
  group_by(Complex, Expression_mut_probs) %>%
  summarise(mean_pct_HET = mean(pct_HET), 
            mean_pct_abundant_HM = mean(pct_abundant_HM),
            mean_pct_low_HM = mean(pct_low_HM),
            mean_pct_abundant_mono = mean(pct_abundant_mono), 
            mean_pct_low_mono = mean(pct_low_mono),
            num_reps = n())

avg_outcome_reps_gene_exp <- avg_outcome_reps_gene_exp %>%
  pivot_longer(cols = c(mean_pct_HET, mean_pct_abundant_HM, mean_pct_low_HM,
                        mean_pct_abundant_mono, mean_pct_low_mono),
               names_to = 'assembly', values_to = 'mean_value') %>%
  mutate(assembly = ifelse(assembly == 'mean_pct_HET', 'HET',
                           ifelse(assembly == 'mean_pct_abundant_HM', 'High HM',
                                  ifelse(assembly == 'mean_pct_abundant_mono', 'Monomer of high HM',
                                         ifelse(assembly == 'mean_pct_low_HM', 'Low HM',
                                                ifelse(assembly == 'mean_pct_low_mono', 'Monomer of low HM', NA)
                                         )
                                  )
                           )
  )
  ) %>%
  filter(assembly %in% c('HET', 'High HM', 'Low HM'))


p_fig4F_all_lines <- avg_outcome_reps_gene_exp %>%
  ggplot(aes(x = as.factor(Expression_mut_probs), y = mean_value)) +
  geom_hline(yintercept = 50, linetype = 'dashed') +
  geom_line(aes(group = interaction(assembly, Complex), colour = assembly), 
            alpha = 0.05, linewidth = 0.5) +
  stat_summary(fun = 'mean', geom = 'line', aes(group = assembly, colour = assembly), 
               linewidth = 1) +
  scale_fill_manual(values = c('#4575B4', '#FC8D59', '#ff6666')) +
  scale_colour_manual(values = c('#4575B4', '#FC8D59', '#ff6666')) +
  xlab('Probability of sampled mutations affecting synthesis rate') +
  ylab('Percentage of complex') +
  labs(fill = '', colour = '') +
  theme(legend.position = 'top', legend.justification = 'center', 
        panel.grid.major = element_blank(), 
        axis.title = element_text(size = 9), 
        axis.text = element_text(size = 8)) +
  ylim(0, 100)
p_fig4F_all_lines

#### Draw figure 4  ####

p_fig4 <- plot_grid(draw_CHeatmap(p_fig4A), draw_CHeatmap(p_fig4B), 
                    draw_CHeatmap(p_fig4C_activity), draw_CHeatmap(p_fig4D_logw), 
                    draw_CHeatmap(p_fig4E), p_fig4F_all_lines,
                    nrow = 3, labels = c('A', 'B', 'C', 'D', 'E', 'F'), 
                    label_size = panel_label_size, label_fontface = 'bold'
                    )
ggsave(p_fig4, device = cairo_pdf, width = 20, height = 20, units = 'cm', dpi = 300, 
       filename = 'Figures/Main_figures/4.Fig4.pdf')
## Cartoons for the two types of simulations were added manually

#### Figure EV4: Similar to figure 4 but with alpha = 60 ####

## Figure EV4A: Fitness when alpha = 60 ##

alpha <-60
beta <- 0.5

## Look at the logw
# Organize as a matrix for a heatmap
data_syn_ratio_bind_wide <- data_syn_ratio_bind %>%
  mutate(total_activity = 0.1*cA + 0.1*cB + 1*cAA + 1*cBB + 1*cAB) %>% 
  mutate(logw = round(log2(beta^(log2(total_activity / alpha)^2)), 5)) %>%
  select(log2_val, binding_diff, logw) %>%
  arrange(log2_val, desc(binding_diff)) %>%
  pivot_wider(names_from = log2_val, values_from = logw)

# Move the first column to the index
first_column <- data_syn_ratio_bind_wide$binding_diff
data_syn_ratio_bind_wide <- as.matrix(data_syn_ratio_bind_wide[, 2:ncol(data_syn_ratio_bind_wide)])
rownames(data_syn_ratio_bind_wide) <- first_column

## Draw heatmaps for effect of synthesis rate

# Define column labels
column_labels = seq(from = -2, to = 1.95, by = 0.05)
new_column_labels = column_labels
for(i in 1:length(column_labels)){
  
  cur_val = column_labels[i]
  new_val = ifelse(mod(round(cur_val*10, 2), 5) == 0, toString(cur_val), '')
  
  new_column_labels[i] = new_val
  
}

# Define row labels
row_labels = rev(seq(from = -3, to = 2.9, by = 0.1))
new_row_labels = row_labels
for(i in 1:length(row_labels)){
  
  cur_val = row_labels[i]
  new_val = ifelse(mod(cur_val*10, 10) == 0, toString(cur_val), '')
  
  new_row_labels[i] = new_val
  
}


p_figEV4A_logw <- Heatmap(
  data_syn_ratio_bind_wide,
  cluster_columns = F, cluster_rows = F,
  col = colorRamp2(
    breaks = c(-1, -0.75, -0.5, -0.25, 0),
    colors = cividis(5)
  ),
  show_column_names = T, column_names_side = 'bottom',
  column_names_rot = 90, column_names_centered = T,
  column_names_gp = gpar(fontsize=column_name_size, fontface='bold', align = 'bottom', colour = 'black'),
  column_labels = new_column_labels,

  show_row_names = T, row_names_side = 'left',
  row_names_rot = 90,  row_names_centered = T,
  row_names_gp = gpar(fontsize=row_name_size, fontface = 'bold', align = 'left', colour = 'black'),
  row_labels = new_row_labels,

  width=unit(heatmap_width, 'cm'), height = unit(heatmap_height, 'cm'),
  border = T,
  row_title = expression(paste(
    bold('\u0394'),
    bold(G['bind,HET AB']), 
    bold(' - \u0394'),
    bold(G['bind,HM AA']), 
    bold(' (kcal/mol)'),
    sep = '')),
  column_title = expression(paste(bold('log2('),
                                  bold(s['B']),
                                  bold(' / '),
                                  bold(s['A']), 
                                  bold(')'))),
  row_title_gp = gpar(fontsize=row_title_size, fontface = 'bold'),
  column_title_gp = gpar(fontsize=column_title_size, fontface = 'bold'),
  column_title_side = 'bottom',
  
  
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(
    title = "log2(Fitness)",
    title_gp = gpar(fontsize = legend_title_size),
    legend_height = unit(legend_height, "cm"),
    grid_width = unit(legend_width, "cm"),
    border='black',
    lwd=1.7,
    labels_gp = gpar(fontsize = legend_label_size),
    title_position = "leftcenter-rot"
  ),
  rect_gp = gpar(col = 'transparent',
                 lwd = 0
  ), use_raster = T
)
p_figEV4A_logw

#### Figure EV4B: Parametric simulations with alpha = 60 ####

## Load data for alpha = 60
all_results_parametric_gene_exp_dup_60opt <- read_delim('Results_simulations/Results_pdup1/Parametric_simulations/008_simulations_parametric_gene_exp_pdup1_60opt/all_results_all_sims.tsv', delim = '\t')

all_results_parametric_gene_exp_dup_60opt %<>% 
  select(-dA, -aA, -dAA, -aAA, -dB, -aB, -dBB, -aBB, -dAB, -t)

## Let's see in how many replicates and conditions we do see the duplication
all_results_parametric_dup_check_60opt <- all_results_parametric_gene_exp_dup_60opt %>% ungroup() %>%
  mutate(notNA = !(is.na(sB))) %>%
  group_by(Replicate, Expression_mut_probs, intHM_param) %>%
  summarise(dup_mut_count = sum(notNA))

all_results_param_dup_check_final_60opt <- all_results_parametric_dup_check_60opt %>% ungroup() %>%
  mutate(dup_bool = (dup_mut_count > 0)) %>%
  group_by(Expression_mut_probs, intHM_param) %>%
  summarise(rep_dup_count = sum(dup_bool))

## The only parameter that changed was the mean of the distribution of mutational effects for HM
table(all_results_parametric_gene_exp_dup_60opt$intHM_param)
table(all_results_parametric_gene_exp_dup_60opt$intHET_param)
table(all_results_parametric_gene_exp_dup_60opt$stab_param)
table(all_results_parametric_gene_exp_dup_60opt$Expression_mut_probs)

all_results_parametric_gene_exp_dup_60opt %<>% 
  mutate(fixed_mut = rep(1:200, nrow(all_results_parametric_gene_exp_dup_60opt) / 200))

# Calculate the percentage of HET
all_results_parametric_gene_exp_dup_60opt %<>% filter(!(is.na(aAB))) %>%
  mutate(Total_complexes = cAA + cAB + cBB + cA + cB) %>% # Add the total concentration of complexes
  mutate(pct_HET = round((cAB *100 / Total_complexes), 2), 
         pct_HMs = round((cAA + cBB) * 100 / Total_complexes, 2), 
         pct_mono = round((cA + cB) * 100 / Total_complexes, 2)
  )

# Label the different replicates depending on where they finish
all_results_parametric_gene_exp_dup_final_points_60opt <- all_results_parametric_gene_exp_dup_60opt %>%
  filter(fixed_mut == 200) %>%
  mutate(Outcome = ifelse(pct_HET >= 70, 'HET dominant', 
                          ifelse(pct_HMs >= 70, 'HM dominant', 
                                 ifelse(pct_mono >= 70, 'Monomers', 
                                        ifelse((pct_HET + pct_HMs) >= 70, 'Both HM and HET', 
                                               'Ambiguous')
                                 )
                          )
  )
  )

# Add the info on the final points
all_results_parametric_gene_exp_dup_outcome_60opt <- left_join(x = all_results_parametric_gene_exp_dup_60opt,
                                                               y = all_results_parametric_gene_exp_dup_final_points_60opt %>% 
                                                                 select(Replicate, intHM_param, Expression_mut_probs, Outcome), 
                                                               by = c('Replicate' = 'Replicate', 'intHM_param' = 'intHM_param', 
                                                                      'Expression_mut_probs' = 'Expression_mut_probs'))


#### Convert the results of the simulations above to a heatmap ####
#### x-axis can be probability of changes to expression, y axis can be effect on HM ####

data_fig_EV4B <- all_results_parametric_gene_exp_dup_outcome_60opt %>% 
  filter(fixed_mut == 200) %>% ungroup() %>%
  group_by(Expression_mut_probs, intHM_param) %>%
  summarise(median_HET = median(pct_HET), 
            mean_HET = mean(pct_HET)) %>%
  select(-median_HET)

param_gene_exp_final_points_60opt <- data_fig_S10B %>%
  pivot_wider(names_from = intHM_param, values_from = mean_HET)

needed_rownames <- param_gene_exp_final_points_60opt$Expression_mut_probs
param_gene_exp_final_points_60opt <- as.matrix(param_gene_exp_final_points_60opt %>%
                                                 ungroup() %>% select(-Expression_mut_probs))
rownames(param_gene_exp_final_points_60opt) <- needed_rownames


p_figEV4B <- Heatmap(t(param_gene_exp_final_points_60opt),
                    cluster_columns = F, cluster_rows = F, 
                    col = colorRamp2(
                      breaks = c(0, 20, 40, 60, 80, 100),
                      # colors = magma(6)),
                      colors = viridis(6)),
                    show_column_names = T, row_names_side = 'left',
                    width=unit(heatmap_width, 'cm'), height = unit(heatmap_height, 'cm'),
                    border = T,
                    column_title = "Probability of mutations affecting gene expression",
                    # row_title = "Mean ddG binding on the homomer",
                    row_title = expression(paste(
                      bold('Mean \u0394\u0394'), 
                      bold(G['bind,HM']), 
                      sep = ''
                    )),
                    row_title_gp = gpar(fontsize=row_title_size, fontface = 'bold'),
                    column_title_gp = gpar(fontsize=column_title_size, fontface = 'bold'),
                    column_title_side = 'bottom',
                    row_names_rot = 0, 
                    row_names_centered = T,
                    row_names_gp = gpar(fontsize=row_name_size, fontface = 'bold', align = 'left'),
                    column_names_gp = gpar(fontsize=column_name_size, fontface='bold', align = 'left'),
                    show_heatmap_legend = TRUE,
                    heatmap_legend_param = list(
                      at = c(0, 20, 40, 60, 80, 100),
                      title = "Percentage of HET AB (%)", 
                      title_gp = gpar(fontsize = legend_title_size),
                      legend_height = unit(legend_height, "cm"),
                      grid_width = unit(legend_width, "cm"),
                      border='black',
                      lwd=1.7,
                      labels_gp = gpar(fontsize = legend_label_size),
                      title_position = "leftcenter-rot"
                    ), 
)
p_figEV4B

#### Figure EV4C: General trajectories for simulations with structures with changes in synth. rates (opt60) ####

## Load the data for simulations with alpha = 60
all_results_gene_exp_60opt <- read_delim('Results_simulations/Results_pdup1/Simulations_with_structures/008_simulations_gene_expression_pdup1_60opt/all_results_all_sims.tsv', delim = '\t')

all_results_gene_exp_60opt %<>% 
  select(-dA, -aA, -dAA, -aAA, -dB, -aB, -dBB, -aBB, -dAB, -aAB, -t)

# Add a column indicating how many mutations have fixed up to that point
all_results_gene_exp_60opt %<>% group_by(Complex, Replicate, Start_stab, Start_aff, Expression_mut_probs) %>%
  mutate(fixed_mut = row_number())

## Look at the final points to check the outcomes of simulations
all_data_final_points_gene_exp_60opt <- all_results_gene_exp_60opt %>% ungroup() %>%
  mutate(Complex = toupper(Complex)) %>%
  filter(!(is.na(sB))) %>% # Seems like in some replicates the duplication did not fix
  mutate(Total_complexes = cAA + cAB + cBB + cA + cB) %>% # Add the total concentration of complexes
  filter(fixed_mut == 200, Start_stab == -5, Start_aff == -10) %>%
  mutate(outcome = ifelse(cAB / Total_complexes >= 0.7, 'HET dominant',
                          ifelse( (cAA + cBB) / Total_complexes >= 0.7, 'HM dominant',
                                  ifelse( (cA + cB) / Total_complexes >= 0.7, 'Monomers',
                                          ifelse( (cAA + cAB + cBB) / Total_complexes >= 0.7, 'Both HM and HET',
                                                  'Ambiguous')
                                  )
                          )
  )
  )


# Add a column indicating how many mutations have fixed up to that point
gene_exp_dup_check_60opt <- all_results_gene_exp_60opt %>%
  filter(fixed_mut == 200) %>%
  mutate(dup_check = !(is.na(sB))) %>%
  select(Complex, Replicate, Expression_mut_probs, dup_check)

gene_exp_all_final_points_60opt <- all_results_gene_exp_60opt %>%
  filter(fixed_mut == 200) %>% 
  mutate(Total_complexes = cAA + cAB + cBB + cA + cB) %>% # Add the total concentration of complexes
  mutate(pct_HET = round((cAB *100 / Total_complexes), 2), 
         pct_HMs = round((cAA + cBB) * 100 / Total_complexes, 2), 
         pct_mono = round((cA + cB) * 100 / Total_complexes, 2))

gene_exp_long_filtered_final_points_60opt <- inner_join(x = gene_exp_all_final_points_60opt, 
                                                        y = gene_exp_dup_check_60opt, 
                                                        by = c('Complex' = 'Complex', 'Replicate' = 'Replicate', 
                                                               'Expression_mut_probs' = 'Expression_mut_probs'))

avg_outcome_reps_gene_exp_60opt <- gene_exp_long_filtered_final_points_60opt %>% ungroup() %>%
  filter(dup_check) %>%
  mutate(pct_abundant_HM = ifelse(cAA > cBB, cAA * 100 / Total_complexes, cBB * 100 / Total_complexes), 
         pct_low_HM = ifelse(cAA > cBB, cBB * 100 / Total_complexes, cAA * 100 / Total_complexes), 
         pct_abundant_mono = ifelse(cAA > cBB, cA * 100 / Total_complexes, cB * 100 / Total_complexes), 
         pct_low_mono = ifelse(cAA > cBB, cB * 100 / Total_complexes, cA * 100 / Total_complexes)) %>%
  group_by(Complex, Expression_mut_probs) %>%
  summarise(mean_pct_HET = mean(pct_HET), 
            mean_pct_abundant_HM = mean(pct_abundant_HM),
            mean_pct_low_HM = mean(pct_low_HM),
            mean_pct_abundant_mono = mean(pct_abundant_mono), 
            mean_pct_low_mono = mean(pct_low_mono),
            num_reps = n())

avg_outcome_reps_gene_exp_60opt <- avg_outcome_reps_gene_exp_60opt %>%
  pivot_longer(cols = c(mean_pct_HET, mean_pct_abundant_HM, mean_pct_low_HM,
                        mean_pct_abundant_mono, mean_pct_low_mono),
               names_to = 'assembly', values_to = 'mean_value') %>%
  mutate(assembly = ifelse(assembly == 'mean_pct_HET', 'HET',
                           ifelse(assembly == 'mean_pct_abundant_HM', 'High HM',
                                  ifelse(assembly == 'mean_pct_abundant_mono', 'Monomer of high HM',
                                         ifelse(assembly == 'mean_pct_low_HM', 'Low HM',
                                                ifelse(assembly == 'mean_pct_low_mono', 'Monomer of low HM', NA)
                                         )
                                  )
                           )
  )
  ) %>%
  filter(assembly %in% c('HET', 'High HM', 'Low HM'))

p_figEV4C_all_lines <- avg_outcome_reps_gene_exp_60opt %>%
  ggplot(aes(x = as.factor(Expression_mut_probs), y = mean_value)) +
  geom_hline(yintercept = 50, linetype = 'dashed') +
  geom_line(aes(group = interaction(assembly, Complex), colour = assembly), 
            alpha = 0.1) +
  stat_summary(fun = 'mean', geom = 'line', aes(group = assembly, colour = assembly), 
               linewidth = 2.5) +
  scale_fill_manual(values = c('#4575B4', '#FC8D59', '#ff6666')) +
  scale_colour_manual(values = c('#4575B4', '#FC8D59', '#ff6666')) +
  xlab('Probability of sampled mutations affecting synthesis rate') +
  ylab('Percentage of complex') +
  labs(fill = '', colour = '') +
  theme(legend.position = 'top', legend.justification = 'center', 
        panel.grid.major = element_blank()) +
  ylim(0, 100)
p_figEV4C_all_lines

p_figEV4_top <- plot_grid(draw_CHeatmap(p_figEV4A_logw), draw_CHeatmap(p_figEV4B),
                         labels = c('A', 'B'), label_fontface = 'bold', 
                         label_size = panel_label_size, ncol = 2)
p_figEV4 <- plot_grid(p_figEV4_top, p_figEV4C_all_lines, labels = c('', 'C'), 
                     label_fontface = 'bold', label_size = panel_label_size, nrow = 2, 
                     rel_heights = c(0.8, 1))
ggsave(p_figEV4, device = cairo_pdf, width = 18, height = 15, units = 'cm',  dpi = 300, 
       filename = 'Figures/ExpandedView_figures/FigEV4.Simulations_gene_exp_alpha60.pdf')
## Cartoons were added manually

#### Figure 5 ####

#### Figures with the space of solutions with different specific activities ####

#### Figure 5A ####

## Load the data
data_space_activity <- read_delim('Results_solution_space/solutionSpace_hetBias_diff_binding.tsv', 
                                  delim = '\t')

## Need a heatmap with the percentage of HET
# Organize as a matrix for a heatmap
data_space_activity_wide <- data_space_activity %>%
  select(hetBias, binding_diff, pct_AB) %>%
  arrange(hetBias, desc(binding_diff)) %>%
  pivot_wider(names_from = hetBias, values_from = pct_AB)

# Move the first column to the index
first_column <- data_space_activity_wide$binding_diff
data_space_activity_wide <- as.matrix(data_space_activity_wide[, 2:ncol(data_space_activity_wide)])
rownames(data_space_activity_wide) <- first_column

## Draw heatmaps for effect of specific activity

# Define column labels
column_labels = seq(from = -80, to = 80, by = 1)
new_column_labels = column_labels
for(i in 1:length(column_labels)){
  
  cur_val = column_labels[i]
  new_val = ifelse(mod(round(cur_val, 2), 20) == 0, toString(cur_val), '')
  
  new_column_labels[i] = new_val
  
}

# Define row labels
row_labels = rev(seq(from = -3, to = 2.9, by = 0.1))
new_row_labels = row_labels
for(i in 1:length(row_labels)){
  
  cur_val = row_labels[i]
  new_val = ifelse(mod(cur_val*10, 10) == 0, toString(cur_val), '')
  
  new_row_labels[i] = new_val
  
}

p_fig5A <- Heatmap(
  data_space_activity_wide,
  cluster_columns = F, cluster_rows = F,
  col = colorRamp2(
    breaks = c(0, 20, 40, 60, 80, 100),
    colors = viridis(6)),
  show_column_names = T, column_names_side = 'bottom',
  column_names_rot = 90, column_names_centered = T,
  column_names_gp = gpar(fontsize=column_name_size, fontface='bold', align = 'bottom', colour = 'black'),
  column_labels = new_column_labels,

  show_row_names = T, row_names_side = 'left',
  row_names_rot = 90,  row_names_centered = T,
  row_names_gp = gpar(fontsize=row_name_size, fontface = 'bold', align = 'left', colour = 'black'),
  row_labels = new_row_labels,

  width=unit(heatmap_width, 'cm'), height = unit(heatmap_height, 'cm'),
  border = T,
  row_title = expression(paste(
    bold('\u0394'),
    bold(G['bind,HET AB']), 
    bold(' - \u0394'),
    bold(G['bind,HM AA']), 
    bold(' (kcal/mol)'),
    sep = '')),
  column_title = 'HET activity bias (%)',
  row_title_gp = gpar(fontsize=row_title_size, fontface = 'bold'),
  column_title_gp = gpar(fontsize=column_title_size, fontface = 'bold'),
  column_title_side = 'bottom',
  
  
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(
    at = c(0, 20, 40, 60, 80, 100),
    title = "Percentage of HET AB (%)",
    title_gp = gpar(fontsize = legend_title_size),
    legend_height = unit(legend_height, "cm"),
    grid_width = unit(legend_width, "cm"),
    border='black',
    lwd=1.7,
    labels_gp = gpar(fontsize = legend_label_size),
    title_position = "leftcenter-rot"
  ),
  rect_gp = gpar(col = 'transparent',
                 lwd = 0
  ), use_raster = T
)
p_fig5A

## Figure 5B: Similar to figure 6A but with HM BB ##

## Need a heatmap with the percentage of HET
# Organize as a matrix for a heatmap
data_space_activity_wide <- data_space_activity %>%
  mutate(pct_BB = 100 * cBB / (cA + cB + cAA + cAB + cBB)) %>%
  select(hetBias, binding_diff, pct_BB) %>%
  arrange(hetBias, desc(binding_diff)) %>%
  pivot_wider(names_from = hetBias, values_from = pct_BB)

# Move the first column to the index
first_column <- data_space_activity_wide$binding_diff
data_space_activity_wide <- as.matrix(data_space_activity_wide[, 2:ncol(data_space_activity_wide)])
rownames(data_space_activity_wide) <- first_column

# Define column labels
column_labels = seq(from = -80, to = 80, by = 1)
new_column_labels = column_labels
for(i in 1:length(column_labels)){
  
  cur_val = column_labels[i]
  new_val = ifelse(mod(round(cur_val, 2), 20) == 0, toString(cur_val), '')
  
  new_column_labels[i] = new_val
  
}

# Define row labels
row_labels = rev(seq(from = -3, to = 2.9, by = 0.1))
new_row_labels = row_labels
for(i in 1:length(row_labels)){
  
  cur_val = row_labels[i]
  new_val = ifelse(mod(cur_val*10, 10) == 0, toString(cur_val), '')
  
  new_row_labels[i] = new_val
  
}

p_fig5B <- Heatmap(
  data_space_activity_wide,
  cluster_columns = F, cluster_rows = F,
  col = colorRamp2(
    breaks = c(0, 20, 40, 60, 80, 100),
    colors = viridis(6)),
  show_column_names = T, column_names_side = 'bottom',
  column_names_rot = 90, column_names_centered = T,
  column_names_gp = gpar(fontsize=column_name_size, fontface='bold', align = 'bottom', colour = 'black'),
  column_labels = new_column_labels,

  show_row_names = T, row_names_side = 'left',
  row_names_rot = 90,  row_names_centered = T,
  row_names_gp = gpar(fontsize=row_name_size, fontface = 'bold', align = 'left', colour = 'black'),
  row_labels = new_row_labels,

  width=unit(heatmap_width, 'cm'), height = unit(heatmap_height, 'cm'),
  border = T,
  row_title = expression(paste(
    bold('\u0394'),
    bold(G['bind,HET AB']), 
    bold(' - \u0394'),
    bold(G['bind,HM AA']), 
    bold(' (kcal/mol)'),
    sep = '')),
  column_title = 'HET activity bias (%)',
  row_title_gp = gpar(fontsize=row_title_size, fontface = 'bold'),
  column_title_gp = gpar(fontsize=column_title_size, fontface = 'bold'),
  column_title_side = 'bottom',
  
  
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(
    at = c(0, 20, 40, 60, 80, 100),
    title = "Percentage of HM BB (%)",
    title_gp = gpar(fontsize = legend_title_size),
    legend_height = unit(legend_height, "cm"),
    grid_width = unit(legend_width, "cm"),
    border='black',
    lwd=1.7,
    labels_gp = gpar(fontsize = legend_label_size),
    title_position = "leftcenter-rot"
  ),
  rect_gp = gpar(col = 'transparent',
                 lwd = 0
  ), use_raster = T
)
p_fig5B

## Figure 5C: A heatmap of the total activity in the solution space ##

## Look at the total activity
# Organize as a matrix for a heatmap
data_space_activity_wide <- data_space_activity %>%
  select(hetBias, binding_diff, total_activity) %>%
  arrange(hetBias, desc(binding_diff)) %>%
  pivot_wider(names_from = hetBias, values_from = total_activity)

# Move the first column to the index
first_column <- data_space_activity_wide$binding_diff
data_space_activity_wide <- as.matrix(data_space_activity_wide[, 2:ncol(data_space_activity_wide)])
rownames(data_space_activity_wide) <- first_column

# Define column labels
column_labels = seq(from = -80, to = 80, by = 1)
new_column_labels = column_labels
for(i in 1:length(column_labels)){
  
  cur_val = column_labels[i]
  new_val = ifelse(mod(round(cur_val, 2), 20) == 0, toString(cur_val), '')
  
  new_column_labels[i] = new_val
  
}

# Define row labels
row_labels = rev(seq(from = -3, to = 2.9, by = 0.1))
new_row_labels = row_labels
for(i in 1:length(row_labels)){
  
  cur_val = row_labels[i]
  new_val = ifelse(mod(cur_val*10, 10) == 0, toString(cur_val), '')
  
  new_row_labels[i] = new_val
  
}

p_fig5C <- Heatmap(
  data_space_activity_wide,
  cluster_columns = F, cluster_rows = F,
  col = colorRamp2(
    breaks = c(0, 20, 40, 60, 80),
    colors = magma(5)
  ),
  show_column_names = T, column_names_side = 'bottom',
  column_names_rot = 90, column_names_centered = T,
  column_names_gp = gpar(fontsize=column_name_size, fontface='bold', align = 'bottom', colour = 'black'),
  column_labels = new_column_labels,

  show_row_names = T, row_names_side = 'left',
  row_names_rot = 90,  row_names_centered = T,
  row_names_gp = gpar(fontsize=row_name_size, fontface = 'bold', align = 'left', colour = 'black'),
  row_labels = new_row_labels,

  width=unit(heatmap_width, 'cm'), height = unit(heatmap_height, 'cm'),
  border = T,
  row_title = expression(paste(
    bold('\u0394'),
    bold(G['bind,HET AB']), 
    bold(' - \u0394'),
    bold(G['bind,HM AA']), 
    bold(' (kcal/mol)'),
    sep = '')),
  column_title = 'HET activity bias (%)',
  row_title_gp = gpar(fontsize=row_title_size, fontface = 'bold'),
  column_title_gp = gpar(fontsize=column_title_size, fontface = 'bold'),
  column_title_side = 'bottom',
  
  
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(
    title = "Total activity",
    title_gp = gpar(fontsize = legend_title_size),
    legend_height = unit(legend_height, "cm"),
    grid_width = unit(legend_width, "cm"),
    border='black',
    lwd=1.7,
    labels_gp = gpar(fontsize = legend_label_size),
    title_position = "leftcenter-rot"
  ),
  rect_gp = gpar(col = 'transparent',
                 lwd = 0
  ), use_raster = T
)
p_fig5C

#### Figure 5D: A similar heatmap but with log fitness ####

alpha <- 80
beta <- 0.5

## Look at the logw
# Organize as a matrix for a heatmap
data_space_activity_wide <- data_space_activity %>%
  mutate(logw = round(log2(beta^(log2(total_activity / alpha)^2)), 5)) %>%
  select(hetBias, binding_diff, logw) %>%
  arrange(hetBias, desc(binding_diff)) %>%
  pivot_wider(names_from = hetBias, values_from = logw)

# Move the first column to the index
first_column <- data_space_activity_wide$binding_diff
data_space_activity_wide <- as.matrix(data_space_activity_wide[, 2:ncol(data_space_activity_wide)])
rownames(data_space_activity_wide) <- first_column

# Define column labels
column_labels = seq(from = -80, to = 80, by = 1)
new_column_labels = column_labels
for(i in 1:length(column_labels)){
  
  cur_val = column_labels[i]
  new_val = ifelse(mod(round(cur_val, 2), 20) == 0, toString(cur_val), '')
  
  new_column_labels[i] = new_val
  
}

# Define row labels
row_labels = rev(seq(from = -3, to = 2.9, by = 0.1))
new_row_labels = row_labels
for(i in 1:length(row_labels)){
  
  cur_val = row_labels[i]
  new_val = ifelse(mod(cur_val*10, 10) == 0, toString(cur_val), '')
  
  new_row_labels[i] = new_val
  
}


p_fig5D_logw <- Heatmap(
  data_space_activity_wide,
  cluster_columns = F, cluster_rows = F,
  col = colorRamp2(
    breaks = c(-1, -0.75, -0.5, -0.25, 0),
    colors = cividis(5)
  ),
  show_column_names = T, column_names_side = 'bottom',
  column_names_rot = 90, column_names_centered = T,
  column_names_gp = gpar(fontsize=column_name_size, fontface='bold', align = 'bottom', colour = 'black'),
  column_labels = new_column_labels,

  show_row_names = T, row_names_side = 'left',
  row_names_rot = 90,  row_names_centered = T,
  row_names_gp = gpar(fontsize=row_name_size, fontface = 'bold', align = 'left', colour = 'black'),
  row_labels = new_row_labels,

  width=unit(heatmap_width, 'cm'), height = unit(heatmap_height, 'cm'),
  border = T,
  row_title = expression(paste(
    bold('\u0394'),
    bold(G['bind,HET AB']), 
    bold(' - \u0394'),
    bold(G['bind,HM AA']), 
    bold(' (kcal/mol)'),
    sep = '')),
  column_title = 'HET activity bias (%)',
  row_title_gp = gpar(fontsize=row_title_size, fontface = 'bold'),
  column_title_gp = gpar(fontsize=column_title_size, fontface = 'bold'),
  column_title_side = 'bottom',
  
  
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(
    title = "log2(Fitness)",
    title_gp = gpar(fontsize = legend_title_size),
    legend_height = unit(legend_height, "cm"),
    grid_width = unit(legend_width, "cm"),
    border='black',
    lwd=1.7,
    labels_gp = gpar(fontsize = legend_label_size),
    title_position = "leftcenter-rot"
  ),
  rect_gp = gpar(col = 'transparent',
                 lwd = 0
  ), use_raster = T
)
p_fig5D_logw

#### Figure 5E: Overview of the parametric simulations with biased activity ####

## Load data from the run with pdup = 1, alpha = 80
all_results_parametric_biased <- read_delim('Results_simulations/Results_pdup1/Parametric_simulations/008_simulations_parametric_hetBias_pdup1_80opt/all_results_all_sims.tsv', delim = '\t')

all_results_parametric_biased %<>% 
  select(-dA, -dAA, -dB, -dBB, -dAB, -t, -Start_aff, -Start_stab)


all_results_parametric_biased %<>% 
  mutate(aAB = as.numeric(aAB), 
         aAA = as.numeric(aAA)) %>%
  mutate(HETbias = 100*(aAB - aAA)) %>%
  mutate(HETbias = as.numeric(HETbias)) %>%
  mutate(fixed_mut = rep(1:200, nrow(all_results_parametric_biased) / 200), 
         total_activity = aA*cA + aB*cB + aAA*cAA + aAB*cAB + aBB*cBB) 

# Calculate the percentage of HET
all_results_parametric_biased %<>% filter(!(is.na(aAB))) %>%
  mutate(Total_complexes = cAA + cAB + cBB + cA + cB) %>% # Add the total concentration of complexes
  mutate(pct_HET = round((cAB *100 / Total_complexes), 2), 
         pct_HMs = round((cAA + cBB) * 100 / Total_complexes, 2), 
         pct_mono = round((cA + cB) * 100 / Total_complexes, 2)
  )

# Label the different replicates depending on where they finish
all_results_parametric_biased_final_points <- all_results_parametric_biased %>% 
  filter(fixed_mut == 200) %>%
  mutate(Outcome = ifelse(pct_HET >= 70, 'HET dominant', 
                          ifelse(pct_HMs >= 70, 'HM dominant', 
                                 ifelse(pct_mono >= 70, 'Monomers', 
                                        ifelse((pct_HET + pct_HMs) >= 70, 'Both HM and HET', 
                                               'Ambiguous')
                                 )
                          )
  )
  )

# Add the info on the final points
all_results_parametric_biased_outcome <- left_join(x = all_results_parametric_biased,
                                                   y = all_results_parametric_biased_final_points %>% 
                                                     select(Replicate, intHM_param, HETbias, Outcome), 
                                                   by = c('Replicate' = 'Replicate', 
                                                          'intHM_param' = 'intHM_param', 
                                                          'HETbias' = 'HETbias'))


#### Figure S10: Individual trajectories ####

p_figS10 <- all_results_parametric_biased_outcome %>%
  mutate(intHM_param = str_c('Mean ddG HM = ', intHM_param, sep = ''), 
         HETbias = str_c('HETbias = ', HETbias, sep = '')) %>%
  mutate(HETbias = factor(HETbias, 
                          levels = c('HETbias = -80', 'HETbias = -70', 'HETbias = -60',
                                     'HETbias = -50', 'HETbias = -40', 'HETbias = -30',
                                     'HETbias = -20', 'HETbias = -10', 'HETbias = 0',
                                     'HETbias = 10', 'HETbias = 20', 'HETbias = 30',
                                     'HETbias = 40', 'HETbias = 50', 'HETbias = 60',
                                     'HETbias = 70', 'HETbias = 80'))) %>%
  mutate(Outcome = factor(Outcome, levels = c('HET dominant', 'Both HM and HET', 'HM dominant', 'Monomers', 'Ambiguous'))) %>%
  ggplot(aes(x = fixed_mut, y = pct_HET, group = Replicate, colour = Outcome)) +
  geom_line() +
  scale_colour_manual(values = c('#4575B4', '#008000', '#FC8D59', 'black', '#B300B3')) +
  facet_grid(intHM_param~HETbias) +
  xlab('Fixed mutations') + ylab('HET percentage (%)') +
  theme(legend.position = 'top', legend.justification = 'center', 
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 18, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 20), 
        legend.title = element_text(size = 20), legend.text = element_text(size = 18), 
        strip.text = element_text(size = 18)
        ) +
  geom_hline(yintercept = 50, linetype = 'dashed')
p_figS10
ggsave(p_figS10, device = cairo_pdf, width = 45, height = 28, dpi = 300,
       filename = 'Figures/Supplementary_figures/Appendix_FigS10.Trajectories_parametric_biased_simulations_80opt.pdf')

param_biased_final_points <- all_results_parametric_biased_outcome %>% 
  filter(fixed_mut == 200) %>% ungroup() %>%
  group_by(HETbias, intHM_param) %>%
  arrange(HETbias, intHM_param) %>%
  summarise(
    mean_HET = mean(pct_HET)) %>%
  pivot_wider(names_from = intHM_param, values_from = mean_HET)

needed_rownames <- param_biased_final_points$HETbias
param_biased_final_points <- as.matrix(param_biased_final_points %>% ungroup() %>% select(-HETbias))
rownames(param_biased_final_points) <- needed_rownames

p_fig5E <- Heatmap(t(param_biased_final_points),
                   cluster_columns = F, cluster_rows = F, 
                   col = colorRamp2(
                     breaks = c(0, 20, 40, 60, 80, 100),
                     colors = viridis(6)
                   ),
                   show_column_names = T, row_names_side = 'left',
                   width=unit(heatmap_width, 'cm'), height = unit(heatmap_height, 'cm'),
                   border = T,
                   column_title = "HET activity bias (%)",
                   row_title = expression(paste(
                     bold('Mean \u0394\u0394'), 
                     bold(G['bind,HM']), 
                     sep = ''
                   )),
                   row_title_gp = gpar(fontsize=row_title_size, fontface = 'bold'),
                   column_title_gp = gpar(fontsize=column_title_size, fontface = 'bold'),
                   column_title_side = 'bottom',
                   row_names_rot = 0, 
                   row_names_centered = T,
                   row_names_gp = gpar(fontsize=row_name_size, fontface = 'bold', align = 'left'),
                   column_names_gp = gpar(fontsize=column_name_size, fontface='bold', align = 'left'),
                   show_heatmap_legend = TRUE,
                   heatmap_legend_param = list(
                     at = c(0, 20, 40, 60, 80, 100),
                     title = "Percentage of HET AB (%)", 
                     title_gp = gpar(fontsize = legend_title_size),
                     legend_height = unit(legend_height, "cm"),
                     grid_width = unit(legend_width, "cm"),
                     border='black',
                     lwd=1.7,
                     labels_gp = gpar(fontsize = legend_label_size),
                     title_position = "leftcenter-rot"
                   ), 
)
p_fig5E

#### Figure 5F ####

## Simulations with different specific activities for HET and HM

all_results_biased %<>% 
  mutate(fixed_mut = rep(1:200, nrow(all_results_biased) / 200)) %>%
  filter(aAA != 0.1, aAB != 0.1) %>%
  filter(Complex %in% all_results$Complex)

all_biased_data_summary <- all_results_biased %>% 
  ungroup() %>% filter(!(is.na(aAB))) %>%
  group_by(Complex, aAA, aAB, aBB, fixed_mut) %>%
  mutate(Total_complexes = cAA + cAB + cBB + cA + cB) %>% # Add the total concentration of complexes
  mutate(pct_HET = round((cAB *100 / Total_complexes), 2), 
         pct_HMs = round(100*(cAA + cBB) / Total_complexes, 2)
  ) %>%
  summarise(
    aAA = mean(aAA), aAB = mean(aAB),
    med_pct_HET = median(pct_HET), 
    mean_pct_HET = mean(pct_HET), 
    mean_pct_HMs = mean(pct_HMs),
    sem_pct_HET = sd(pct_HET) / sqrt(n())
  )

all_biased_data_summary_final_points <- all_biased_data_summary %>% 
  filter(fixed_mut == 200, aAA == 1, aAB == 1) %>% 
  arrange(desc(mean_pct_HET))

# Relevel according to the means from the simulation with aAA = 1 and aAB = 1
all_biased_data_summary %<>% 
  mutate(Complex = toupper(Complex)) %>%
  ## To have the same order for the structures as in the previous figures
  mutate(Complex = factor(Complex, levels = unique(all_data_summary_final_points$Complex)))

all_results_biased_final_points <- all_results_biased %>% ungroup() %>% 
  filter(!(is.na(aAB)), !(Complex %in% structures_remove), fixed_mut == 200) %>%
  group_by(Complex, aAA, aAB, aBB, fixed_mut) %>%
  mutate(Total_complexes = cAA + cAB + cBB + cA + cB) %>% # Add the total concentration of complexes
  mutate(pct_HET = cAB * 100 / Total_complexes, 
         pct_HMs = (cAA + cBB) * 100 / Total_complexes, 
         pct_monomers = (cA + cB) * 100 / Total_complexes, 
         HET_activity_bias = case_when(
           aAA < 1 ~ 100 - aAA*100, # When AA is less active, this is a positive bias for HET
           aAB < 1 ~ -(100 - aAB*100), # When AB is less active, this is a negative bias for HET (against HET)
           TRUE ~ 0
         )
  )

avg_outcome_reps_biased <- all_results_biased_final_points %>% ungroup() %>%
  mutate(pct_abundant_HM = ifelse(cAA > cBB, cAA * 100 / Total_complexes, cBB * 100 / Total_complexes), 
         pct_low_HM = ifelse(cAA > cBB, cBB * 100 / Total_complexes, cAA * 100 / Total_complexes), 
         pct_abundant_mono = ifelse(cAA > cBB, cA * 100 / Total_complexes, cB * 100 / Total_complexes), 
         pct_low_mono = ifelse(cAA > cBB, cB * 100 / Total_complexes, cA * 100 / Total_complexes)) %>%
  group_by(Complex, HET_activity_bias) %>%
  summarise(mean_pct_HET = mean(pct_HET), 
            mean_pct_abundant_HM = mean(pct_abundant_HM),
            mean_pct_low_HM = mean(pct_low_HM),
            mean_pct_abundant_mono = mean(pct_abundant_mono), 
            mean_pct_low_mono = mean(pct_low_mono),
            num_reps = n())

p_fig5F_all_lines <- avg_outcome_reps_biased %>%
  pivot_longer(cols = c(mean_pct_HET, mean_pct_abundant_HM, mean_pct_low_HM, 
                        mean_pct_abundant_mono, mean_pct_low_mono), 
               names_to = 'assembly', values_to = 'mean_value') %>%
  mutate(assembly = ifelse(assembly == 'mean_pct_HET', 'HET', 
                           ifelse(assembly == 'mean_pct_abundant_HM', 'High HM', 
                                  ifelse(assembly == 'mean_pct_abundant_mono', 'Monomer of high HM',
                                         ifelse(assembly == 'mean_pct_low_HM', 'Low HM',
                                                ifelse(assembly == 'mean_pct_low_mono', 'Monomer of low HM', NA)
                                         )
                                  )
                           )
  )
  ) %>%
  filter(assembly %in% c('HET', 'High HM', 'Low HM')) %>%
  ggplot(aes(x = as.factor(HET_activity_bias), y = mean_value, fill = assembly)) +
  geom_line(aes(group = interaction(assembly, Complex), colour = assembly), 
            alpha = 0.05, linewidth = 0.5) +
  geom_hline(yintercept = 50, linetype = 'dashed') +
  stat_summary(geom = 'line', fun = 'mean', aes(group = assembly, colour = assembly), 
               linewidth = 1) +
  scale_fill_manual(values = c('#4575B4', '#FC8D59', '#ff6666')) +
  scale_colour_manual(values = c('#4575B4', '#FC8D59', '#ff6666')) +
  xlab('HET activity bias (%)') +
  ylab('Percentage of complex') +
  labs(fill = '', colour = '') +
  theme(legend.position = 'top', legend.justification = 'center', 
        panel.grid.major = element_blank(), 
        axis.title = element_text(size = 9), 
        axis.text = element_text(size = 8)) +
  ylim(0, 100)
p_fig5F_all_lines

#### Draw figure 5 ####

p_fig5 <-  plot_grid(draw_CHeatmap(p_fig5A), draw_CHeatmap(p_fig5B),
                     draw_CHeatmap(p_fig5C), draw_CHeatmap(p_fig5D_logw),
                     draw_CHeatmap(p_fig5E), p_fig5F_all_lines,
                     nrow = 3, labels = c('A', 'B', 'C', 'D', 'E', 'F'), 
                     label_size = panel_label_size, label_fontface = 'bold')

ggsave(p_fig5, device = cairo_pdf, width = 20, height = 20, units = 'cm', dpi = 300, 
       filename = 'Figures/Main_figures/5.Fig5.pdf')

#### Figure EV5: Similar to figure 5 but with alpha = 60 ####

## Figure EV5A ##

alpha <- 60
beta <- 0.5

## Load the data
data_space_activity <- read_delim('Results_solution_space/solutionSpace_hetBias_diff_binding.tsv', 
                                  delim = '\t')

## Look at the logw
# Organize as a matrix for a heatmap
data_space_activity_wide <- data_space_activity %>%
  mutate(logw = round(log2(beta^(log2(total_activity / alpha)^2)), 5)) %>%
  select(hetBias, binding_diff, logw) %>%
  arrange(hetBias, desc(binding_diff)) %>%
  pivot_wider(names_from = hetBias, values_from = logw)

# Move the first column to the index
first_column <- data_space_activity_wide$binding_diff
data_space_activity_wide <- as.matrix(data_space_activity_wide[, 2:ncol(data_space_activity_wide)])
rownames(data_space_activity_wide) <- first_column

## Draw heatmaps for effect of different specific activities

# Define column labels
column_labels = seq(from = -80, to = 80, by = 1)
new_column_labels = column_labels
for(i in 1:length(column_labels)){
  
  cur_val = column_labels[i]
  new_val = ifelse(mod(round(cur_val, 2), 20) == 0, toString(cur_val), '')
  
  new_column_labels[i] = new_val
  
}

# Define row labels
row_labels = rev(seq(from = -3, to = 2.9, by = 0.1))
new_row_labels = row_labels
for(i in 1:length(row_labels)){
  
  cur_val = row_labels[i]
  new_val = ifelse(mod(cur_val*10, 10) == 0, toString(cur_val), '')
  
  new_row_labels[i] = new_val
  
}

p_figEV5A_logw <- Heatmap(
  data_space_activity_wide,
  cluster_columns = F, cluster_rows = F,
  col = colorRamp2(
    breaks = c(-1, -0.75, -0.5, -0.25, 0),
    colors = cividis(5)
  ),
  show_column_names = T, column_names_side = 'bottom',
  column_names_rot = 90, column_names_centered = T,
  column_names_gp = gpar(fontsize=column_name_size, fontface='bold', align = 'bottom', colour = 'black'),
  column_labels = new_column_labels,

  show_row_names = T, row_names_side = 'left',
  row_names_rot = 90,  row_names_centered = T,
  row_names_gp = gpar(fontsize=row_name_size, fontface = 'bold', align = 'left', colour = 'black'),
  row_labels = new_row_labels,

  width=unit(heatmap_width, 'cm'), height = unit(heatmap_height, 'cm'),
  border = T,
  row_title = expression(paste(
    bold('\u0394'),
    bold(G['bind,HET AB']), 
    bold(' - \u0394'),
    bold(G['bind,HM AA']), 
    bold(' (kcal/mol)'),
    sep = '')),
  column_title = 'HET activity bias (%)',
  row_title_gp = gpar(fontsize=row_title_size, fontface = 'bold'),
  column_title_gp = gpar(fontsize=column_title_size, fontface = 'bold'),
  column_title_side = 'bottom',
  
  
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(
    title = "log2(Fitness)",
    title_gp = gpar(fontsize = legend_title_size),
    legend_height = unit(legend_height, "cm"),
    legend_width = unit(2, "cm"),
    border='black',
    lwd=1.7,
    labels_gp = gpar(fontsize = legend_label_size),
    title_position = "leftcenter-rot"
  ),
  rect_gp = gpar(col = 'transparent',
                 lwd = 0
  ), use_raster = T
)
p_figEV5A_logw

## Figure EV5B: Load parametric simulations with the optimum at 60 ##

## Data from the run with pdup = 1, alpha = 60
all_results_parametric_biased_opt60 <- read_delim('Results_simulations/Results_pdup1/Parametric_simulations/008_simulations_parametric_hetBias_pdup1_60opt/all_results_all_sims.tsv', delim = '\t')

all_results_parametric_biased_opt60 %<>% 
  select(-sA, -dA, -dAA, -sB, -dB, -dBB, -dAB)


all_results_parametric_biased_opt60 %<>% 
  mutate(aAB = as.numeric(aAB), 
         aAA = as.numeric(aAA)) %>%
  mutate(HETbias = 100*(aAB - aAA)) %>%
  mutate(HETbias = as.numeric(HETbias)) %>%
  mutate(fixed_mut = rep(1:200, nrow(all_results_parametric_biased_opt60) / 200), 
         total_activity = aA*cA + aB*cB + aAA*cAA + aAB*cAB + aBB*cBB) 

### Some tests to look at the data for each simulation ####

# Calculate the percentage of HET
all_results_parametric_biased_opt60 %<>% filter(!(is.na(aAB))) %>%
  mutate(Total_complexes = cAA + cAB + cBB + cA + cB) %>% # Add the total concentration of complexes
  mutate(pct_HET = round((cAB *100 / Total_complexes), 2), 
         pct_HMs = round((cAA + cBB) * 100 / Total_complexes, 2), 
         pct_mono = round((cA + cB) * 100 / Total_complexes, 2)
  )

# Label the different replicates depending on where they finish
all_results_parametric_biased_final_points_opt60 <- all_results_parametric_biased_opt60 %>% 
  filter(fixed_mut == 200) %>%
  mutate(Outcome = ifelse(pct_HET >= 70, 'HET dominant', 
                          ifelse(pct_HMs >= 70, 'HM dominant', 
                                 ifelse(pct_mono >= 70, 'Monomers', 
                                        ifelse((pct_HET + pct_HMs) >= 70, 'Both HM and HET', 
                                               'Ambiguous')
                                 )
                          )
  )
  )

# Add the info on the final points
all_results_parametric_biased_outcome_opt60 <- left_join(x = all_results_parametric_biased_opt60,
                                                         y = all_results_parametric_biased_final_points_opt60 %>% 
                                                           select(Replicate, intHM_param, HETbias, Outcome), 
                                                         by = c('Replicate' = 'Replicate', 
                                                                'intHM_param' = 'intHM_param', 
                                                                'HETbias' = 'HETbias'))


#### Convert the results of the simulations above to a heatmap ####
#### x-axis can be probability of changes to expression, y axis can be effect on HM ####

param_biased_final_points_opt60 <- all_results_parametric_biased_outcome_opt60 %>% 
  filter(fixed_mut == 200) %>% ungroup() %>%
  group_by(HETbias, intHM_param) %>%
  arrange(HETbias, intHM_param) %>%
  summarise(
    mean_HET = mean(pct_HET)) %>%
  pivot_wider(names_from = intHM_param, values_from = mean_HET)
needed_rownames <- param_biased_final_points_opt60$HETbias
param_biased_final_points_opt60 <- as.matrix(param_biased_final_points_opt60 %>%
                                               ungroup() %>% select(-HETbias))
rownames(param_biased_final_points_opt60) <- needed_rownames

p_figEV5B <- Heatmap(t(param_biased_final_points_opt60),
                     cluster_columns = F, cluster_rows = F, 
                     col = colorRamp2(
                       breaks = c(0, 20, 40, 60, 80, 100),
                       colors = viridis(6)
                     ),
                     show_column_names = T, row_names_side = 'left',
                     width=unit(heatmap_width, 'cm'), height = unit(heatmap_height, 'cm'),
                     border = T,
                     column_title = "HET activity bias (%)",
                     row_title = expression(paste(
                       bold('Mean \u0394\u0394'), 
                       bold(G['bind,HM']), 
                       sep = ''
                     )),
                     row_title_gp = gpar(fontsize=row_title_size, fontface = 'bold'),
                     column_title_gp = gpar(fontsize=column_title_size, fontface = 'bold'),
                     column_title_side = 'bottom',
                     row_names_rot = 0, 
                     row_names_centered = T,
                     row_names_gp = gpar(fontsize=row_name_size, fontface = 'bold', align = 'left'),
                     column_names_gp = gpar(fontsize=column_name_size, fontface='bold', align = 'left'),
                     show_heatmap_legend = TRUE,
                     heatmap_legend_param = list(
                       at = c(0, 20, 40, 60, 80, 100),
                       title = "Percentage of HET AB (%)", 
                       title_gp = gpar(fontsize = legend_title_size),
                       legend_height = unit(legend_height, "cm"),
                       legend_width = unit(2, "cm"),
                       border='black',
                       lwd=1.7,
                       labels_gp = gpar(fontsize = legend_label_size),
                       title_position = "leftcenter-rot"
                     ), 
)
p_figEV5B

#### Figure EV5C-D: Simulations with structures and activity biases ####

all_results_biased_60opt <- read_delim('Results_simulations/Results_pdup1/Simulations_with_structures/008_biased_simulations_60opt/all_results_all_sims.tsv', delim = '\t')

all_results_biased_60opt %<>% 
  select(-sA, -dA, -aA, -dAA, -sB, -dB, -aB, -dBB, -dAB)

all_results_biased_60opt %<>%
  mutate(fixed_mut = rep(1:200, nrow(all_results_biased_60opt) / 200)) %>%
  filter(aAA != 0.1, aAB != 0.1) %>%
  filter(Complex %in% all_results$Complex)


all_biased_data_summary_60opt <- all_results_biased_60opt %>% 
  ungroup() %>% filter(!(is.na(aAB))) %>%
  group_by(Complex, aAA, aAB, aBB, fixed_mut) %>%
  mutate(Total_complexes = cAA + cAB + cBB + cA + cB) %>% # Add the total concentration of complexes
  mutate(pct_HET = round((cAB *100 / Total_complexes), 2), 
         pct_HMs = round(100*(cAA + cBB) / Total_complexes, 2)
  ) %>%
  summarise(
    aAA = mean(aAA), aAB = mean(aAB),
    med_pct_HET = median(pct_HET), 
    mean_pct_HET = mean(pct_HET), 
    mean_pct_HMs = mean(pct_HMs),
    sem_pct_HET = sd(pct_HET) / sqrt(n())
  )

all_biased_data_summary_final_points_60opt <- all_biased_data_summary_60opt %>% 
  filter(fixed_mut == 200, aAA == 1, aAB == 1) %>% 
  arrange(desc(mean_pct_HET))


# Relevel according to the means from the simulation with aAA = 1 and aAB = 1
all_biased_data_summary_60opt %<>% 
  mutate(Complex = toupper(Complex)) %>%
  ## To have the same order for the structures as in the previous figures
  mutate(Complex = factor(Complex, levels = unique(all_data_summary_final_points$Complex)))

all_results_biased_final_points_60opt <- all_results_biased_60opt %>% ungroup() %>% 
  filter(!(is.na(aAB)), !(Complex %in% structures_remove), fixed_mut == 200) %>%
  group_by(Complex, aAA, aAB, aBB, fixed_mut) %>%
  mutate(Total_complexes = cAA + cAB + cBB + cA + cB) %>% # Add the total concentration of complexes
  mutate(pct_HET = cAB * 100 / Total_complexes, 
         pct_HMs = (cAA + cBB) * 100 / Total_complexes, 
         pct_monomers = (cA + cB) * 100 / Total_complexes, 
         HET_activity_bias = case_when(
           aAA < 1 ~ 100 - aAA*100, # When AA is less active, this is a positive bias for HET
           aAB < 1 ~ -(100 - aAB*100), # When AB is less active, this is a negative bias for HET (against HET)
           TRUE ~ 0
         )
  )

#### Figure EV5C: Show facets of values of delta(G_fold) ####

all_results_biased_final_points_new_60opt <- all_results_biased_final_points_60opt %>% ungroup() %>%
  mutate(pct_both = pct_HET + pct_HMs) %>%
  mutate(Outcome = case_when(
    pct_HET >= 70 ~ 'HET dominant', 
    pct_HMs >= 70 ~ 'HM dominant', 
    pct_both >= 70 ~ 'Both HM and HET',
    pct_monomers >= 70 ~ 'Monomers',
    TRUE ~ 'Ambiguous'
  )) %>%
  select(Complex, Replicate, HET_activity_bias, curr_stab_A, curr_stab_B, Outcome,
         pct_HET, pct_HMs, pct_monomers) %>%
  mutate(more_stable_subunit = case_when(
    curr_stab_A <= curr_stab_B ~ curr_stab_A,
    curr_stab_B < curr_stab_A ~ curr_stab_B
  ),
  less_stable_subunit = case_when(
    curr_stab_A > curr_stab_B ~ curr_stab_A,
    curr_stab_B >= curr_stab_A ~ curr_stab_B
  )
  )

p_figEV5C <- all_results_biased_final_points_new_60opt %>%
  pivot_longer(cols = c('more_stable_subunit', 'less_stable_subunit'), names_to = 'subunit',
               values_to = 'dG_fold') %>%
  mutate(Outcome = factor(Outcome, levels = c('HET dominant', 'Both HM and HET', 'HM dominant',
                                              'Ambiguous', 'Monomers'))) %>%
  mutate(subunit = case_when(
    subunit == 'more_stable_subunit' ~ 'More stable subunit', 
    subunit == 'less_stable_subunit' ~ 'Less stable subunit'
  )) %>%
  filter(Outcome %in% c('HET dominant', 'Both HM and HET', 'HM dominant')) %>%
  ggplot(aes(x = as.factor(HET_activity_bias), y = dG_fold, fill = subunit)) +
  geom_point(aes(colour = subunit), position = position_jitterdodge(jitter.width = 0.3)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  ylim(-10, 5) +
  facet_wrap(~Outcome) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  theme(panel.grid.major =  element_blank(),
        legend.position = 'top', legend.justification = 'center', 
        axis.text.x = element_text(angle = 90, size = 8), 
        axis.text.y = element_text(size = 8), 
        axis.title = element_text(size = 9), 
        strip.text = element_text(size = 9)
        ) +
  xlab('HET activity bias (%)') + 
  labs(y = expression(paste(bold('\u0394'),
                            bold(G[fold]), 
                            bold(' (kcal/mol)'), sep = '')
  ), fill = '', colour = ''
  )
p_figEV5C

#### Figure EV5D: Percentage of heteromers throughout ####

avg_outcome_reps_biased_60opt <- all_results_biased_final_points_60opt %>% ungroup() %>%
  mutate(pct_abundant_HM = ifelse(cAA > cBB, cAA * 100 / Total_complexes, cBB * 100 / Total_complexes), 
         pct_low_HM = ifelse(cAA > cBB, cBB * 100 / Total_complexes, cAA * 100 / Total_complexes), 
         pct_abundant_mono = ifelse(cAA > cBB, cA * 100 / Total_complexes, cB * 100 / Total_complexes), 
         pct_low_mono = ifelse(cAA > cBB, cB * 100 / Total_complexes, cA * 100 / Total_complexes)) %>%
  group_by(Complex, HET_activity_bias) %>%
  summarise(mean_pct_HET = mean(pct_HET), 
            mean_pct_abundant_HM = mean(pct_abundant_HM),
            mean_pct_low_HM = mean(pct_low_HM),
            mean_pct_abundant_mono = mean(pct_abundant_mono), 
            mean_pct_low_mono = mean(pct_low_mono),
            num_reps = n())

## Draw figure EV5D
p_figEV5D_all_lines <- avg_outcome_reps_biased_60opt %>%
  pivot_longer(cols = c(mean_pct_HET, mean_pct_abundant_HM, mean_pct_low_HM, 
                        mean_pct_abundant_mono, mean_pct_low_mono), 
               names_to = 'assembly', values_to = 'mean_value') %>%
  mutate(assembly = ifelse(assembly == 'mean_pct_HET', 'HET', 
                           ifelse(assembly == 'mean_pct_abundant_HM', 'High HM', 
                                  ifelse(assembly == 'mean_pct_abundant_mono', 'Monomer of high HM',
                                         ifelse(assembly == 'mean_pct_low_HM', 'Low HM',
                                                ifelse(assembly == 'mean_pct_low_mono', 'Monomer of low HM', NA)
                                         )
                                  )
                           )
  )
  ) %>%
  filter(assembly %in% c('HET', 'High HM', 'Low HM')) %>%
  ggplot(aes(x = as.factor(HET_activity_bias), y = mean_value, fill = assembly)) +
  geom_line(aes(group = interaction(assembly, Complex), colour = assembly), 
            alpha = 0.2) +
  geom_hline(yintercept = 50, linetype = 'dashed') +
  stat_summary(geom = 'line', fun = 'mean', aes(group = assembly, colour = assembly), 
               linewidth = 2) +
  scale_fill_manual(values = c('#4575B4', '#FC8D59', '#ff6666')) +
  scale_colour_manual(values = c('#4575B4', '#FC8D59', '#ff6666')) +
  xlab('HET activity bias (%)') +
  ylab('Percentage of complex') +
  labs(fill = '', colour = '') +
  theme(legend.position = 'top', legend.justification = 'center', 
        panel.grid.major = element_blank(), 
        axis.title = element_text(size = 9), 
        axis.text = element_text(size = 8)
        ) +
  ylim(0, 100)
p_figEV5D_all_lines

## Draw figure EV5 ##

p_figEV5_top <- plot_grid(draw_CHeatmap(p_figEV5A_logw), draw_CHeatmap(p_figEV5B), ncol = 2, 
                          labels = c('A', 'B'), label_size = panel_label_size, label_fontface = 'bold')

p_figEV5 <- plot_grid(p_figEV5_top, p_figEV5C, p_figEV5D_all_lines, nrow = 3, 
                      labels = c('', 'C', 'D'), label_size = panel_label_size, 
                      label_fontface = 'bold')

ggsave(p_figEV5, device = cairo_pdf, width = 20, height = 20, units = 'cm', dpi = 300, 
       filename = 'Figures/ExpandedView_figures/FigEV5.Biased_sims_60opt_new.pdf')
## Cartoons were added manually

