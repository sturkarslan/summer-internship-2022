library(readxl)
library(tidyverse)
library(synergyfinder)
library(ggplot2)
library(pheatmap)
library(grid)

## Concentrations----
btz_concs <- rev(c("0 nM",
                   "2 nM",
                   "3 nM",
                   "3.25 nM",
                   "4 nM",
                   "4.25 nM",
                   "5 nM",
                   "6 nM"))

btz_new_concs <- rev(c("0 nM",
                       "0.5 nM",
                       "0.75 nM",
                       "1 nM",
                       "2 nM",
                       "3 nM",
                       "4 nM",
                       "6 nM"))

flu_concs <- rev(c("max",
                   "min",
                   "0 ug/ml",
                   "5 ug/ml",
                   "10 ug/ml",
                   "15 ug/ml",
                   "20 ug/ml",
                   "25 ug/ml",
                   "30 ug/ml",
                   "35 ug/ml"))

flu_new_concs <- rev(c("max",
                       "min",
                       "0 ug/ml",
                       "2.5 ug/ml",
                       "7.5 ug/ml",
                       "10 ug/ml",
                       "40 ug/ml",
                       "50 ug/ml",
                       "60 ug/ml",
                       "70 ug/ml"))

lut_concs <- rev(c("max",
                   "min",
                   "0 ug/ml",
                   "5 ug/ml",
                   "6.5 ug/ml",
                   "8 ug/ml",
                   "9.5 ug/ml",
                   "11 ug/ml",
                   "12.5 ug/ml",
                   "14 ug/ml"))

lut_new_concs <- rev(c("max",
                   "min",
                   "0 ug/ml",
                   "1 ug/ml",
                   "3 ug/ml",
                   "4 ug/ml",
                   "5 ug/ml",
                   "8 ug/ml",
                   "12.5 ug/ml",
                   "14 ug/ml"))
## Experiements----
assay1_file <- "data/checkerboard/Checkerboard Assay051922.xlsx"
assay2_file <- "data/checkerboard/Checkerboard Assay062322.xlsx"
assay3_file <- "data/checkerboard/Checkerboard Assay070722.xlsx"
assay4_file <- "data/checkerboard/Checkerboard Assay072822.xlsx"

# Combvine experiemtns
assays <- c(assay2_file,assay3_file,assay4_file)


# Data processing----
## Read all experiment results and worksheets into a common fil
all_assays_df <- tibble()

for(assay in assays){
  print(assay)

  if(assay == assay2_file ) {
    exp_id <- "062322"
  }
  if(assay == assay3_file ) {
    exp_id <- "070722"
  }
  if(assay == assay4_file ) {
    exp_id <- "072822"
  }

  # read sample sheet names
  samples <- excel_sheets(assay)

  # loop through each worksheet and combine
  for(sample in samples){
    var_name <- strsplit(sample, split = " ")[[1]][1]
    #print(var_name)
    drug_name <- strsplit(var_name, split = "+", fixed=T)[[1]][1]
    #print(drug_name)
    tmp1 <- read_excel(assay, sheet=sample, col_names = F) #%>%
    #dplyr::select(-"<>")
    if(assay == assay4_file ) {
      rownames(tmp1) <- btz_new_concs
    } else{
      rownames(tmp1) <- btz_concs
    }


    # assign colnames
    if(var_name %in% c("Flu+Bortezomib_MM.1S_only","Flu+Bortezomib_MM.1S+hFOB")){
      if(assay == assay4_file ) {
        colnames(tmp1) <- flu_new_concs
      } else{
        colnames(tmp1) <- flu_concs
      }

    } else{
      if(assay == assay4_file ) {
        colnames(tmp1) <- lut_new_concs
      } else {
        colnames(tmp1) <- lut_concs
      }
    }

    tmp2 <- tmp1 %>%
      dplyr::select(!c(min,max))

    if(assay == assay4_file ) {
      rownames(tmp2) <- btz_new_concs
    } else{
      rownames(tmp2) <- btz_concs
    }

    # for each condition, create a heatmap and write to png
    png(paste0("output/synergy/exp_",exp_id,"_",sample, "_checkerboard.png"))
    setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
    pheatmap(tmp2, cluster_rows = F, cluster_cols = F, main = paste0(sample,"\n","Exp:", exp_id),display_numbers = T,fontsize = 14)
    setHook("grid.newpage", NULL, "replace")
    grid.text(paste0(drug_name), y=-0.07, gp=gpar(fontsize=16))
    grid.text("BTZ", x=-0.07, rot=90, gp=gpar(fontsize=16))
    dev.off()

    assign(var_name, tmp1)

    # format data for synergyfinder

    z1 <- tmp1 %>%
      rownames_to_column("BTZ_conc") %>%
      pivot_longer(cols = !BTZ_conc,
                   names_to = c("conc_c", "conc_c_unit"),
                   names_sep = " ",
                   values_to = "response_value"
      )

    # get the mean of the miniaml # of cells
    min_mean <- z1 %>%
      filter(conc_c == "min") %>%
      pull(response_value) %>%
      mean()

    # get the mean of the maximal # of cells
    max_mean <- z1 %>%
      filter(conc_c == "max") %>%
      pull(response_value) %>%
      mean()

    # get the % viability
    z2 <- z1 %>%
      mutate(response = (response_value* 100)/max_mean) %>%
      filter(conc_c != "min") %>%
      filter(conc_c != "max") %>%
      separate(BTZ_conc, into = c("conc_r", "conc_r_unit"), sep = " ") %>%
      mutate(block_id = var_name) %>%
      mutate(drug_row = paste("BTZ")) %>%
      mutate(drug_col = paste(drug_name)) %>%
      mutate(conc_c = as.numeric(conc_c)) %>%
      mutate(conc_r = as.numeric(conc_r)) %>%
      mutate(experiment = exp_id) %>%
      mutate(condition = sample)

    # run synergy analysis
    synergy_results <- synergy_function(input_df = z2, block = var_name, exp = exp_id)
    #dev.off()

    # combine all the data
    all_assays_df <- bind_rows(all_assays_df, z2)

  }

}



# Synergy calculation function----
synergy_function <- function(input_df, block,exp){
  z <- input_df
  res <- NULL

  # reshape the data
  res <- ReshapeData(
    data = z,
    data_type = "viability",
    impute = TRUE,
    impute_method = NULL,
    noise = TRUE,
    seed = 1)

  # calculate synergy
  res <- CalculateSynergy(
    data = res,
    method = c("ZIP", "HSA", "Bliss", "Loewe"),
    Emin = NA,
    Emax = NA,
    correct_baseline = "all")

  # calculate sensitivity
  # res <- CalculateSensitivity(
  #   data = res,
  #   correct_baseline = "all"
  # )
  #
  PlotDoseResponse(data = res,
                   heatmap_plot_title = paste0(block),
                   block_ids = block,
                   drugs = c(1,2),
                   save_file = T,
                   file_type = "png",
                   file_name = c(paste0("output/synergy/exp_",exp,"_",as.character(block),"_DoseResponse.png")),
                   high_value_color = "red",
                   low_value_color= "blue"
                   )
  dev.off()

  pp1 <- Plot2DrugHeatmap(data = res,
                          plot_title = paste0(block," Bliss Synergy"),
                          plot_block = block,
                          drugs = c(1, 2),
                          plot_value = "Bliss_synergy",
                          dynamic = F,
                          summary_statistic = c("mean","quantile_75"),
                          high_value_color = "red",
                          low_value_color= "blue"
                          )
  print(pp1)

  ggsave(scale = 2, plot = pp1, device = "png",filename = paste0("output/synergy/exp_",exp,"_",block,"_Heatmap_synergy.png"))
  dev.off()

  return(res)
}





