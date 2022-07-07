# Function to plot Single patient single regulon pliot
patient_regulon_activity_plot <- function(patient_id, regulon, regulons,patient_exp, label){
  # get regulon genes
  my.regulon.genes <- regulons[regulon][[1]]

  # Get single patient value to combine with cohort for plotting
  patients_exp_single <- patient_exp %>%
    dplyr::select(X1, paste(patient_id))

  # Get TPM values
  patients_exp_data <- patient_exp %>%
    dplyr::select(X1, paste(patient_id)) %>%
    dplyr::filter(X1 %in% my.regulon.genes) %>%
    dplyr::rename("zscore" = paste(patient_id)) %>%
    dplyr::mutate(patient= patient_id)

  # filter cohort expression data for regulon genes
  exp_data_cohort <- cohort_expr %>%
    dplyr::filter(X %in% my.regulon.genes) %>%
    unique()

  ## Add patient data to cohort data
  exp_data_with_patient <- inner_join(exp_data_cohort, patients_exp_single, by=c("X" = "X1") )

  ## melt
  exp_data_cohort_melt <- reshape2::melt(exp_data_with_patient)

  ## Add median and color
  exp_data_final <-  exp_data_cohort_melt %>%
    group_by(variable) %>%
    dplyr::mutate(median=median(value)) %>%
    dplyr::mutate(Q4=quantile(value)["100%"]) %>%
    dplyr::mutate(Color=if_else(variable == patient_id, "red", "gray"))


  ### GENE PLOT
  # plot of cohort expression for regulon genes across patients
  # q <- ggplot(exp_data_cohort_melt, aes(x=reorder(X,value,na.rm = TRUE), y=value, group=X))
  # q <- q + geom_violin(alpha=0.7, show.legend = F, fill="gray")
  # q <- q + geom_point(data=patients_exp_data, aes(x=X1, y=zscore), color="red", inherit.aes = FALSE)
  # q <- q + geom_line(data=patients_exp_data, aes(x=X1, y=zscore, group=2),color="red", inherit.aes = FALSE)
  # q <- q + theme(axis.text.x = element_text(size=7, angle=45))
  # q <- q + labs(x="Gene(s)", y="Gene Expression z-score", title=paste0(patient_id, ", ", "Regulon: ", regulon, " - ",label))
  # print(q)

  ### PATIENT PLOT
  # plot of regulon gene exopression across patients
  p <- ggplot(exp_data_final,
              aes(x=reorder(variable,value,na.rm = TRUE),
                  y=value, group=variable,fill=Color,color=Color)
  )
  p <- p + geom_violin(show.legend = F)
  p <- p + geom_vline(xintercept = patient_id, color="red", linetype="dotted")
  p <- p + geom_point(data=exp_data_final, aes(x=reorder(variable,value,na.rm = TRUE), y=median, group=variable),color="blue", inherit.aes = FALSE)
  p <- p + scale_fill_identity(aesthetics = c("fill","color"))
  p <- p + theme(axis.text.x = element_blank())
  p <- p + labs(x="Patients", y="Gene Expression z-score", title=paste0(patient_id, ", ", "Regulon: ", regulon, " - ",label))
  print(p)

}
