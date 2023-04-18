#' Create Panel Meta Data
#'
#' @param flowset flowSet containing FCS files
#'
#' @return data frame containing channel, antigen, and class information
#' @export
#'
create_panel <- function(flowset){

  param_data <- flowset[[1]]@parameters

  panel_df <- data.frame("channel" = param_data[rownames(param_data)]$name,
                         "antigen" = param_data[rownames(param_data)]$desc,
                         "class" = "none")

  none_params <- c("Time", "SSC-H", "SSC-A", "FSC-H", "FSC-A", "SSC-B-H",
                   "SSC-B-A", "EvID")

  for (i in 1:nrow(panel_df)) {

    if (!is.na(panel_df[i,"antigen"])){

      if (!panel_df[i,"channel"] %in% none_params){

        cat(paste0(rep("\n",4)))
        prompt <- paste(paste0(panel_df[i, 1],":    ", panel_df[i, 2]),
                        "Choose parameter class:",
                        "","s: state","t: type","[]: none",
                        sep="\n")
        answer <- readline(prompt)

        param_class <- switch(answer,
                              s = "state",
                              t = "type",
                              "none"
        )
        panel_df[i, 3] <- param_class
      }

    }


  }

  panel_df

}
