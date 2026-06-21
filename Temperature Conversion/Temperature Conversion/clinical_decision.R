clinical_decision <- function(disease, severity){
    switch(disease,
           
           "malaria" = switch(severity,
                              "mild"   = "Prescribe oral artemether",
                              "severe" = "Admit and give IV artesunate"
           ),
           
           "cholera" = switch(severity,
                              "mild"   = "Oral rehydration therapy",
                              "severe" = "Admit and give IV fluids"
           ),
           
           "tuberculosis" = switch(severity,
                                   "mild"   = "Start standard DOTS therapy",
                                   "severe" = "Admit and start DOTS therapy"
           ),
           
           "Disease not in database"  # default
    )
}

clinical_decision("malaria", "mild")
