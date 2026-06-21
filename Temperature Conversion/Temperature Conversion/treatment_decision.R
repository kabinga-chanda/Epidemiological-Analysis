treatment_decision <- function(disease, severity, age, temperature){
    
    switch(disease,
           
           "malaria" = {
               switch(severity,
                      
                      "severe" = {
                          if(age > 5){
                              if(temperature > 38.5){
                                  "Admit - IV artesunate + cooling"
                              } else {
                                  "Admit - IV artesunate"
                              }
                          } else {
                              "Admit - Pediatric IV artesunate"
                          }
                      },
                      
                      "mild" = "Oral artemether and rest"
               )
           },
           
           "cholera" = {
               switch(severity,
                      
                      "severe" = {
                          if(age > 60){
                              "Admit - IV fluids elderly protocol"
                          } else {
                              "Admit - IV fluids standard protocol"
                          }
                      },
                      
                      "mild" = "Oral rehydration therapy"
               )
           },
           
           "Consult senior clinician"
    )
}