    epi_recommendation <- function(region, setting, disease, evidence){
    
    switch(region,
           
           "africa" = {
               switch(setting,
                      "rural" = {
                          switch(disease,
                                 "malaria" = {
                                     switch(evidence,
                                            "strong" = "Deploy rural malaria campaign with strong evidence base",
                                            "weak" = "Deploy rural malaria campaign with caution",
                                            "Not in database"
                                         
                                     )
                                 },
                                 
                                 
                                 "cholera" = {
                                     switch(evidence,
                                            "strong" = "Deploy rural WASH intervention with strong evidence",
                                            "weak" = "Deploy rural WASH intervention with caution",
                                            "Not in database"
                                            
                                     )
                                 }
                                 
                              
                          )
                      },
                      
                      "urban" = {
                          switch(disease,
                                 "malaria" = {
                                     switch(evidence,
                                            "strong" = "Deploy urban malaria campaign with strong evidence",
                                            "weak" = "Deploy urban malaria campaign with caution",
                                            "Not in database"
                                            
                                     )
                                 },
                                 
                                 
                                 "cholera" = {
                                     switch(evidence,
                                            "strong" = "Deploy urban WASH intervention with strong evidence",
                                            "weak" = "Deploy urban WASH intervention with caution",
                                            "Not in database"
                                            
                                     )
                                 }
                                 
                                 
                          )
                      }
                      
                      
                      
                   
               )
           },
           
           
           
           "europe" = {
               switch(setting,
                      "rural" = {
                          switch(disease,
                                 "malaria" = {
                                     switch(evidence,
                                            "strong" = "Rare case - initiate rural malaria protocol",
                                            "weak" = "Rare case - monitor rural malaria situation"
                                            
                                     )
                                 },
                                 
                                 
                                 "cholera" = {
                                     switch(evidence,
                                            "strong" = "Rare case - initiate rural cholera protocol",
                                            "weak" = "Rare case - monitor rural cholera situation"
                                            
                                     )
                                 }
                                 
                                 
                          )
                      },
                      
                      "urban" = {
                          switch(disease,
                                 "malaria" = {
                                     switch(evidence,
                                            "strong" = "Rare case - initiate urban malaria protocol",
                                            "weak" = "Rare case - monitor urban malaria situation",
                                            "Not in database"
                                            
                                     )
                                 },
                                 
                                 
                                 "cholera" = {
                                     switch(evidence,
                                            "strong" = "Rare case - initiate urban cholera protocol",
                                            "weak" = "Rare case - monitor urban cholera situation",
                                            "Not in database"
                                            
                                     )
                                 }
                                 
                                 
                          )
                      }
                      
                      
                      
                      
               )
           },
           
           
           "Not in database"
           
           
           
        
        
    )
    
    
    
}


epi_recommendation("africa", "rural", "malaria", "strong")



regions   <- c("africa", "africa", "europe", "europe")
settings  <- c("rural",  "urban",  "rural",  "urban")
diseases  <- c("malaria","cholera","malaria","cholera")
evidences <- c("strong", "weak",   "strong", "weak")


sapply(1:8, function(i) epi_recommendation(
    regions[i],
    settings[i],
    diseases[i],
    evidences[i]
))



sap <- sapply(1:2, function(i) epi_recommendation(
           regions[i],
           settings [i],
           diseases [i],
           evidences[i]
       ))



lapply(1:4, function(i) epi_recommendation(
    regions[i], settings[i], diseases[i], evidences[i]
))
