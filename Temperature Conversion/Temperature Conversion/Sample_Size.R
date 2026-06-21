
sample_size <- function(study_type, outcome, groups){
    
    switch(study_type,                  # Level one
           "observational" = { 
               switch(outcome,          # Level two
                      "continuous" = {
                          switch(groups, # Level three
                                 "single" = "One sample mean formula",
                                 "multiple" = "Two sample mean formula",
                                 "Not in database"
                          )
                      },
                      
                      "binary" = {
                          switch(groups,
                                 "single" = "One sample proportion formula",
                                 "multiple" = "Two sample proportion formula",
                                 "Not in database"
                          )
                      }
                      
                      
               )
               
           },
         
           
           
           "experimental" = {
               switch(outcome,
                      "continuous" = {
                          switch(groups,
                                 
                                 "single" = "One sample RCT formula",
                                 "multiple" = "Two sample RCT formula",
                                 "Not in database"
                          )
                      },
                      
                      "binary" = {
                          switch(groups,
                                 
                                 "single" = "One sample RCT proportion",
                                 "multiple" = "Two sample RCT proportion",
                                 "Not in database"
                          )
                      }
                      
                   
               )
               
           },
           
           
           "Not in database"
           
           
    )
    
    
    
   

}

sample_size("observational", "binary", "multiple")
sample_size("experimental", "continuous", "single")
