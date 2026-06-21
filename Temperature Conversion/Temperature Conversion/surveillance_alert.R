surveillance_alert <- function(disease, cases, location){
    
    switch(disease,
           
           "measles" = {
               if(cases < 5){
                   "No alert"
               } else {
                   switch(location,
                          "rural" = "Rural measles outbreak - deploy team",
                          "urban" = "Urban measles outbreak - notify ministry"
                   )
               }
           },
           
           "cholera" = {
               if(cases < 2){
                   "No alert"
               } else {
                   switch(location,
                          "rural" = "Rural cholera outbreak - deploy WASH team",
                          "urban" = "Urban cholera outbreak - notify ministry"
                   )
               }
           },
           
           "ebola" = {
               if(cases == 0){
                   "No alert"
               } else {
                   "IMMEDIATE EBOLA ALERT - notify WHO"
               }
           },
           
           "Disease not in database"
    )
}


surveillance_alert("measles", 6, "urban")
