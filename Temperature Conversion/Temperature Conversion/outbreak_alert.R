outbreak_alert <- function(disease, cases){
    
    
    
    switch(disease,
           
           "measles" = {
               if (cases < 5){
                   "Watch"
               } else if(cases <= 10){
                   "Watch"
               } else{
                   "Outbreak alert!"
               }
           }
        
    )
    
    
    
    switch(disease,
           
           "cholera" = {
               if (cases < 2){
                   "Watch"
               } else if(cases > 2 & cases <= 5 ){
                   "Watch"
               } else{
                   "Outbreak alert!"
               }
           }
           
    )
    
    
    switch(disease,
           
           "ebola" = {
               if (cases  == 0){
                   "No alert"
               } else if(cases >= 1){
                   "Immediate outbreak alert!"
               }
           }
           
    )
    
    "Condition not a prority"
    
    
    
    
}

outbreak_alert("dysentary", 0)

