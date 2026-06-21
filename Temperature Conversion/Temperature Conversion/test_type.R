test_type <- function(sample){
    
    switch (sample,
    "one" = "one t.test",
    "two" = "two sample t.test",
    "three" = "anova",
    "could be a non-parametric test"
    )
}


test_type("regression")


clinical_desc <- function(disease, severity){
    
    switch(disease,
           
           "malaria" = switch( severity,
                                "mild" = "Prescribe oral artemether",
                               "severe" = "Admit and give IV artesunate"
           ),
           
           
           "cholera" = switch( severity,
                               "mild" = "Oral rehydration therapy",
                               "severe" = "Admit and give IV fluids"
           ),
           
           "tuberculosis" = switch( severity,
                                    "mild" = "Start standard DOTS therapy",
                                    "severe" = "Admit and start DOTS therapy"
               
           ),
           
           "condition not in the dataset"
        
        
        
    )
    
}

clinical_desc("tubercu", "very very")
