risk_classifier <- function(age, bmi, smoker){
    
    if(age > 60 & smoker == "yes" & bmi > 30){
        "Critical risk"
        
    } else if(age > 60 & smoker == "yes"){
        "High risk"
        
    } else if(age > 60 & bmi > 30){
        "High risk"
        
    } else if(bmi > 30 & smoker == "yes"){
        "High risk"
        
    } else if(age > 60){
        "Moderate risk"
        
    } else if(bmi > 30){
        "Moderate risk"
        
    } else if(smoker == "yes"){
        "Moderate risk"
        
    } else {
        "Low risk"
    }
}


risk_classifier(age = 60, smoker = "20", bmi = 20)
