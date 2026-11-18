
# Load Libraries ----------------------------------------------------------

pacman::p_load(
    tidyverse,
    rio,
    modelsummary,
    broom,
    gtsummary,
    sjPlot
)


# Load Data ---------------------------------------------------------------

data <- import("https://evalf20.classes.andrewheiss.com/data/injury.csv")
data_raw <- data |> 
    filter(
       mi == 0
    )


# EDA ---------------------------------------------------------------------



# DiD by hand -------------------------------------------------------------

diff <- data_raw |> 
    summarise(
        avg_dur = mean( ldurat, na.rm = T),
        .by = c(afchnge, highearn)
    )
    

bef_cont <-  diff |> 
    filter(
        afchnge == 0,
        highearn == 0
    ) |> 
    pull(avg_dur)


bef_trt <-  diff |> 
    filter(
        afchnge == 0,
        highearn == 1
    ) |> 
    pull(avg_dur)


aft_cont <-  diff |> 
    filter(
        afchnge == 1,
        highearn == 0
    ) |> 
    pull(avg_dur)


aft_trt <-  diff |> 
    filter(
        afchnge == 1,
        highearn == 1
    ) |> 
    pull(avg_dur)



diff_matrix <- matrix(c(
    
    bef_cont, aft_cont, 
    bef_trt, aft_trt 
), byrow = T, nrow = 2
)

colnames(diff_matrix) <- c("Pre", "Post")
row.names(diff_matrix) <- c("Control", "Treatment")

diff_matrix

diff_df_matrix <- as.data.frame.matrix(diff_matrix) |> 
    mutate(
        diff = Post - Pre
    )

did <- (diff_matrix[4] -  diff_matrix[2]) - (diff_matrix[3] -  diff_matrix[1])

diff_in_diffs <- diff_df_matrix[2,3] - diff_df_matrix[1,3]

# DiD Regression ----------------------------------------------------------


data_raw |> 
    group_by(afchnge, highearn) |> 
    count()


model <- lm(
    ldurat ~ highearn +  afchnge + (highearn * afchnge),
    data = data_raw
)

tbl_regression(
    model,
    intercept = T,
    label = list(
        highearn ~"High Earners",
        afchnge ~ "After Change"
    )
) |> 
    bold_p() |> 
    bold_labels() |> 
    modify_caption("**Difference-in-Differences Model: Kentucky Workers' Compensation (1980)**")



tab_model(
    model,
    show.p = F
    )
