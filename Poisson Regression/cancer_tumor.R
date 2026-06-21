df <- data.frame(
    time = c(2, 3, 6, 8, 9, 10, 11, 13, 14, 16, 21, 22, 24, 26, 27,
             7, 13, 15, 18, 23, 20, 24,
             1, 5, 17, 18, 25,
             18, 25,
             4, 19),
    
    x = c(rep(0, 22), rep(1, 9)),
    
    n = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
          2,2,2,2,2,3,4,
          1,1,1,1,1,
          2,2,
          3,4)
)


poi_muv <- glm(n ~ x+offset(log(time)),
               family = poisson,
               data = df)

summary(poi_muv)
check_overdispersion(poi_muv)

tbl_regression(
    poi_muv,
    exponentiate = T,
    label = list(
        x = md("Size of the Tumor")
    ),
    pvalue_fun = label_style_pvalue(digits = 3)
) |> 
    bold_labels() |> 
    as_gt() |> 
    tab_header(
        md("**Predictors of Tumors**")
    )
