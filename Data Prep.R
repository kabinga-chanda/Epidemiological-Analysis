# Load the libraries
pacman::p_load(epitools,
               tidyverse,
               rio, 
               janitor, 
               gt, 
               gtsummary,
               naniar,
               epikit,
               summarytools,
               broom,
               epiR,
               brms,
               cmdstanr,
               BH
               )
# Load the data
data(oswego)

data <- oswego


# Re-code the variables
data_clean <- data |>
    mutate(
        across(8:21,
               ~ if_else(. == "Y", 1,
                         if_else(. == "N", 0, NA_real_)))
    )

### Recode the ill var
data_clean <- data |>
    mutate(
        ill = if_else(ill == "Y", 1,
                      if_else(ill == "N", 0, NA_real_)),
        sex = if_else(sex == "F", "Female", "Male"),
        onset.datetime = mdy_hm(paste(onset.date, "2024", onset.time))
    )

# mutateage groups
data_clean <- data_clean |> 
    mutate(age_grp = age_categories(age,
                                    by = 15,
                                    lower = 3,
                                    upper = 77))

    
# Create the Epidemic Curve
data_clean |>
    mutate(ill = if_else(ill == "1", "Ill", "Not Ill")) |> 
ggplot(aes(x = onset.datetime)) +
    geom_histogram(aes(fill = ill), bins = 30, color = "black") +
    scale_fill_manual(
        values = c("Not Ill" = "lightblue", "Ill" = "lightpink"),
        labels = c("Not Ill", "Ill"),
        name = "Illness Status"
    ) +
    scale_x_datetime(
        date_breaks = "3 hours",
        date_labels = "%m/%d %I:%M %p",
        expand = expansion(mult = c(0.1, 0.1))
    ) +
    labs(
        title = "Epidemic Curve of Foodborne Illness Outbreak",
        x = "Date and Time of Onset",
        y = "Number of Cases"
    ) +
    theme(
        panel.background = element_blank(),
        axis.line.x = element_line(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.y = element_line(linewidth = 0.3, color = "grey")
    )

# Characteristics and Demographics
data_clean |>
    mutate(ill = if_else(ill == "1", "Ill", "Not Ill")) |> 
    select(-c(id,
              onset.date,
              onset.time,
              age_grp,meal.time, 
              onset.datetime,  
              8:21)) |> 
    tbl_summary(by = ill,
                label = list(
                    age = "Age (In Years)",
                    sex = "Sex of the Participant")
    ) |> 
    bold_labels() |> 
    add_overall(
        col_label = "**Overall**  \nN = {style_number(N)}",
    ) |> 
    add_p() |> 
    as_gt() |>
    tab_header(
        title = md("**Outbreak of Gastrointestinal Illness in Oswego County, 1940**")
    ) |>
    tab_source_note(
        source_note = md("**Source:** Oswego An Outbreak of Gastrointestinal Illness Following a Church Supper<br>
                         (updated 2003): *S. aureus* outbreak among church picnic attendees, 1940; the classic,<br>
                         straightforward outbreak investigation in a defined population.<br>
                         Training modules available at<br>
                         *https://www.cdc.gov/eis/casestudies/xoswego.401-303.student.pdf*")
    ) |>
    opt_align_table_header(align = "left") |>
    tab_options(
        table.border.top.style = "solid",
        table.border.bottom.style = "solid",
        column_labels.border.bottom.style = "solid",
        table_body.hlines.style = "none",
        table_body.vlines.style = "none"
    ) 
    
# Logisitic Regression


food_vars <- c("baked.ham", "spinach", "mashed.potato", "cabbage.salad", "jello",
               "rolls", "brown.bread", "milk", "coffee", "water",
               "cakes", "vanilla.ice.cream", "chocolate.ice.cream", "fruit.salad")


tbl_uvregression(
    data = data_clean,
    method = glm,
    y = ill,
    method.args = list(family = binomial),
    include = c(all_of(food_vars), age, sex),
    exponentiate = TRUE,
    show_single_row = -sex,
    label = list(
        baked.ham ~ "Baked Ham",
        spinach ~ "Spinach",
        mashed.potato ~ "Mashed Potato",
        cabbage.salad ~ "Cabbage Salad",
        jello ~ "Jello",
        rolls ~ "Rolls",
        brown.bread ~ "Brown Bread",
        milk ~ "Milk",
        coffee ~ "Coffee",
        water ~ "Water",
        cakes ~ "Cakes",
        vanilla.ice.cream ~ "Vanilla Ice Cream",
        chocolate.ice.cream ~ "Chocolate Ice Cream",
        fruit.salad ~ "Fruit Salad",
        age ~ "Age (years)",
        sex ~ "Sex"
    ))|> 
    as_gt() |> 
    tab_options(
        table.border.top.style = "solid",
        table.border.bottom.style = "solid",
        column_labels.border.bottom.style = "solid",
        table_body.hlines.style = "none",
        table_body.vlines.style = "none"
    ) 


#Multivariate Logistic Regression
model_all <- glm(
    reformulate(c(food_vars, "age", "sex"), response = "ill"),
    data = data_clean,
    family = binomial
)

tbl_regression(
    model_all,
    show_single_row = -sex,
    exponentiate = TRUE,
    label = list(
        baked.ham ~ "Baked Ham",
        spinach ~ "Spinach",
        mashed.potato ~ "Mashed Potato",
        cabbage.salad ~ "Cabbage Salad",
        jello ~ "Jello",
        rolls ~ "Rolls",
        brown.bread ~ "Brown Bread",
        milk ~ "Milk",
        coffee ~ "Coffee",
        water ~ "Water",
        cakes ~ "Cakes",
        vanilla.ice.cream ~ "Vanilla Ice Cream",
        chocolate.ice.cream ~ "Chocolate Ice Cream",
        fruit.salad ~ "Fruit Salad",
        age ~ "Age (years)",
        sex ~ "Sex"
    )
) |> 
    bold_labels() |> 
    bold_p() |> 
    as_gt() |> 
    tab_header(title = md("**Multivariate Logistic Regression**")) |> 
    opt_align_table_header(align = "left") |> 
    tab_options(
        table.border.top.style = "solid",
        table.border.bottom.style = "solid",
        column_labels.border.bottom.style = "solid",
        table_body.hlines.style = "none",
        table_body.vlines.style = "none"
    ) |> 
    cols_label(
        estimate = md("**AOR**"),
        conf.low = md("**95% CI**")
    )

# Attack Rates

att_rat_data <-  data_clean |>
    rename(
        "Baked Ham" = baked.ham,
        "Spinach" = spinach,
        "Mashed Potato" = mashed.potato,
        "Cabbage Salad" = cabbage.salad,
        "Jello" = jello,
        "Rolls" = rolls,
        "Brown Bread" = brown.bread,
        "Milk" = milk,
        "Coffee" = coffee,
        "Water" = water,
        "Cakes" = cakes,
        "Vanilla Ice Cream" = vanilla.ice.cream,
        "Chocolate Ice Cream" = chocolate.ice.cream,
        "Fruit Salad" = fruit.salad
    )


att_rat_data_long <- att_rat_data |>
    pivot_longer(
        cols = 8:21,
        names_to = "food_item",
        values_to = "exposed"
    ) |> select(food_item, exposed, ill)

att_rat_data_long |> 
    group_by(food_item) |> 
    summarise("Ill" = sum(ill == "1"),
              "Not Ill"  = sum(ill == "0"),
              "Attack Rate" = Ill/nrow(data_clean))



att_rat_summary <- att_rat_data_long |>
    group_by(food_item) |>
    summarise(
        Ill_exposed = sum(exposed == "Y" & ill == "1", na.rm = TRUE),
        NotIll_exposed = sum(exposed == "Y" & ill == "0", na.rm = TRUE),
        Ill_not_exposed = sum(exposed == "N" & ill == "1", na.rm = TRUE),
        NotIll_not_exposed = sum(exposed == "N" & ill == "0", na.rm = TRUE),
        .groups = "drop"
    )


att_rat_summary <- att_rat_summary |>
    mutate(
        Total_exposed = Ill_exposed + NotIll_exposed,
        Total_not_exposed = Ill_not_exposed + NotIll_not_exposed,
        AR_exposed = round(Ill_exposed / Total_exposed * 100,2),
        AR_not_exposed = round(Ill_not_exposed / Total_not_exposed * 100,2),
        Attack_Rate_Ratio = round(AR_exposed / AR_not_exposed, 2)
    )


att_rat_summary <- att_rat_summary |> 
    arrange(desc(Attack_Rate_Ratio))


att_rat_summary |>
    gt() |>
    tab_header(
        title = md("**Food-specific Attack Rate Table**"),
        subtitle = "Oswego County Outbreak Investigation, 1940"
    ) |>
    # Add spanners
    tab_spanner(
        label = md("**Ate Food**"),
        columns = c(Ill_exposed, NotIll_exposed, Total_exposed, AR_exposed)
    ) |>
    tab_spanner(
        label = md("**Did Not Eat Food**"),
        columns = c(Ill_not_exposed, NotIll_not_exposed, Total_not_exposed, AR_not_exposed)
    ) |>
    cols_label(
        food_item = md("**Food Item**"),
        Ill_exposed = md("**Ill**"),
        NotIll_exposed = md("**Not Ill**"),
        Total_exposed = md("**Total**"),
        AR_exposed = md("**Attack<br> Rate (%)**"),
        Ill_not_exposed = md("**Ill**"),
        NotIll_not_exposed = md("**Not Ill**"),
        Total_not_exposed = md("**Total**"),
        AR_not_exposed = md("**Attack<br> Rate (%)**"),
        Attack_Rate_Ratio = md("**Attack Rate<br> Ratio**")
    ) |>
    fmt_number(
        columns = c(AR_exposed, AR_not_exposed, Attack_Rate_Ratio),
        decimals = 1
    ) |>
    tab_source_note(
        source_note = "AR = Attack Rate; Exposure = Ate the food (Y/N)"
    ) |> 
    opt_align_table_header(align = "left") |> 
    tab_options(
        table.border.top.style = "solid",
        table.border.bottom.style = "solid",
        column_labels.border.bottom.style = "solid",
        table_body.hlines.style = "none",
        table_body.vlines.style = "none"
    )



### Stratied by Sex

# Step 1: Pivot foods to long format
att_rat_data_long <- data_clean |>
    pivot_longer(
        cols = baked.ham:fruit.salad,   # all food columns
        names_to = "food_item",
        values_to = "exposed"
    )




# Step 2: Calculate stratified attack rates by sex
att_rat_sex <- att_rat_data_long |>
    group_by(food_item, sex) |>
    summarise(
        Ill_exposed = sum(exposed == "Y" & ill == 1, na.rm = TRUE),
        NotIll_exposed = sum(exposed == "Y" & ill == 0, na.rm = TRUE),
        Ill_not_exposed = sum(exposed == "N" & ill == 1, na.rm = TRUE),
        NotIll_not_exposed = sum(exposed == "N" & ill == 0, na.rm = TRUE),
        Total_exposed = Ill_exposed + NotIll_exposed,
        Total_not_exposed = Ill_not_exposed + NotIll_not_exposed,
        AR_exposed = round(Ill_exposed / Total_exposed * 100, 1),
        AR_not_exposed = round(Ill_not_exposed / Total_not_exposed * 100, 1),
        Attack_Rate_Ratio = round(AR_exposed / AR_not_exposed, 2),
        .groups = "drop"
    ) |>
    arrange(food_item, sex)


# Create a named vector for nicer food names
food_names <- c(
    "baked.ham" = "Baked Ham",
    "brown.bread" = "Brown Bread",
    "cabbage.salad" = "Cabbage Salad",
    "cakes" = "Cakes",
    "chocolate.ice.cream" = "Chocolate Ice Cream",
    "coffee" = "Coffee",
    "fruit.salad" = "Fruit Salad",
    "jello" = "Jello",
    "mashed.potato" = "Mashed Potato",
    "milk" = "Milk",
    "rolls" = "Rolls",
    "spinach" = "Spinach",
    "vanilla.ice.cream" = "Vanilla Ice Cream",
    "water" = "Water"
)

# Apply nicer names to your table
att_rat_sex <- att_rat_sex |>
    mutate(food_item = food_names[food_item])


# Step 3: Create a publication-ready table with gt
att_rat_sex |> 
    gt() |>
    tab_header(
        title = md("**Stratified Food-specific Attack Rates by Sex**"),
        subtitle = "Oswego-style Outbreak Analysis"
    ) |>
    tab_spanner(
        label = md("**Ate Food**"),
        columns = c(Ill_exposed, NotIll_exposed, Total_exposed, AR_exposed)
    ) |>
    tab_spanner(
        label = md("**Did Not Eat Food**"),
        columns = c(Ill_not_exposed, NotIll_not_exposed, Total_not_exposed, AR_not_exposed)
    ) |>
    cols_label(
        food_item = md("**Food Item**"),
        sex = "Sex",
        Ill_exposed = md("**Ill**"),
        NotIll_exposed =md("**Not Ill**"),
        Total_exposed = md("**Total**"),
        AR_exposed = md("**Attack<br> Rate (%)**"),
        Ill_not_exposed = md("**Ill**"),
        NotIll_not_exposed = md("**Not Ill**"),
        Total_not_exposed = md("**Total**"),
        AR_not_exposed = md("**Attack <br>Rate (%)**"),
        Attack_Rate_Ratio = md("**Attack <br>Rate Ratio**")
    ) |>
    fmt_number(
        columns = c(AR_exposed, AR_not_exposed, Attack_Rate_Ratio),
        decimals = 1
    ) |>
    tab_source_note(
        source_note = "AR = Attack Rate; Exposure = Ate the food (Y/N)"
    ) |> 
    tab_options(
        table.border.top.style = "solid",
        table.border.bottom.style = "solid",
        column_labels.border.bottom.style = "solid",
        table_body.hlines.style = "none",
        table_body.vlines.style = "none"
    )



### Attack Rate Plot

att_rat_summary |>
    ggplot(aes(x = reorder(food_item, Attack_Rate_Ratio), y = Attack_Rate_Ratio)) +
    geom_col(fill = "thistle",
             color = "black") +
    coord_flip() +
    labs(y = "Attack Rate Ratio", x = "Food Item")+
    theme(panel.background = element_blank(),
          axis.line.x = element_line(),
          panel.grid.major.x = element_line(color = "grey",
                                            size = 0.04)
          )



##

brm_model <- brm(ill ~ baked.ham + spinach + mashed.potato,
                 family = bernoulli(), data = data_clean)
summary(brm_model)



# Install dependencies one at a time
install.packages("bayesplot")
install.packages("bridgesampling")
install.packages("nleqslv")

# Then install brms
install.packages("brms")







---
    title: "Oswego Foodborne Outbreak Investigation, 1940"
subtitle: "Analysis of Gastrointestinal Illness Following a Church Supper"
author: "Kainga Chanda"
date: today
format:
    html:
    toc: true
toc-depth: 3
toc-location: left
number-sections: true
theme: cosmo
code-fold: true
code-tools: true
embed-resources: true
fig-width: 8
fig-height: 6
smooth-scroll: true
editor: visual
execute:
    warning: false
message: false
echo: false
---
