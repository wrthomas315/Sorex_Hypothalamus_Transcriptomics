# Creating the data frame
data <- data.frame(
  Value1 = c(0.036663, 0.02633, 0.029997, 0.03333, 0.031663, 0.026997, 
             0.011997, 0.006997, 0.004997, 0.009663, 0.012997, 0.008663),
  Group = c('ControlOverexpression', 'ControlOverexpression', 'ControlOverexpression',
            'ControlScrambled', 'ControlScrambled', 'ControlScrambled',
            'HeatOverexpression', 'HeatOverexpression', 'HeatOverexpression',
            'HeatScrambled', 'HeatScrambled', 'HeatScrambled')
)

# Extracting the values for each group
control_overexpression <- data %>% filter(Group == 'ControlOverexpression') %>% pull(Value1)
control_scrambled <- data %>% filter(Group == 'ControlScrambled') %>% pull(Value1)
heat_overexpression <- data %>% filter(Group == 'HeatOverexpression') %>% pull(Value1)
heat_scrambled <- data %>% filter(Group == 'HeatScrambled') %>% pull(Value1)

# Performing t-tests
ttest_results <- list(
  'ControlOverexpression vs ControlScrambled' = t.test(control_overexpression, control_scrambled),
  'ControlOverexpression vs HeatOverexpression' = t.test(control_overexpression, heat_overexpression),
  'ControlOverexpression vs HeatScrambled' = t.test(control_overexpression, heat_scrambled),
  'ControlScrambled vs HeatOverexpression' = t.test(control_scrambled, heat_overexpression),
  'ControlScrambled vs HeatScrambled' = t.test(control_scrambled, heat_scrambled),
  'HeatOverexpression vs HeatScrambled' = t.test(heat_overexpression, heat_scrambled)
)

# Print the t-test results
print(ttest_results)

# Combine all control groups and all heat groups
all_control <- c(control_overexpression, control_scrambled)
all_heat <- c(heat_overexpression, heat_scrambled)

# Perform t-test comparing all control values with all heat values
overall_ttest <- t.test(all_control, all_heat)

# Print the overall t-test result
print(overall_ttest)
