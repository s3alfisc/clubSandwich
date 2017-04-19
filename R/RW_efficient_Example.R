
library(clubSandwich)
library(robumeta)

data(dropoutPrevention)

##Build a "NULL" model with robu() function
robu_result_null <- robu(LOR1 ~ study_design + service_hrs+ duration, 
                    data = dropoutPrevention, studynum = studyID, var.eff.size = varLOR, 
                    modelweights = "HIER",small = T)

user_weights <- robu_result_null$data.full$r.weights

##Build a "MAIN" model with robu() function. (For comparison)
robu_result_main <- robu(LOR1 ~ study_design + service_hrs + duration+average_age, 
                    data = dropoutPrevention, studynum = studyID, var.eff.size = varLOR,
                    userweights = robu_result_null$data.full$r.weights,small = F)

robu_result_main$reg_table[6,]

##Constraints is a formula.
constraints<- as.formula("~average_age-1")
obj = robu_result_null
constraints= constraints
cluster = NULL
residual_adjustment = "CR2"
test_adjustment = "CR2"
random_seed = 2252
bootstraps = 2000
##Wild bootstrap result. (NULL HYPOTHESIS: average_age has no effect)
result<- Wild_bootstrap(obj = robu_result_null,
               constraints= constraints,
               cluster = NULL,
               residual_adjustment = "CR2", 
               test_adjustment = "CR2", 
               random_seed = 2252)

result$res_vec
