# eMethod2. Multiple imputation (R code)
# NOTE:
# All binary variables in the original dataset were coded as 0/1,	
# where 0 indicates absence and 1 indicates presence of the condition or event.
# No recoding from character or alternative numeric formats was required.
# Binary variables were explicitly converted to factors with levels (0, 1)
# to ensure correct interpretation as categorical predictors in the imputation models.
library(mice)
set.seed(123)

data1 <- read.csv("Mice_data.csv", stringsAsFactors = FALSE)
names(data1) <- make.names(names(data1), unique = TRUE)

keep_exist <- function(x, nms) intersect(x, nms)

## ---- 1) Define variable types (only for variables you will use) ----
# outcome (0/1)
stopifnot("antitussive" %in% names(data1))
data1$antitussive <- as.integer(data1$antitussive)
stopifnot(all(is.na(data1$antitussive) | data1$antitussive %in% c(0L, 1L)))

# exposure (factor)
stopifnot("exposure" %in% names(data1))
data1$exposure <- factor(data1$exposure)

# binary covariates (already 0/1) -> factor(0,1)
bin_vars <- c("Sex", "History_of_hypertension", "Diabetes_mellitus", "Chronic_obstructive_pulmonary_disease", "Coronary_heart_disease", "History_of_tumor", "History_of_tuberculosis", "Preoperative_cough", "Preoperative_pulmonary_infection", "Preoperative_hypoxemia", "History_of_smoking", "Alcohol_consumption", "Hypercapnia", "Mediastinal_swing", "Difficult_intubation", "Conversion_to_endotracheal_intubation", "Patient_controlled_intravenous_analgesia", "Intraoperative_transfusion", "Chest_tube_placement", "Urinary_catheterization", "Right_lung_surgery", "Upper_lobe_resection", "Hilar_lymph_node_dissection", "Subcarinal_lymph_node_dissection", "Paratracheal_lymph_node_dissection", "Acute_cough", "Chronic_cough", "Hoarseness", "Sore_throat", "Postoperative_pulmonary_complications", "Respiratory_infection", "Respiratory_failure", "Pleural_effusion", "Atelectasis", "Pneumothorax", "Bronchospasm", "Aspiration_pneumonia", "Antibiotic_use_beyond_48_hours", "Bronchoscopy", "Constipation", "Venous_thromboembolism", "Sedative_use", "Postoperative_rescue_opioid_use", "New_onset_arrhythmia", "Readmission_due_to_complications", "Acute_kidney_injury", "Acute_liver_injury" )
bin_vars <- keep_exist(bin_vars, names(data1))
for (v in bin_vars) {
  data1[[v]] <- as.integer(data1[[v]])
  stopifnot(all(is.na(data1[[v]]) | data1[[v]] %in% c(0L, 1L)))
  data1[[v]] <- factor(data1[[v]], levels = c(0, 1))
}

# ordinal covariates
ord_vars <- c("ASA","NYHA") 
ord_vars <- keep_exist(ord_vars, names(data1))
for (v in ord_vars) data1[[v]] <- ordered(data1[[v]])

# nominal categorical covariates
nom_vars <- c("Type_of_lung_resection")
nom_vars <- keep_exist(nom_vars, names(data1))
for (v in nom_vars) data1[[v]] <- factor(data1[[v]])

# continuous covariates to be imputed with pmm
vars_cont <- c("MVV","FEV1","Predicted_FEV1","FEV1_FVC","DLCO",
               "PaO2_FiO2","PaCO2","Lactate","SpO2","WBC","RBC",
               "HCT","Hb","PLT","ANC","Neutrophil_percentage")
stopifnot(all(vars_cont %in% names(data1)))
for (v in vars_cont) data1[[v]] <- suppressWarnings(as.numeric(data1[[v]]))

## ---- 2) Specify imputation models ----
pred <- make.predictorMatrix(data1)

# exclude ID from prediction (if present)
if ("NO" %in% names(data1)) {
  pred[, "NO"] <- 0
  pred["NO", ] <- 0
}

# impute only continuous covariates
meth <- make.method(data1)
meth[] <- ""                # default: no imputation
meth[vars_cont] <- "pmm"    # only continuous variables
meth["antitussive"] <- ""  # outcome not imputed
meth["exposure"] <- ""     # exposure not imputed
if ("NO" %in% names(data1)) meth["NO"] <- ""

imp <- mice(
  data1,
  m = 5,
  method = meth,
  predictorMatrix = pred,
  maxit = 20,
  seed = 123,
  printFlag = FALSE
)
## ---- 3) Analysis model + pooling ----
fit <- with(imp, glm(
  antitussive ~ exposure + MVV + FEV1 + Predicted_FEV1 + FEV1_FVC + DLCO +
    PaO2_FiO2 + PaCO2 + Lactate + SpO2 + WBC + RBC + HCT + Hb +
    PLT + ANC + Neutrophil_percentage,
  family = binomial()
))

res <- summary(pool(fit), conf.int = TRUE, exponentiate = TRUE)
write.csv(res, "MI_logistic_results.csv", row.names = FALSE)
