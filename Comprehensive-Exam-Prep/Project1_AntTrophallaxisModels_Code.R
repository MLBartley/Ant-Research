####
####
#### Project One Final* Code
#### February 2017
####


# *Wishful thinking; who am I kidding? 


###Outline
    # Ant data visualization - 4 hour High/Low density
            #histograms of events per ant
            #scatterplots of each density (and each location in low density)
            #scatter plots of each density with entrance times     
    # Ant data simple model - want to show motivation - doesn't work!
    # Simulated data - visualization 
            # scatter plots of each simulated density?   
    # Simulated data penalized model - show that it can be effective
    # Simulated data penalized model with covariates - include the biology
  

#Call in Trophallaxis Data
    # Currently the .cvs files also load a bunch of columns that are empty 
    # In low density, chamber 4 is entrance and chamber 1 has queen

col1_high4 <- read.csv("./Data/Colony1_trophallaxis_high_density_4hr.csv")
col1_low4 <- read.csv("./Data/Colony1_trophallaxis_low_density_4hr.csv")

col2_high4 <- read.csv("./Data/Colony2_trophallaxis_high_density_4hr.csv")
col2_low4 <- read.csv("./Data/Colony2_trophallaxis_low_density_4hr.csv")

col3.high4 <- read.csv("./Data/Colony3_trophallaxis_high_density_4hr.csv")
col3_low4 <- read.csv("./Data/Colony3_trophallaxis_low_density_4hr.csv")


#Call in Foraging Data

inout_high4 <- read.csv("./Data/Colony1_in&out_high_density_4hr.csv")
inout_low4 <- read.csv("./Data/Colony1_in&out_low_density_4hr.csv")



#Visualize the Trophallaxis Data - NEED TO UPDATE SUMVIS FUNCTION
    #Want these to save to .pdf (in new folder?)

sumvis_troph(data = col1_low4, entrance = FALSE, hours = 4, density = "low")

sumvis_troph(data = col2_low4, entrance = F, hours = 4, density = "low")

sumvis_troph(data = col3_low4, entrance = F, hours = 4, density = "low")


#Prep Trophallaxis Data - note decided to keep prep.torph.data function 


col1_high4_5 = prep_troph_data(col1_high4, hours = 4, delta_t =  5)
col1_low4_5 = prep_troph_data(col1_low4, hours = 4, delta_t = 5)

col1.high4.30 = prep_troph_data(col1.high4, 30)
col1.low4.30 = prep_troph_data(col1.low4, 30)


col2.high4.5 = prep_troph_data(col2.high4, 5)
col2_low4_5 = prep_troph_data(col2_low4, hours = 4, delta_t =  5)

col2.high4.30 = prep_troph_data(col2.high4, 30)
col2.low4.30 = prep_troph_data(col2.low4, 30)


col3.high4.5 = prep_troph_data(col3.high4, 5)
col3_low4_5 = prep_troph_data(col3_low4, hours = 4, 5)

col3.high4.30 = prep_troph_data(col3.high4, 30)
col3.low4.30 = prep_troph_data(col3.low4, 30)


#Prep In & Out Data

col1_hi4_inout_1 <- prep_inout_data(data = inout_high4, delta_t = 1, hours = 4)
col1_lo4_inout_1 <- prep_inout_data(data = inout_low4, delta_t = 1, hours = 4)

col1_hi4_inout_5 <- prep_inout_data(data = inout_high4, delta_t = 5, hours = 4)
col1_lo4_inout_5 <- prep_inout_data(data = inout_low4, delta_t = 5, hours = 4)




# Variables needed for all 

path <- "./Comprehensive-Exam-Prep/output_images/"
states <- 2
n_mcmc <- 5000
hours <- 4
X <- sample(x = c(1, 2), size = hours*60*60, replace = T)
lambda <- c(.01, .08)
registerDoMC(cores=5)


#Ant Data - Simple Model 

# theta <- matrix(data = c(5000, 1, 1, 5000), nrow = 2, ncol = 2, byrow = T) 


P <- matrix(c(.997, .003, .003, .997), nrow = 2, byrow = T)

param_start <- list(X = X, lambda = lambda, P = P)

simple_col1hi4bin1 <- DT_mcmc_troph(starts_data = col1_high4_5$starts_persec, 
                                    ant_file = col1_high4_5$data, title = "Test", 
                                    a = .005, b = .001, c = .005, d = .001, 
                                    theta = theta, states = states, 
                                    n_mcmc = n_mcmc, delta_t = 1, hours = hours,
                                    param_start = param_start, fig_save = TRUE, 
                                    fig_path = path, fig_name = "simp_col1hi4bin1") 

simple_col1lo4tbin1 <- DT_mcmc_troph(starts_data = col1_low4_5$starts_persec, 
                                    ant_file = col1_low4_5$data, title = "Test", 
                                    a = .005, b = .001, c = .005, d = .001, 
                                    theta = theta, states = states, 
                                    n_mcmc = n_mcmc, delta_t = 1, hours = hours,
                                    param_start = param_start, fig_save = TRUE, 
                                    fig_path = path, fig_name = "simp_col1lo4tbin1") 

simple_col1lo4qbin1 <- DT_mcmc_troph(starts_data = col1_low4_5$queen_starts_persec, 
                                    ant_file = col1_high4_5$data, title = "Test", 
                                    a = .005, b = .001, c = .005, d = .001, 
                                    theta = theta, states = states, 
                                    n_mcmc = n_mcmc, delta_t = 1, hours = hours,
                                    param_start = param_start, fig_save = TRUE, 
                                    fig_path = path, fig_name = "simp_col1lo4qbin1") 

simple_col1lo4ebin1 <- DT_mcmc_troph(starts_data = col1_low4_5$entrance_start_persec, 
                                    ant_file = col1_low4_5$data, title = "Test", 
                                    a = .005, b = .001, c = .005, d = .001, 
                                    theta = theta, states = states, 
                                    n_mcmc = n_mcmc, delta_t = 1, hours = hours,
                                    param_start = param_start, fig_save = TRUE, 
                                    fig_path = path, fig_name = "simp_col1lo4ebin1") 


simple_col2lo4qbin1 <- foreach (i = seq(15000, 25000, by =  1000) ,
                                .errorhandling="remove") %dopar% 
                                    DT_mcmc_troph(starts_data = col2_low4_5$queen_starts_persec, 
                                                 ant_file = col2_low4_5$data, chamber = "queen", 
                                                 title = "Test", 
                                                 a = .005, b = .001, c = .005, d = .001, 
                                                 theta = matrix(c(i, 1, 1, i), 2, 2), states = states, 
                                                 n_mcmc = n_mcmc, delta_t = 1, hours = hours,
                                                 param_start = param_start, fig_save = TRUE, 
                                                 fig_path = path, fig_name = "simp_col2lo4qbin1") 

sumtable_model(results = simple_col2lo4qbin1, compare = seq(5000, 15000, by =  1000), 
               file_path = "./Comprehensive-Exam-Prep/output_tables/", 
               file_name = "sim_col2lo4qbin1", model = "simple")


simple_col3lo4ebin1 <- foreach (i = seq(5000, 5000, by =  1000) ,
                                .errorhandling="remove") %dopar% DT_mcmc_troph(starts_data = col3_low4_5$entrance_start_persec, 
                                                                               ant_file = col3_low4_5$data, title = "Test", 
                                                                               a = .005, b = .001, c = .005, d = .001, 
                                                                               theta = matrix(c(i, 1, 1, i), 2, 2), states = states, 
                                                                               n_mcmc = n_mcmc, delta_t = 1, hours = hours,
                                                                               param_start = param_start, fig_save = TRUE, 
                                                                               fig_path = path, fig_name = "simp_col3lo4ebin1") 

sumtable_model(results = simple_col3lo4ebin1, compare = seq(5000, 5000, by =  1000), 
               file_path = "./Comprehensive-Exam-Prep/output_tables/", 
               file_name = "sim_col3lo4ebin1", model = "simple")




#Ant Data - Penalized Model 
penalty <- exp(seq(-25, 1, by =  3)) 
tau <- matrix( c(.00001, 0, 
    0, .00001), nrow = 2, ncol = 2)
gamma <- c(.005, .005)
start <- list(X = X, lambda = lambda, gamma = gamma)
delta_t <- 1


# penalize_col1hi4bin1 <- lapply(penalty, FUN = DT_pen_mcmc,
#                               starts_data = col1_high4_5$starts_persec,
#                               states = states, ant_file = col1_high4_5$data,
#                               hours = hours,
#                               a = .005, b = .001, c = .005, d = .001,
#                               tau = tau, tau.pen = 0, n_mcmc = n_mcmc,
#                               delta_t = delta_t, start = start, fig_save = TRUE,
#                               fig_path = path,
#                               fig_name = "pen_col1hi4bin1_")

penalize_col1hi4bin1 <- foreach (i = exp(seq(-25, -5, by =  2)) ,
                                 .errorhandling="remove") %dopar% 
                                      DT_pen_mcmc(penalty = i, starts_data = col1_high4_5$starts_persec, 
                                                  states = states, ant_file = col1_high4_5$data,
                                                  hours = hours, 
                                                  a = .005, b = .001, c = .005, d = .001,
                                                  tau = tau, tau.pen = 0, n_mcmc = n_mcmc, 
                                                  delta_t = delta_t, start = start, fig_save = TRUE,
                                                  fig_path = path, 
                                                  fig_name = "pen_col1hi4bin1_")

sumtable_model(results = penalize_col1hi4bin1, compare = penalty, 
               file_path = "./Comprehensive-Exam-Prep/output_tables/", 
               file_name = "pen_col1hi4bin1", model = "penalized")


penalize_col1lo4tbin1 <- foreach (i = exp(seq(-25, -5, by =  1)) ,
                                 .errorhandling="remove") %dopar% 
                          DT_pen_mcmc(penalty = i, starts_data = col1_low4_5$starts_persec, 
                                      states = states, ant_file = col1_low4_5$data,
                                      hours = hours, 
                                      a = .005, b = .001, c = .005, d = .001,
                                      tau = tau, tau.pen = 0, n_mcmc = n_mcmc, 
                                      delta_t = delta_t, start = start, fig_save = TRUE,
                                      fig_path = path, 
                                      fig_name = "pen_col1lo4tbin1_")

sumtable_model(results = penalize_col1lo4tbin1, compare = penalty, 
               file_path = "./Comprehensive-Exam-Prep/output_tables/", 
               file_name = "pen_col1lo4tbin1", model = "penalized")


penalize_col1lo4qbin1 <- foreach (i = exp(seq(-25, -5, by =  1)) ,
                                 .errorhandling="remove") %dopar% 
                        DT_pen_mcmc(penalty = i, starts_data = col1_low4_5$queen_starts_persec, 
                                    states = states, ant_file = col1_low4_5$data,
                                    hours = hours, 
                                    a = .005, b = .001, c = .005, d = .001,
                                    tau = tau, tau.pen = 0, n_mcmc = n_mcmc, 
                                    delta_t = delta_t, start = start, fig_save = TRUE,
                                    fig_path = path, 
                                    fig_name = "pen_col1lo4qbin1_")

sumtable_model(results = penalize_col1lo4qbin1, compare = penalty, 
               file_path = "./Comprehensive-Exam-Prep/output_tables/", 
               file_name = "pen_col1lo4qbin1", model = "penalized")


penalize_col1lo4ebin1 <- foreach (i = exp(seq(-25, -15, by =  .5)) ,
                                 .errorhandling="remove") %dopar% 
                                DT_pen_mcmc(penalty = i, starts_data = col1_low4_5$entrance_start_persec, 
                                states = states, ant_file = col1_low4_5$data,
                                hours = hours, 
                                a = .005, b = .001, c = .005, d = .001,
                                tau = tau, tau.pen = 0, n_mcmc = n_mcmc, 
                                delta_t = delta_t, start = start, fig_save = TRUE,
                                fig_path = path, 
                                fig_name = "pen_col1lo4ebin1_")

sumtable_model(results = penalize_col1lo4ebin1, compare = exp(seq(-25, -15, by =  .5)), 
               file_path = "./Comprehensive-Exam-Prep/output_tables/", 
               file_name = "pen_col1lo4ebin1", model = "penalized")


penalize_col2lo4qbin1 <- foreach (i = exp(seq(-25, -15, by =  1)) ,
                                  .errorhandling="remove") %dopar% 
  DT_pen_mcmc(penalty = i, starts_data = col2_low4_5$queen_starts_persec, 
              states = states, ant_file = col2_low4_5$data, chamber = "queen",
              hours = hours, 
              a = .005, b = .001, c = .005, d = .001,
              tau = tau, tau.pen = 0, n_mcmc = n_mcmc, 
              delta_t = delta_t, start = start, fig_save = TRUE,
              fig_path = path, 
              fig_name = "pen_col2lo4qbin1_")

sumtable_model(results = penalize_col2lo4qbin1, compare = exp(seq(-25, -15, by =  1)), 
               file_path = "./Comprehensive-Exam-Prep/output_tables/", 
               file_name = "pen_col2lo4qbin1", model = "penalized")

#check exp(-19) over different tuning 
penalize_col2lo4qbin1 <- foreach (i = (seq(.000001, .0001, by =  .00001)) ,
                                  .errorhandling="remove") %dopar% 
  DT_pen_mcmc(tau = matrix( c(i, 0, 
                              0, i), nrow = 2, ncol = 2), penalty = exp(-19), starts_data = col2_low4_5$queen_starts_persec, 
              states = states, ant_file = col2_low4_5$data, chamber = "queen",
              hours = hours, 
              a = .005, b = .001, c = .005, d = .001,
              tau.pen = 0, n_mcmc = n_mcmc, 
              delta_t = delta_t, start = start, fig_save = TRUE,
              fig_path = path, 
              fig_name = "pen_col2lo4qbin1_")

sumtable_model(results = penalize_col2lo4qbin1, compare = seq(.000001, .0001, by =  .00001), 
               file_path = "./Comprehensive-Exam-Prep/output_tables/", 
               file_name = "pen_col2lo4qbin1", model = "penalized")




#Ant Data - Penalized Model with Covariates
penalty <- seq(0, 10, by =  1)
tau <- matrix( c(.0001, 0, 0 , 0,
                  0, .0001, 0, 0,
                  0, 0, .0001, 0, 
                  0, 0, 0, .00001), nrow = 4, ncol = 4)


alpha.beta = c(.005, .001, .005, .001)
start <- list(X = X, lambda = lambda, alpha.beta = alpha.beta)
covariate <- col1_hi4_inout_1$cov


pencov_col1hi4bin1 <- lapply(penalty,
                            FUN = DT_pencov_mcmc,
                            covariate = covariate, title = "Test",
                            starts_data = col1_high4_5$starts_persec, 
                            states = states, ant_file = col1_high4_5$data,
                            hours = hours, start = start,
                            a = .005, b = .001, c = .005, d = .001, 
                            tau = tau, n_mcmc = n_mcmc, delta_t = delta_t,
                            fig_save = TRUE, fig_path = path, 
                            fig_name = "pencov_col1hi4bin1.")

covariate <- col1_lo4_inout_1$cov

pencov_col1lo4tbin1 <- lapply(penalty, FUN = DT_pencov_mcmc,
                              covariate = covariate, title = "Test",
                              starts_data = col1_low4_5$starts_persec, 
                              states = states, ant_file = col1_low4_5$data,
                              hours = hours, start = start,
                              a = .005, b = .001, c = .005, d = .001, 
                              tau = tau, n_mcmc = n_mcmc, delta_t = delta_t,
                              fig_save = TRUE, fig_path = path, 
                              fig_name = "pencov_col1lo4tbin1.")

pencov_col1lo4qbin1 <- lapply(penalty, FUN = DT_pencov_mcmc,
                              covariate = covariate, title = "Test",
                              starts_data = col1_low4_5$queen_starts_persec, 
                              states = states, ant_file = col1_low4_5$data,
                              hours = hours, start = start,
                              a = .005, b = .001, c = .005, d = .001, 
                              tau = tau, n_mcmc = n_mcmc, delta_t = delta_t,
                              fig_save = TRUE, fig_path = path, 
                              fig_name = "pencov_col1lo4qbin1.")

pencov_col1lo4ebin1 <- lapply(penalty, FUN = DT_pencov_mcmc,
                             covariate = covariate, title = "Test",
                              starts_data = col1_low4_5$entrance_start_persec, 
                              states = states, ant_file = col1_low4_5$data,
                              hours = hours, start = start,
                              a = .005, b = .001, c = .005, d = .001, 
                              tau = tau, n_mcmc = n_mcmc, delta_t = delta_t,
                              fig_save = TRUE, fig_path = path, 
                              fig_name = "pencov_col1lo4ebin1.")
