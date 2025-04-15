# calc_intxn_effSize_meta(x_a, s_a, n_a, 
#                         x_b, s_b, n_b, 
#                         x_c, s_c, n_c, 
#                         x_ab, s_ab, n_ab) 
# 
# Function calculates the individual effect size (Hedges g), variance, and
# weight of treatments A, B, and AB. Then, function calculates main effects
# of treatments A, B, and AB. Main effects of treatments A and B differ from
# individual effects as they compare the net effect of a  stressor in presence 
# and absence of a second stressor, while individual effects compare the 
# response in the presence of a since stressor alone.
#
# Function inputs:
#
#   - x_a     = treatment mean of driver A
#   - s_a     = treatment standard deviation of driver A
#   - n_a     = treatment sample size of driver A
#   - x_b     = treatment mean of driver B
#   - s_b     = treatment standard deviation of driver B
#   - n_b     = treatment sample size of driver B
#   - x_c     = treatment mean of control
#   - s_c     = treatment standard deviation of control
#   - n_c     = treatment sample size of control
#   - x_ab    = treatment mean of driver A + B
#   - s_ab     = treatment standard deviation of driver A + B
#   - n_ab    = treatment sample size of driver A + B
# 
# Function returns:
#
# A list containing the individual, main, and interactive effect
# sizes along with effect size variances and weights
calc_intxn_effSize_meta <- function(x_a, s_a, n_a,
                                    x_b, s_b, n_b,
                                    x_c, s_c, n_c,
                                    x_ab, s_ab, n_ab) {
  
  # Degrees of freedom
  m = (n_c + n_a + n_b + n_ab - 4)
  
  ####################
  # Individual effects
  ####################
  
  # Individual standard deviation
  s_ind_a = sqrt(((n_c - 1)*(s_c^2) + (n_a - 1)*(s_a)^2) / (n_a + n_c - 2))
  s_ind_b = sqrt(((n_c - 1)*(s_c^2) + (n_b - 1)*(s_b)^2) / (n_b + n_c - 2))
  s_ind_ab = sqrt(((n_c - 1)*(s_c^2) + (n_ab - 1)*(s_ab)^2) / (n_ab + n_c - 2))
  
  # Individual correction term
  j_ind_a = 1 - ( 3 / (4 * (n_a + n_c - 2) - 1))
  j_ind_b = 1 - ( 3 / (4 * (n_b + n_c - 2) - 1))
  j_ind_ab = 1 - ( 3 / (4 * (n_ab + n_c - 2) - 1))
  
  # Individual effect size
  g_a = (x_a - x_c) / s_ind_a * j_ind_a
  g_b = (x_b - x_c) / s_ind_b * j_ind_b
  g_ab = (x_ab - x_c) / s_ind_ab * j_ind_ab
  
  # Individual effect size variance (v)
  v_a = ((n_a + n_c) / (n_a * n_c)) + ((g_a)^2 / (2 * (n_a + n_c)))
  v_b = ((n_b + n_c) / (n_b * n_c)) + ((g_b)^2 / (2 * (n_b + n_c)))
  v_ab = ((n_ab + n_c) / (n_ab * n_c)) + ((g_ab)^2 / (2 * (n_ab + n_c)))
  
  # Individual effect size weight (w)
  w_a = 1 / v_a
  w_b = 1 / v_b
  w_ab = 1 / v_ab
  
  ####################
  # Interaction effects
  #################### 

  # Pooled standard deviation
  s_int = sqrt( 
    ((n_c - 1)*(s_c^2) + 
       (n_a - 1)*(s_a^2) +
       (n_b - 1)*(s_b^2) +
       (n_ab - 1)*(s_ab^2)) / m)
  
  # Correction term for small sample size
  j_m = 1 - (3 / (4 * m - 1))
  
  # Main effects of treatment A and B
  dA = ((x_a + x_ab) - (x_b + x_c)) / (2 * s_int) * j_m
  dB = ((x_b + x_ab) - (x_a + x_c)) / (2 * s_int) * j_m
  
  # Interaction effect size
  dAB = ((x_ab - x_b) - (x_a - x_c)) / (2 * s_int) * j_m
  
  # Variance of the main and interaction effect sizes
  v_a_main = 0.25 * (
    (1 / n_c) + (1 / n_a) + (1 / n_b) + (1 / n_ab) + 
      ((dA^2) / (2 * (n_c + n_a + n_b + n_ab))))
  
  v_b_main = 0.25 * (
    (1 / n_c) + (1 / n_a) + (1 / n_b) + (1 / n_ab) + 
      ((dB^2) / (2 * (n_c + n_a + n_b + n_ab))))
  
  v_ab_int = 0.25 * (
    (1 / n_c) + (1 / n_a) + (1 / n_b) + (1 / n_ab) + 
      ((dAB^2) / (2 * (n_c + n_a + n_b + n_ab))))
  
  # Weights of the main and interaction effect sizes
  w_a_main = 1 / v_a_main
  w_b_main = 1 / v_b_main
  w_ab_int = 1/ v_ab_int
  
  return(
    data.frame(
      g_a, s_ind_a, j_ind_a, v_a, w_a,
      g_b, s_ind_b, j_ind_b, v_b, w_b,
      g_ab, s_ind_ab, j_ind_ab, v_ab, w_ab,
      
      s_int, j_m,
      
      dA, v_a_main, w_a_main,
      dB, v_b_main, w_b_main,
      dAB, v_ab_int, w_ab_int))
  
}

# Test fxn
# calc_effectSz_intxn(x_a = 2, s_a = 2, n_a = 9,
#                     x_b = 4, s_b = 2, n_b = 9,
#                     x_c = 1, s_c = 1, n_c = 9,
#                     x_ab = 8, s_ab = 4, n_ab = 4)



