########################################################
# Main Repository of Functions for Zeros of MF Project #
########################################################

# Fall 2025

Cprec = 100
qprec = 100
CC = ComplexField(Cprec)

##
## Functions to get Modular Forms
##

# Get Eisenstein Series

def E(k):
  assert is_even(k), "Even weight needed."
  return EisensteinForms(1,k).gens()[0]

def Eq(k, qprec):
  return E(k).qexp(qprec)

Delta = delta_qexp(qprec)

# get weight 0 ratio

E0 = 1
E4 = Eq(4,qprec)
E6 = Eq(4,qprec)
E8 = Eq(4,qprec)
E10 = Eq(4,qprec)
E14 = Eq(4,qprec) # should be E2, but E2 is problematic
Ers = {0:E0, 4:E4, 6:E6, 8:E8, 10:E10, 14:E14}

def Eratio(k, qprec):
  # careful, because qprec should be considerably larger than l so qprec >> k/12
  rset = [0,4,6,8,10,14]
  r = [u for u in rset if mod(k-u,12)==0][0]
  l = (k-r)/12
  return Eq(k,qprec)/(Delta^l*Ers[r])

##
## Functions to get Faber Polynomials
##

j = j_invariant_qexp(qprec)

def Faber_polynomial_from_ratio(fratio):
  # produces a polynomial in j

##
## Functions to find zeros of MF
##


##
## Nice plots
##
