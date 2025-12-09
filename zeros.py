########################################################
# Main Repository of Functions for Zeros of MF Project #
########################################################

# Fall 2025

Cprec = 1000
qprec = 1000
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
E6 = Eq(6,qprec)
E8 = Eq(8,qprec)
E10 = Eq(10,qprec)
E14 = Eq(14,qprec) # should be E2, but E2 is problematic
Ers = {0:E0, 4:E4, 6:E6, 8:E8, 10:E10, 14:E14}

def Eratio(k, qprec):
  # careful, because qprec should be considerably larger than l so qprec >> k/12
  rset = [0,4,6,8,10,14]
  r = [u for u in rset if mod(k-u,12)==0][0]
  l = Integer((k-r)/12)
  return Eq(k,qprec)/(Delta**l*Ers[r])

def fratio(f,k):
  # careful, because qprec should be considerably larger than l so qprec >> k/12
  rset = [0,4,6,8,10,14]
  r = [u for u in rset if mod(k-u,12)==0][0]
  l = Integer((k-r)/12)
  return f/(Delta**l*Ers[r])

##
## Functions to get Faber Polynomials
##

j = j_invariant_qexp(qprec)

P = PolynomialRing(QQ, 'J')
J = P.gen()

def Faber_polynomial_from_ratio(fratio):
    v=fratio.valuation()
    av= fratio[v]
    if v == 0:
        return av
    v*=-1 # v is a negative power 
    return av*J**v + Faber_polynomial_from_ratio(fratio - av*j**v)
##
## Functions to find zeros of MF
##

## inverse of j
def j_inverse_wiki(j_val):
  # j-value is complex number
  if j_val == 0:
      return e**(2*pi*I/3)
  x = polygen(CC)
  roots = (x**2-x+(1728/(4*j_val))).roots(CC)
  a = roots[0][0]
  x = I * (hypergeometric([1/6, 5/6], [1], 1-a) / hypergeometric([1/6, 5/6], [1], a))
  return x

# test find z such that j(z) = 732.545683918438*I
  
## move to fundamental domain
def translate_to_fundamental_domain(z):
  """ Translate a point in the complex plane to the fundamental domain of the upper half plane. """
  # REAL PART
  if not (z.real() >= -0.5 and z.real() < 0.5):  # if R{z} in (-1/2, 1/2]
    # translate to (-1/2, 1/2]
    a = z.real()
    return translate_to_fundamental_domain(z - floor(a + 0.5))
  # NORM
  if z.norm() < 1 or (z.norm()==1 and z.real()>0):
    return translate_to_fundamental_domain(-1/z)
  return z

  
## nice j inverse
def j_inverse(j_val):
  z = j_inverse_wiki(j_val)
  z_fund = translate_to_fundamental_domain(z)
  return CC(z_fund)


##
## Nice plots
##

def fundamental_domain(H, xmin=-0.6, xmax=0.6, ymin=0.7):
  # fundamental domain for SL_2(ZZ) action on upper half plane
  # up to imaginary part H
  return arc((0,0), 1, figsize=[20,10],sector=(pi/3,pi/2), aspect_ratio=1, xmin=xmin, xmax=xmax, ymin=ymin, ymax=H, ticks=[[],[]])+arc((0,0), 1, sector=(pi/2,2*pi/3), linestyle='dashed')+line([(0.5, sqrt(3)/2), (0.5, H)])+line([(-0.5, sqrt(3)/2), (-0.5, H)], linestyle='dashed')



def find_roots(fratio):
  j_zeroes = Faber_polynomial_from_ratio(fratio).roots(CC) # George Added .roots

  j_roots = [x for (x,_) in j_zeroes]
  z_roots = []

  for x in j_roots:
    z_root = j_inverse_wiki(x)
    z_roots.append(z_root)
  return fundamental_domain(2)+points(z_roots,color='red')
